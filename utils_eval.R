score_submission <- function(results, truths=NULL, dir=compdir, do_exemplars=TRUE) {
    datasets <- sort(unique(results$dataset.num))
    prac_est <- min(is.na(results$id.practice))==0
    if (is.null(truths)) {
        truths <- glue('{dir}/_all_SATT.RDS') %>% readRDS() %>%
            select(dataset.num, DGP_name, variable, level, year, id.practice, true = SATT)
    }

    #Filter truth down to only rows we expect
    #That way we can error if any of the expected rows aren't there, isntead of just inner joining
    #result-only estimates (e.g. for C practices) are fine, we just ignore
    truth <- truths %>% filter(dataset.num %in% datasets)
    if (!prac_est) {
        truth <- filter(truth, is.na(id.practice))
    }

    eval <- truth %>%
        inner_join(results, by=c('dataset.num', 'variable', 'level', 'year', 'id.practice')) %>%
        mutate(bias = satt - true,
               abs_bias = abs(bias),
               rmse = bias^2,
               cover = true >= `lower.90` & true <= `upper.90`,
               width = `upper.90` - `lower.90`)

    if (nrow(eval)!=nrow(truth)) warning('Missing expected estimates')

    dgp_level <- eval %>%
        filter(is.na(id.practice)) %>%
        group_by(DGP_name, variable, level, year) %>%
        summarize(sd_bias = sd(bias, na.rm=TRUE),
                  sd_abs_bias = sd(abs_bias, na.rm=TRUE),
                  sd_width = sd(width, na.rm=TRUE),
                  across(c(bias, abs_bias, rmse, cover, width), mean, na.rm=TRUE, .names = "{.col}_mean"),
                  n_datasets = length(unique(dataset.num))) %>%
        ungroup %>%
        mutate(rmse_mean = sqrt(rmse_mean),
               bias_l = bias_mean - qnorm(.95)*sd_bias/sqrt(n_datasets),
               bias_u = bias_mean + qnorm(.95)*sd_bias/sqrt(n_datasets),
               abs_bias_l = abs_bias_mean - qnorm(.95)*sd_abs_bias/sqrt(n_datasets),
               abs_bias_u = abs_bias_mean + qnorm(.95)*sd_abs_bias/sqrt(n_datasets),
               width_l = width_mean - qnorm(.95)*sd_width/sqrt(n_datasets),
               width_u = width_mean + qnorm(.95)*sd_width/sqrt(n_datasets),
               cover_l = cover_mean - qnorm(.95)*sqrt(cover_mean*(1-cover_mean)/n_datasets),
               cover_u = cover_mean + qnorm(.95)*sqrt(cover_mean*(1-cover_mean)/n_datasets)) %>%
        select(DGP_name, variable, level, year,
               bias=bias_mean, abs_bias = abs_bias_mean,
               rmse = rmse_mean, cover=cover_mean, width=width_mean,
               everything()) %>%
        bind_rows(eval %>%
                      filter(!is.na(level)) %>%
                      score_avg_rollup(TRUE) %>%
                      mutate(variable='Subgroups'))

    if (prac_est) {
        dgp_level <- dgp_level %>%
            bind_rows(eval %>%
                          filter(!is.na(id.practice)) %>%
                          score_avg_rollup(TRUE) %>%
                          mutate(variable='Practices'))
    }

    overall <- eval %>%
        filter(is.na(id.practice)) %>%
        group_by(variable, level, year) %>%
        summarize(across(c(bias, abs_bias, rmse, cover, width), mean, na.rm=TRUE)) %>%
        ungroup %>%
        mutate(rmse = sqrt(rmse)) %>%
        select(variable, level, year, bias, abs_bias, rmse, cover, width) %>%
        bind_rows(eval %>%
                      filter(!is.na(level)) %>%
                      score_avg_rollup(FALSE) %>%
                      mutate(variable='Subgroups'))

    if (prac_est) {
        overall <- overall %>%
            bind_rows(eval %>%
                          filter(!is.na(id.practice)) %>%
                          score_avg_rollup(FALSE) %>%
                          mutate(variable='Practices'))
    }

    truth_diffs <- truth %>%
        inner_join(truths %>% filter(variable=='Overall' & is.na(year)) %>% select(dataset.num, true_overall = true),
                   by='dataset.num') %>%
        mutate(true = true - true_overall) %>%
        select(-true_overall)

    est_diffs <- results %>%
        inner_join(results %>% filter(variable=='Overall' & is.na(year)) %>% select(dataset.num, satt_overall = satt),
                   by='dataset.num') %>%
        mutate(satt = satt - satt_overall,
               lower.90 = lower.90 - satt_overall,
               upper.90 = upper.90 - satt_overall) %>%
        select(-satt_overall) %>%
        inner_join(truth_diffs, by=c('dataset.num', 'variable', 'level', 'year', 'id.practice')) %>%
        mutate(rmse_diff = (satt - true)^2,
               abse_diff = abs(satt) - abs(true),
               excl_0 = lower.90 > 0 | upper.90 < 0,
               excl_5 = lower.90 > 5 | upper.90 < -5,
               excl_10 = lower.90 > 10 | upper.90 < -10)

    overall_diffs <- est_diffs %>%
        filter(is.na(id.practice)) %>%
        group_by(variable, level, year) %>%
        summarize(rmse_diff = sqrt(mean(rmse_diff)),
                  abse_diff = mean(abse_diff)) %>%
        ungroup %>%
        bind_rows(est_diffs %>%
                      filter(!is.na(level)) %>%
                      score_diff_rollup(FALSE) %>%
                      mutate(variable='Subgroups'))
    if (prac_est) {
        overall_diffs <- overall_diffs %>%
            bind_rows(est_diffs %>%
                          filter(!is.na(id.practice)) %>%
                          score_diff_rollup(FALSE) %>%
                          mutate(variable='Practices'))
    }

    dgp_diffs <- est_diffs %>%
        filter(is.na(id.practice)) %>%
        group_by(DGP_name, variable, level, year) %>%
        summarize(rmse_diff = sqrt(mean(rmse_diff)),
                  abse_diff = mean(abse_diff)) %>%
        ungroup %>%
        select(DGP_name, variable, level, year, rmse_diff, abse_diff) %>%
        bind_rows(est_diffs %>%
                      filter(!is.na(level)) %>%
                      score_diff_rollup(TRUE) %>%
                      mutate(variable='Subgroups'))

    if (prac_est) {
        dgp_diffs <- dgp_diffs %>%
            bind_rows(est_diffs %>%
                          filter(!is.na(id.practice)) %>%
                          score_diff_rollup(TRUE) %>%
                          mutate(variable='Practices'))
    }

    overall <- inner_join(overall, overall_diffs, by=c('variable', 'level', 'year'))
    dgp_level <- inner_join(dgp_level, dgp_diffs, by=c('DGP_name', 'variable', 'level', 'year'))

    if (do_exemplars) {
        exemplars <- est_diffs %>%
            filter(is.na(variable) | variable!='Overall') %>%
            mutate(type = ifelse(is.na(id.practice), 'Subgroups','Practices'),
                   power_5 = ifelse(abs(true)>5, excl_0 & sign(satt) == sign(true), NA),
                   power_10 = ifelse(abs(true)>10, excl_0 & sign(satt) == sign(true), NA),
                   types_0 = ifelse(excl_0, sign(satt) == sign(true), NA),
                   types_5 = ifelse(excl_5, sign(satt-5) == sign(true-5), NA),
                   types_10 = ifelse(excl_10, sign(satt-10) == sign(true-10), NA)) %>%
            expand_grid(dumb = 1:3) %>%
            mutate(DGP_name = case_when(dumb==1 ~ 'Overall',
                                        dumb==2 & str_detect(DGP_name,'none') ~ 'RCT',
                                        dumb==2 & str_detect(DGP_name,'r2-high_hetero-high') & !str_detect(DGP_name,'none') ~ 'high_het',
                                        dumb==3 ~ DGP_name,
                                        TRUE ~ 'no')) %>%
            filter(DGP_name!='no') %>%
            group_by(type, DGP_name) %>%
            summarize(across(matches('types|power'), mean, na.rm=TRUE)) %>%
            ungroup

        sensspec <- eval %>%
            mutate(excl_0 = lower.90 > 0 | upper.90 < 0,
                   excl_5 = lower.90 > 5 | upper.90 < -5,
                   excl_10 = lower.90 > 10 | upper.90 < -10,
                   power_5 = ifelse(abs(true)>5, excl_0 & sign(satt) == sign(true), NA),
                   power_10 = ifelse(abs(true)>10, excl_0 & sign(satt) == sign(true), NA),
                   types_0 = ifelse(excl_0, sign(satt) != sign(true), NA),
                   types_5 = ifelse(excl_5, sign(satt-5) != sign(true-5), NA),
                   types_10 = ifelse(excl_10, sign(satt-10) != sign(true-10), NA)) %>%
            expand_grid(dumb=1:4,
                        stupid=1:2) %>%
            mutate(variable = case_when(stupid==1 & is.na(id.practice) ~ variable,
                                        stupid==2 & is.na(id.practice) & variable!='Overall' ~ 'Subgroups',
                                        stupid==2 & !is.na(id.practice) ~ 'Practices',
                                        TRUE ~ 'no'),
                   level = case_when(variable %in% c('Subgroups','Practices') ~ NA_character_, TRUE ~ level),
                   DGP_name = case_when(dumb==1 ~ 'Overall',
                                        dumb==2 & str_detect(DGP_name,'r2-high_hetero-high') & str_detect(DGP_name,'rand') ~ 'var_high_het',
                                        dumb==2 & str_detect(DGP_name,'r2-high_hetero-high') & str_detect(DGP_name,'pretrend') ~ 'pre_high_het',
                                        dumb==2 & !str_detect(DGP_name,'r2-high_hetero-high') ~ 'all_low_het',
                                        dumb==3 & str_detect(DGP_name,'r2-high_hetero-high') & !str_detect(DGP_name,'none') ~ 'all_high_het',
                                        dumb==3 & str_detect(DGP_name,'r2-low_hetero-low') & !str_detect(DGP_name,'none') ~ 'all_lowest_het',
                                        dumb==4 ~ DGP_name,
                                        TRUE ~ 'no')) %>%
            filter(variable!='no' & DGP_name!='no') %>%
            group_by(DGP_name, variable, level, year) %>%
            summarize(across(matches('types|power'), ~ sum(!is.na(.x)), .names='n_{.col}'),
                      across(matches('types|power'), mean, na.rm=TRUE)) %>%
            ungroup

        return(list(overall=overall,
                    dgp_level=dgp_level,
                    eval=eval %>% rename(se=rmse),
                    exemplars=exemplars,
                    sensspec=sensspec))
    } else {
        return(list(overall=overall,
                    dgp_level=dgp_level,
                    eval=eval %>% rename(se=rmse)))
    }
}

score_avg_rollup <- function(subresults, by_DGP) {
    dataset_rollup <- subresults %>%
        group_by(dataset.num, DGP_name) %>%
        summarize(across(c(bias, abs_bias, rmse, cover, width), mean, na.rm=TRUE)) %>%
        ungroup
    if (by_DGP) {
        dataset_rollup %>%
            group_by(DGP_name) %>%
            summarize(across(c(bias, abs_bias, rmse, cover, width), mean)) %>%
            ungroup %>%
            mutate(rmse = sqrt(rmse)) %>%
            return
    } else {
        dataset_rollup %>%
            summarize(across(c(bias, abs_bias, rmse, cover, width), mean)) %>%
            ungroup %>%
            mutate(rmse = sqrt(rmse)) %>%
            return
    }
}

score_diff_rollup <- function(subresults, by_DGP) {
    dataset_rollup <- subresults %>%
        group_by(dataset.num, DGP_name) %>%
        summarize(across(c(rmse_diff, abse_diff), mean, na.rm=TRUE)) %>%
        ungroup
    if (by_DGP) {
        dataset_rollup %>%
            group_by(DGP_name) %>%
            summarize(across(c(rmse_diff, abse_diff), mean)) %>%
            ungroup %>%
            mutate(rmse_diff = sqrt(rmse_diff)) %>%
            return
    } else {
        dataset_rollup %>%
            summarize(across(c(rmse_diff, abse_diff), mean)) %>%
            ungroup %>%
            mutate(rmse_diff = sqrt(rmse_diff)) %>%
            return
    }
}

gather_results <- function(subid) {
    subid1 <- str_extract(subid,".+track[12]")
    subid2 <- str_remove(subid, paste0(subid1,'_'))
    res <- read_csv(glue('evaluation/submissions/{subid1}_overall_{subid2}.csv'), col_types = cols('c','c','c','c','d','d','d'))
    prac_file <- glue('evaluation/submissions/{subid1}_practice-level_{subid2}.csv')
    if (file.exists(prac_file)) {
        res <- res %>%
            bind_rows(read_csv(prac_file, col_types=cols('c','d','d','d','d')))
    } else {
        res$id.practice <- NA_real_
    }

    res <- res %>%
        mutate(dataset.num = str_pad(dataset.num,4,'left','0'),
               year = as.numeric(year),
               variable = case_when(variable=='overall' ~ 'Overall',
                                    TRUE ~ variable),
               level = case_when(level=='0.0' ~ '0',
                                 level=='1.0' ~ '1',
                                 TRUE ~ level)) %>%
        rename(lower.90=lower90,
               upper.90=upper90) %>%
        select(dataset.num, variable, level, year, id.practice, satt, lower.90, upper.90)

    return(res)
}

make_sub_file <- function(subid, name) {
    res <- gather_results(subid)
    x <- score_submission(res)

    bind_rows(x$overall %>%
                  mutate(`Confounding Strength` = 'All',
                         `Confounding Source` = 'All',
                         `Impact Heterogeneity` = 'All',
                         `Idiosyncrasy of Impacts` = 'All') %>%
                  select(`Confounding Strength`, `Confounding Source`,
                         `Impact Heterogeneity`,`Idiosyncrasy of Impacts`,
                         variable, level, year, bias, abs_bias, rmse, cover, width),
              x$dgp_level %>%
                  mutate_knobs %>%
                  mutate(`Confounding Strength` = ifelse(conf=='','None',as.character(conf)),
                         `Confounding Source` = case_when(conf=='' ~ 'None',
                                                          pretrend=='confounding by pretrend' ~ 'Scenario A',
                                                          pretrend=='varying confounding' ~ 'Scenario B'),
                         `Impact Heterogeneity` = str_extract(hetero,'Large|Small'),
                         `Idiosyncrasy of Impacts` = case_when(r2=='mostly X-explained' ~ 'Small',
                                                               r2=='mostly idiosyncratic' ~ 'Large')) %>%
                  select(`Confounding Strength`, `Confounding Source`,
                         `Impact Heterogeneity`,`Idiosyncrasy of Impacts`,
                         variable, level, year, bias, abs_bias, rmse, cover, width)) %>%
        mutate(`Confounding Strength` = factor(`Confounding Strength`, levels=c('All','None','Weak','Strong')),
               `Confounding Source` = factor(`Confounding Source`, levels=c('All','None','Scenario A','Scenario B')),
               `Impact Heterogeneity` = factor(`Impact Heterogeneity`, levels=c('All','Small','Large')),
               `Idiosyncrasy of Impacts` = factor(`Idiosyncrasy of Impacts`, levels=c('All','Small','Large')),
               variable = ifelse(is.na(year),variable,'Yearly')) %>%
        mutate(conf_id=name) %>%
        select(Submission = conf_id, everything()) %>%
        arrange(Submission, `Confounding Strength`, `Confounding Source`,
                `Impact Heterogeneity`,`Idiosyncrasy of Impacts`, variable, level, year)
}
