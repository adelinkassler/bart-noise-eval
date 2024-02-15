#' This script creates all the plots and tables we're interested in


# Setup -------------------------------------------------------------------


suppressPackageStartupMessages({
  library(caret)
  library(glue)
  library(testthat)
  library(BART)
  
  # library(tidyverse) Tidyverse doesn't load on AWS, individually loading
  # packages for compatibility with cloud
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(readr)
  library(tibble)
})

version_name <- "comp_ps_sdy"

# Load fixed datasets -----------------------------------------------------

estimand_truths <- read_csv("Data/ACIC_estimand_truths.csv")
n.patients.combined <- readRDS("Save/n_patients_combined.rds")
acic_winners_evals <- readRDS("Save/acic_winners_performance.rds")
all_evals <- readRDS("Save/all_evals.rds")
focus_datasets <- read_csv("Data/curia_data/dataset_nums.csv", show_col_types = FALSE)
mpr_bart_evals.l <- readRDS("Save/mathematica_bart_models.rds")

# all_evals <- rbind(all_evals, mpr_bart_evals)

dgp_levels <- full_join(
  estimand_truths %>% 
    select(dataset.num:`Idiosyncrasy of Impacts`) %>% 
    distinct(),
  acic_winners_evals %>% 
    distinct(DGP_name, dataset.num)
) %>% mutate(
  `Confounding Strength` = factor(`Confounding Strength`,
                                  levels = c("None", "Weak", "Strong"),
                                  labels = c("No confounding (RCT)",
                                             "Weak confounding",
                                             "Strong confounding")
  ))

focus_dgps <- dgp_levels %>% 
  filter(dataset.num %in% focus_datasets$dataset_num) %>% 
  select(-dataset.num) %>% 
  distinct()

dgp_vars <- c("Confounding Strength", "Confounding Source", 
              "Impact Heterogeneity", "Idiosyncrasy of Impacts")

rm(estimand_truths, acic_winners_evals)

#Turn version into method for better plotting
all_evals <- all_evals %>% 
  mutate(
    method_base = method,
    method = ifelse(method %in% c("BART", "DART"), paste0(method, "_", version), method)
  )

# Convenience functions ---------------------------------------------------

xor.names <- function(...) {
  argnames <- as.character(enexprs(...))
  sets <- map(list(...), names)
  result <- map(seq_along(sets), \(i) {
    setdiff(sets[[i]], reduce(sets[-i], intersect))
  })
  names(result) <- argnames
  return(result)
}


# Summarize ---------------------------------------------------------------

# summarize_evals <- function(evals) {
#   evals %>% 
#     mutate(rmse = bias^2) %>% 
#     group_by(across(!c(satt:upper.90, true:width, rmse))) %>% 
#     summarise(.groups = 'drop',
#               across(c(satt:upper.90, true, bias:width, rmse), 
#                      ~ weighted.mean(.x, n.patients)), 
#               n.patients = sum(n.patients),
#               n.practices = n()) %>% 
#     mutate(rmse = sqrt(rmse)) #%>% 
#     # filter(is.na(year)) %>% 
#     # select(-year)
# }
# 
# all_evals_summ <- summarize_evals(all_evals)

all_evals_summ <- all_evals %>%
  filter(dataset.num %in% focus_datasets$dataset_num) %>%
  group_by(DGP_name, variable, level, year, version, noise, method,
           `Confounding Strength`) %>%
  summarise(
    .groups = 'drop',
    across(c(bias, abs_bias, se, cover, width, n.patients),
           \(x) mean(x, na.rm = TRUE))
  ) %>% 
  mutate(rmse = sqrt(se)) %>% 
  filter(is.na(year)) %>% 
  filter(method != "t1_BC-SGN") %>% 
  filter(!(is.na(variable) & is.na(level) & is.na(year)))

# all_evals_summ <- bind_rows(all_evals_summ, mpr_bart_evals)

# all_evals_summ <- all_evals_summ %>% 
#   mutate(method = ifelse(method %in% c("BART", "DART"), paste0(method, "_", "version"), method))


# Add MPR evals -----------------------------------------------------------

mpr_bart_evals <- mpr_bart_evals.l %>% 
  list_rbind(names_to = "method") %>% 
  mutate(method = paste0("mpr_", method)) %>% 
  mutate(year = NA, version = "fixed", noise = FALSE) %>% 
  left_join(focus_dgps[c("Confounding Strength", "DGP_name")])


# Make datasets for plotting ----------------------------------------------

plotting_data <- all_evals_summ %>% 
  filter(version %in% c("fixed", "default")) %>% 
  # Create sensible variables for visualization
  mutate(
    # So we can treat all subgroups with the same variable
    variable_level = paste0(variable, 
                            ifelse(is.na(level), '', paste0('_', level)),
                            ifelse(is.na(year), '', paste0('_y', year))),
    # Sensible method labelling
    method_label = case_when(
      grepl("^(curia|BART|DART)$", method) ~ method,
      grepl("^(BART_|DART_)", method) ~ method,
      grepl("BART", method) ~ "other_BART",
      .default = "other"
    ),
    `Confounding Strength` = factor(
      `Confounding Strength`,
      levels = c("None", "Weak", "Strong"),
      labels = c("No confounding (RCT)", 
                 "Weak confounding", 
                 "Strong confounding")
    )
  )


# data.plot <- all_evals %>% 
#   # Summarize across datasets
#   summarise(
#     .by = c(variable, level, year, method, `Confounding Strength`, noise),
#     subgroup_size = mean(n.patients, na.rm = TRUE),
#     rmse = caret::RMSE(satt, true),
#     across(c(bias, abs_bias, cover, width), ~weighted.mean(., n.patients))
#   ) %>% 
#   filter(is.na(year)) %>% 
#   # Remove unfinished ACIC submissions
#   filter(!(method %in% c('t1_rBARTimpT', 'ofpsbart', 't1_H-TMLE'))) %>%
#   # Create sensible variables for visualization
#   mutate(
#     # So we can treat all subgroups with the same variable
#     variable_level = paste0(variable, 
#                             ifelse(is.na(level), '', paste0('_', level)),
#                             ifelse(is.na(year), '', paste0('_y', year))),
#     # Ensures methods of interest are at the top of the draw order
#     # method = factor(method, levels = rev(c('curia', 'BART', 'DART', 
#     #                                        setdiff(unique(method), 
#     #                                                c('curia', 'BART', 'DART'))))),
#     # method = factor(method, levels = rev(c('curia', 'BART', 'DART', '_ps_BART', '_ps_DART', '_pd_sd_BART', '_ps_sd_DART', '_ps_sd_2l_BART', '_ps_sd_2l_DART', setdiff(unique(method), c('curia', 'BART', 'DART', '_ps_BART', '_ps_DART', '_pd_sd_BART', '_ps_sd_DART', '_ps_sd_2l_BART', '_ps_sd_2l_DART'))))),
#     # Sensible method labelling
#     method_label = ifelse(grepl("^(curia|BART|DART)$", method), method, "other")#,
#   )


method_colors <- c(
  other = 'lightgray',
  curia = 'green4',
  BART_default = 'red',
  DART_default = 'blue3',
  Oracle = 'purple',
  other_BART = 'orange',
  t2_aipwboot = 'turquoise'
)

# Make line plots ---------------------------------------------------------

map(c("rmse", "bias", "abs_bias", "cover", "width"), function(.x) {
  ggplot(mapping = aes(x = n.patients, y = !!sym(.x), group = method, color = method_label)) + 
    geom_line(data = filter(plotting_data, method_label == 'other'), linewidth = 0.2) +
    geom_line(data = filter(plotting_data, method_label != 'other', is.na(noise)), linewidth = 0.2) +
    geom_line(data = filter(plotting_data, method_label != 'other', noise == FALSE)) +
    geom_line(data = filter(plotting_data, method_label != 'other', noise == TRUE), linetype = 2) +
    scale_color_manual(values = method_colors) +
    facet_wrap(~ `Confounding Strength`) +
    {if(.x %in% c("rmse", "abs_bias", "width")) scale_y_log10() else NULL} + 
    labs(title = .x)
  ggsave(
    file = file.path("Plots",
                     version_name,
                     paste0("plot_", .x, "_by_subgroup_size.png")),
    width = 10,
    height = 3.5
  )
})


# Scatter plots -----------------------------------------------------------


# Mare's plot_perf code
perf <- plotting_data %>%
  filter(!(method %in% c('t1_rBARTimpT', 'ofpsbart', 't1_H-TMLE'))) %>%
  filter(variable=='Overall' & is.na(year))

# perf <- rbind(perf, rep(NA, ncol(perf)))
# perf$method[nrow(perf)] <- 'Oracle'
# perf$rmse[nrow(perf)] <- 0
# perf$width[nrow(perf)] <- 0
# perf$cover[nrow(perf)] <- .9

focus_dgps <- dgp_levels %>% 
  filter(dataset.num %in% focus_datasets$dataset_num) %>% 
  distinct(DGP_name, .keep_all = TRUE)

oracle_row <- tibble(
  method = "Oracle",
  rmse = 0,
  width = 0,
  cover = 0.9,
  variable = "Overall",
  bias = 0,
  abs_bias = 0,
  DGP_name = focus_dgps$DGP_name,
  `Confounding Strength` = focus_dgps$`Confounding Strength`
)
perf <- bind_rows(perf, oracle_row) %>% 
  mutate(Method = case_match(method, 
                             't2_aipwboot' ~ 't2_aipwboot',
                             'curia' ~ 'curia', 
                             'Oracle' ~ 'Oracle',
                             'BART' ~ 'BART',
                             'DART' ~ 'DART',
                             .default = 'Other'),
         noise = ifelse(is.na(noise), FALSE, noise))
# mutate(Method = case_match(
#   method,
#   'curia' ~ 'Curia',
#   'Oracle' ~ 'Oracle',
#   'BART' ~ 'BART',
#   'DART' ~ 'DART',
#   .default = "Other"
# ))

# perf <- perf %>%
#   mutate(Method = ifelse(method == 'curia', 'Curia', 
#                          ifelse(method == 'Oracle', 'Oracle', 'Other')))

# method_colors <- c(Curia = "#F8766D", Oracle = "#00BFC4", Other = "darkgrey")
# method_colors <- c(Curia = "#F8766D", Oracle = "#00BFC4", Other = "darkgrey", 
#                    BART = '#4356FF', DART = '#AB54F9')


#meta <- readRDS('real_submissions.RDS')

curia_rmse <- perf %>%
  filter(Method=='Curia') %>%
  pull(rmse)

perf %>%
  summarise(percent_better = mean(rmse < curia_rmse, na.rm=T))

perf %>%
  ggplot(aes(x=rmse, y=cover, color=Method, shape=noise)) + 
  geom_hline(yintercept = 0.9, color='darkgrey') +
  geom_point() + 
  xlab('RMSE, overall SATT') +
  ylab('Coverage of overall SATT 90% interval') +
  facet_wrap(~ `Confounding Strength`) +
  scale_color_manual(values=method_colors)
ggsave(file.path('Plots', version_name, 'rmse_cover.png'), height=3, width=9)

perf %>%
  ggplot(aes(x=width, y=cover, color=Method, shape=noise)) + 
  geom_hline(yintercept = 0.9, color='darkgrey')+
  geom_point() + 
  xlab('Width of overall SATT 90% interval') +
  ylab('Coverage of overall SATT 90% interval') +
  facet_wrap(~ `Confounding Strength`) +
  scale_color_manual(values=method_colors)
ggsave(file.path("Plots", version_name, 'width_cover.png'), height=3, width=9)

# BART method comparison --------------------------------------------------

all_evals %>% 
  filter(variable == "Overall", is.na(year)) %>% 
  group_by(DGP_name, noise, method) %>%
  summarize(across(satt:width, mean)) %>%
  mutate(rmse = sqrt(se)) %>%
  filter(method %in% c("BART", "DART", "curia"))

bart_vers_evals <- all_evals %>% 
  filter(variable == "Overall",
         is.na(year),
         is.na(noise) | noise == FALSE,
         dataset.num %in% focus_datasets$dataset_num | is.na(dataset.num)) %>%
  group_by(DGP_name, method, `Confounding Strength`) %>%
  summarize(across(bias:width, mean), .groups = 'drop') %>%
  mutate(rmse = sqrt(se)) %>%
  filter(grepl("BART|DART|curia", method)) %>% 
  arrange(DGP_name, rmse) %>% 
  mutate(`Confounding Strength` = factor(`Confounding Strength`,
                                         levels = c("None", "Weak", "Strong"),
                                         labels = c("No confounding (RCT)",
                                                    "Weak confounding",
                                                    "Strong confounding")
  )) %>% 
  bind_rows(mpr_bart_evals %>% filter(variable == "Overall"))

ggplot(data = bart_vers_evals,
       aes(
         x = method,
         y = rmse,
         group = DGP_name,
         fill = `Confounding Strength`, 
         color = `Confounding Strength`
       )) +
  geom_line() + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) + 
  scale_x_discrete(limits = rev(arrange(summarise(bart_vers_evals, .by = method, r = mean(rmse)), r)$method)) +
  labs(title = "Comparison of RMSE for BART versions")
  # geom_bar(stat = "identity", position = "dodge") + 
  # coord_flip()
ggsave(
  file.path("Plots", version_name, "BART_methods_rmse.png"),
  width = 6,
  height = 3
)

ggplot(data = bart_vers_evals,
       aes(
         x = method,
         y = cover,
         group = DGP_name,
         fill = `Confounding Strength`, 
         color = `Confounding Strength`
       )) +
  geom_hline(yintercept = 0.9, linetype = 2) +
  geom_line() + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) + 
  scale_x_discrete(limits = rev(arrange(summarise(bart_vers_evals, .by = method, r = mean(rmse)), r)$method)) +
  labs(title = "Comparison of 90% Interval Coverage for BART versions")
# geom_bar(stat = "identity", position = "dodge") + 
# coord_flip()
ggsave(
  file.path("Plots", version_name, "BART_methods_cover.png"),
  width = 6,
  height = 3
)

BART_method_colors = c(
  t2_xvalBART = "t2_xvalBART",
  t2_xvalBARTps = "t2_xvalBART",
  t1_flexBART_1 = "t1_flexBART",
  t1_flexBART_2 = "t1_flexBART",
  t2_flexBART_4 = "t2_flexBART",
  t2_flexBART_2 = "t2_flexBART",
  t2_flexBART_3 = "t2_flexBART",
  t2_flexBART_1 = "t2_flexBART",
  t2_FisherBART = "t2_FisherBART",
  DART_nops = "gDART",
  DART_default = "gDART",
  DART_nosdy = "gDART",
  DART_nops_nosdy = "gDART",
  BART_default = "gBART",
  BART_nosdy = "gBART",
  BART_nops = "gBART",
  BART_nops_nosdy = "gBART",
  mpr_nopsbart = "mpr",
  mpr_psbart = "mpr",
  mpr_wbcf = "mpr"
)

BART_method_shapes = c(
  t2_xvalBART = "no pscores",
  t2_xvalBARTps = "default",
  t1_flexBART_1 = "default",
  t1_flexBART_2 = "default",
  t2_flexBART_4 = "default",
  t2_flexBART_2 = "default",
  t2_flexBART_3 = "default",
  t2_flexBART_1 = "default",
  t2_FisherBART = "default",
  DART_nops = "no pscores",
  DART_default = "default",
  DART_nosdy = "no sdy",
  DART_nops_nosdy = "no pscores/sdy",
  BART_default = "default",
  BART_nosdy = "no sdy",
  BART_nops = "no pscores",
  BART_nops_nosdy = "no pscores/sdy",
  mpr_nopsbart = "no pscores",
  mpr_psbart = "default",
  mpr_wbcf = "wbcf"
)

ggplot(data = bart_vers_evals,
       aes(
         x = rmse,
         y = cover,
         group = DGP_name,
         color = BART_method_colors[method],
         shape = BART_method_shapes[method]
       )) +
  geom_hline(yintercept = 0.9, linetype = 2, color = 'gray') +
  geom_point() +
  facet_wrap(~`Confounding Strength`, scales = 'free') +
  labs(title = "Comparison of BART versions", color = "Method", shape = "Version")
ggsave("Plots/BART_vers_cover_rmse.png")
  
# Scratchpad --------------------------------------------------------------
stop("Reached end of script")

ggplot(mapping = aes(group = method, x = subgroup_size, y = rmse, color = method_label)) + 
  geom_line(data = filter(data.plot, method_label == 'other')) +
  geom_line(data = filter(data.plot, method_label != 'other', noise == FALSE)) +
  geom_line(data = filter(data.plot, method_label != 'other', noise == TRUE), linetype = 2) +
  # geom_line(aes(group = method, color = method_label)) +
  # scale_color_manual(values = c(curia = 'darkblue', BART = 'darkred', DART = 'cyan', other = 'lightgray')) +
  scale_color_manual(values = method_colors) +
  #scale_linetype_identity(labels = c('10x', 'ydiff', 'base'), breaks = c(2, 3, 1), guide = guide_legend()) +
  facet_wrap(~ `Confounding Strength`) +
  scale_y_log10()

ggplot(mapping = aes(group = method, x = subgroup_size, y = bias, color = method_label)) + 
  geom_line(data = filter(data.plot, method_label == 'other')) +
  geom_line(data = filter(data.plot, method_label != 'other', noise == FALSE)) +
  geom_line(data = filter(data.plot, method_label != 'other', noise == TRUE), linetype = 2) +
  # geom_line(aes(group = method, color = method_label)) +
  # scale_color_manual(values = c(curia = 'darkblue', BART = 'darkred', DART = 'cyan', other = 'lightgray')) +
  scale_color_manual(values = method_colors) +
  #scale_linetype_identity(labels = c('10x', 'ydiff', 'base'), breaks = c(2, 3, 1), guide = guide_legend()) +
  facet_wrap(~ `Confounding Strength`) +
  scale_y_log10()

ggplot(mapping = aes(group = method, x = subgroup_size, y = abs_bias, color = method_label)) + 
  geom_line(data = filter(data.plot, method_label == 'other')) +
  geom_line(data = filter(data.plot, method_label != 'other', noise == FALSE)) +
  geom_line(data = filter(data.plot, method_label != 'other', noise == TRUE), linetype = 2) +
  # geom_line(aes(group = method, color = method_label)) +
  # scale_color_manual(values = c(curia = 'darkblue', BART = 'darkred', DART = 'cyan', other = 'lightgray')) +
  scale_color_manual(values = method_colors) +
  #scale_linetype_identity(labels = c('10x', 'ydiff', 'base'), breaks = c(2, 3, 1), guide = guide_legend()) +
  facet_wrap(~ `Confounding Strength`, ) +
  scale_y_log10()

ggplot(mapping = aes(group = method, x = subgroup_size, y = cover, color = method_label)) + 
  geom_line(data = filter(data.plot, method_label == 'other')) +
  geom_line(data = filter(data.plot, method_label != 'other', noise == FALSE)) +
  geom_line(data = filter(data.plot, method_label != 'other', noise == TRUE), linetype = 2) +
  # geom_line(aes(group = method, color = method_label)) +
  # scale_color_manual(values = c(curia = 'darkblue', BART = 'darkred', DART = 'cyan', other = 'lightgray')) +
  scale_color_manual(values = method_colors) +
  #scale_linetype_identity(labels = c('10x', 'ydiff', 'base'), breaks = c(2, 3, 1), guide = guide_legend()) +
  facet_wrap(~ `Confounding Strength`) +
  scale_y_log10()

ggplot(mapping = aes(group = method, x = subgroup_size, y = width, color = method_label)) + 
  geom_line(data = filter(data.plot, method_label == 'other')) +
  geom_line(data = filter(data.plot, method_label != 'other', noise == FALSE)) +
  geom_line(data = filter(data.plot, method_label != 'other', noise == TRUE), linetype = 2) +
  # geom_line(aes(group = method, color = method_label)) +
  # scale_color_manual(values = c(curia = 'darkblue', BART = 'darkred', DART = 'cyan', other = 'lightgray')) +
  scale_color_manual(values = method_colors) +
  #scale_linetype_identity(labels = c('10x', 'ydiff', 'base'), breaks = c(2, 3, 1), guide = guide_legend()) +
  facet_wrap(~ `Confounding Strength`) +
  scale_y_log10()



map(c("rmse", "bias", "abs_bias", "cover", "width"), function(.x) {
  ggplot(mapping = aes(x = subgroup_size, y = !!sym(.x), group = method, color = method_label)) + 
    geom_line(data = filter(data.plot, method_label == 'other')) +
    geom_line(data = filter(data.plot, method_label != 'other', noise == FALSE)) +
    geom_line(data = filter(data.plot, method_label != 'other', noise == TRUE), linetype = 2) +
    # scale_color_manual(values = c(curia = 'darkblue', BART = 'darkred', DART = 'cyan', other = 'lightgray')) +
    scale_color_manual(values = method_colors) +
    facet_wrap(~ `Confounding Strength`) +
    {if(.x %in% c("rmse")) scale_y_log10() else NULL} + 
    labs(title = .x)
    
  #ggsave(file = paste0("Plots/plot_", .x, "_by_subgroup_size_10x.png"))
})
  

