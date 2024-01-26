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

source("utils_eval.R")


# Load fixed datasets -----------------------------------------------------

estimand_truths <- read_csv("Data/ACIC_estimand_truths.csv")
n.patients.combined <- readRDS("Save/n_patients_combined.rds")
acic_winners_evals <- readRDS("Save/acic_winners_performance.rds")
  
# curia_evals <- readRDS("Save/curia_performance.rds")
curia_estimates <- readRDS("Save/all_curia_estimates_1x_10x.rds")

dgp_levels <- full_join(
  estimand_truths %>% 
    select(dataset.num:`Idiosyncrasy of Impacts`) %>% 
    distinct(),
  acic_winners.df %>% 
    distinct(DGP_name, dataset.num)
)

dgp_vars <- c("Confounding Strength", "Confounding Source", 
              "Impact Heterogeneity", "Idiosyncrasy of Impacts")


# Convenience funs --------------------------------------------------------

weighted.se.mean <- function(x, w, na.rm = T){
  # Credit: https://www.alexstephenson.me/post/2022-04-02-weighted-variance-in-r/
  ## Remove NAs 
  if (na.rm) {
    i <- !is.na(x)
    w <- w[i]
    x <- x[i]
  }
  
  ## Calculate effective N and correction factor
  n_eff <- (sum(w))^2/(sum(w^2))
  correction = n_eff/(n_eff-1)
  
  ## Get weighted variance 
  print(paste(length(x), length(w)))
  numerator = sum(w*(x-weighted.mean(x,w))^2)
  denominator = sum(w)
  
  ## get weighted standard error of the mean 
  se_x = sqrt((correction * (numerator/denominator))/n_eff)
  return(se_x)
}

symmetric.setdiff <- function(...) {
  argnames <- as.character(enexprs(...))
  sets <- list(...)
  result <- map(seq_along(sets), \(i) {
    setdiff(sets[[i]], reduce(sets[-i], intersect))
  })
  names(result) <- argnames
  return(result)
}

reverse.symmetric.setdiff <- function(...) {
  argnames <- as.character(enexprs(...))
  sets <- list(...)
  result <- map(seq_along(sets), \(i) {
    setdiff(reduce(sets[-i], union), sets[[i]])
  })
  names(result) <- argnames
  return(result)
}

xor.names <- function(...) {
  argnames <- as.character(enexprs(...))
  sets <- map(list(...), names)
  result <- map(seq_along(sets), \(i) {
    setdiff(sets[[i]], reduce(sets[-i], intersect))
  })
  names(result) <- argnames
  return(result)
}

# Preprocess, combine BART run(s) -----------------------------------------

get_bart_ests <- function(dir_path, pattern = '^.*$') {
  bart_satt_ests.l <- dir(glue(dir_path), full.names = TRUE) %>%
    keep(\(x) grepl(pattern, x)) %>% 
    #keep(~ !file.info(.x)$isdir) %>% 
    #keep(~ file.info(.x)$size > 10) %>% 
    #grep('[0-9]{4}b?(_ydiff_w)?_ps', ., value = TRUE) %>% 
    map(function(x) {
      df <- readRDS(x)
      df$path <- x
      df$level <- as.character(df$level)
      return(df)
    }) 
  bart_ests <- bart_satt_ests.l %>%
    list_rbind()
  return(bart_ests)
}

add_cols_to_bart_ests <- function(bart_ests) {
  bart_ests %>% 
    mutate(
      noise = case_when(
        grepl('curia_10x', path) ~ TRUE,
        grepl('acic_base', path) ~ FALSE,
        .default = NA
      ),
      method = case_when(
        grepl('BART', path) ~ 'BART',
        grepl('DART', path) ~ 'DART'
      ),
      dataset.num = sub('^.*([0-9]{4}).*$', '\\1', basename(path))
    ) %>% 
    select(-path, -se) #score_submission has rename(se=rmse)
}

bart_ests_refactor <- get_bart_ests("Data/estimates/default/", pattern = "/ests") %>% 
  mutate(version = 'refactor')
bart_ests <- add_cols_to_bart_ests(bart_ests_refactor)

# Preprocess, combine curia -----------------------------------------------

curia_ests <- curia_estimates %>% 
  rename(satt = model_satt, se = model_se) %>% 
  mutate(method = 'curia', version = 'fixed') %>% 
  select(-n_columns, -SATT, -any_of(dgp_vars)) %>% 
  mutate(lower.90 = satt - qnorm(.95) * se,
         upper.90 = satt + qnorm(.95) * se) %>% 
  select(-se) #score_submission has rename(se=rmse)

#' The curia estimates are now in the same form as the bart estimates, which we
#' can confirm with:
expect_equal(length(unlist(xor.names(curia_ests, bart_ests))), 0)


# Evals for BART and Curia ------------------------------------------------

truths <- estimand_truths %>%
  left_join(distinct(acic_winners.df, DGP_name, dataset.num)) %>% 
  select(dataset.num, DGP_name, variable, level, year, id.practice, true = SATT)

bart_curia_ests <- rbind(bart_ests, curia_ests)
.tmp <- 

bart_curia_evals.l <- bart_curia_ests %>% 
  group_split(version, noise, method) %>% 
  map(score_submission, truths=truths)# %>% 
  # list_rbind()
bart_curia_evals.df <- bart_curia_evals.l %>% 
  map("eval") %>% 
  list_rbind()

.tmp <- bart_curia_ests %>% 
  nest_by(version, noise, method) %>% 
  mutate(score = map(data, \(d) score_submission(d, truths=truths)))
  

# > names(bart_curia_evals.l[[1]])
# [1] "overall" "dgp_level" "eval" "exemplars" "sensspec"

# Scratch -----------------------------------------------------------------

# Why can't I rbind?
View(map(bart_curia_evals.l, "dgp_level")[[1]])




# score_submission(test_results, truths)


# make_evals <- function(ests_df) {
#   ests_df %>% 
#     # Merge on the truth
#     left_join(estimand_truths %>% 
#                 filter(is.na(year)), 
#               by = c("dataset.num", "variable", "level", "year", "id.practice"),
#               suffix = c(".ests", "")) %>% 
#     select(!ends_with(".ests")) %>% 
#     rename(true = SATT) %>% 
#     left_join(n.patients.combined, by = c("dataset.num", "variable", "level", "year", "id.practice")) %>% 
#     # filter(dataset.num %in% dataset_nums) %>% 
#     mutate(bias = true - satt,
#            abs_bias = abs(bias),
#            cover = between(true, lower.90, upper.90),
#            width = upper.90 - lower.90)
# }
# 
# bart_curia_ests <- rbind(bart_ests, curia_ests)
# bart_curia_evals <- make_evals(bart_curia_ests)


# Combine w/ acic, summarize ----------------------------------------------

acic_winners_evals_ <- acic_winners_evals %>% 
  left_join(n.patients.combined) %>% 
  left_join(dgp_levels) %>% 
  mutate(version = "fixed") %>% 
  filter(!(method %in% c('t1_rBARTimpT', 'ofpsbart', 't1_H-TMLE')))

all_evals <- rbind(bart_curia_evals, acic_winners_evals_)

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

# Make datasets for plotting ----------------------------------------------

plotting_data <- all_evals_summ %>% 
  # Create sensible variables for visualization
  mutate(
    # So we can treat all subgroups with the same variable
    variable_level = paste0(variable, 
                            ifelse(is.na(level), '', paste0('_', level)),
                            ifelse(is.na(year), '', paste0('_y', year))),
    # Sensible method labelling
    method_label = ifelse(grepl("(curia)|(BART)|(DART)", method), method, "other")#,
    #variant = ifelse(noise, '10x', ifelse(outcome == 'ydiff', 'ydiff', 'baseline'))
  )

method_colors <- c(other = 'lightgray', curia = 'green4', BART = 'red2', DART = 'blue3')

# Make line plots ---------------------------------------------------------

map(c("rmse", "bias", "abs_bias", "cover", "width"), function(.x) {
  ggplot(mapping = aes(x = n.practices, y = !!sym(.x), group = method, color = method_label)) + 
    geom_line(data = filter(plotting_data, method_label == 'other')) +
    geom_line(data = filter(plotting_data, method_label != 'other', noise == FALSE)) +
    geom_line(data = filter(plotting_data, method_label != 'other', noise == TRUE), linetype = 2) +
    scale_color_manual(values = method_colors) +
    facet_wrap(~ `Confounding Strength`) +
    {if(.x %in% c("rmse")) scale_y_log10() else NULL} + 
    labs(title = .x)
  ggsave(file = paste0("Plots/default/plot_", .x, "_by_subgroup_size.png"))
})


# Scratchpad --------------------------------------------------------------

truths <- estimand_truths %>%
  left_join(distinct(acic_winners.df, DGP_name, dataset.num)) %>% 
  select(dataset.num, DGP_name, variable, level, year, id.practice, true = SATT)

test_results <- {bart_curia_ests %>% 
  group_split(version, noise, method)}[[1]]

score_submission(test_results, truths)


