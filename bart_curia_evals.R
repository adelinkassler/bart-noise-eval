#' This script creates all the plots and tables we're interested in


# Setup -------------------------------------------------------------------


suppressPackageStartupMessages({
  library(beepr)
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
base_error_handler <- options("error")
loud_error_handler <- function() {
  beep(9)
  # Calling base error handler
  if (!is.null(old_error_handler) && is.function(old_error_handler$error)) {
    old_error_handler$error()
  }
}
options(error = loud_error_handler)

# Load fixed datasets -----------------------------------------------------

focus_datasets <- read_csv("Data/curia_data/dataset_nums.csv", show_col_types = FALSE)
estimand_truths <- read_csv("Data/ACIC_estimand_truths.csv", show_col_types = FALSE)
n.patients.combined <- readRDS("Save/n_patients_combined.rds")
acic_winners_evals <- readRDS("Save/acic_winners_performance.rds")
mpr_bart_evals.l <- readRDS("Save/mathematica_bart_models.rds")
  
# curia_evals <- readRDS("Save/curia_performance.rds")
curia_estimates <- readRDS("Save/all_curia_estimates_1x_10x.rds")

dgp_levels <- full_join(
  estimand_truths %>% 
    select(dataset.num:`Idiosyncrasy of Impacts`) %>% 
    distinct(),
  acic_winners_evals %>% 
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

bart_ests_default <- get_bart_ests("Data/estimates/default/", pattern = "/ests") %>% mutate(version = 'default')
bart_ests_nops <- get_bart_ests("Data/estimates/nops/", pattern = "/ests") %>% mutate(version = 'nops')
bart_ests_nosdy <- get_bart_ests("Data/estimates/nosdy/", pattern = "/ests") %>% mutate(version = 'nosdy')
bart_ests_nops_nosdy <- get_bart_ests("Data/estimates/nops_nosdy/", pattern = "/ests") %>% mutate(version = 'nops_nosdy')

bart_ests_rbind <- bind_rows(bart_ests_default, bart_ests_nops, bart_ests_nosdy, bart_ests_nops_nosdy)
bart_ests <- add_cols_to_bart_ests(bart_ests_rbind)

# Preprocess, combine curia -----------------------------------------------

curia_ests <- curia_estimates %>% 
  rename(satt = model_satt, se = model_se) %>% 
  mutate(method = 'curia', version = 'fixed') %>% 
  select(-n_columns, -SATT, -any_of(dgp_vars)) %>% 
  mutate(lower.90 = satt - qnorm(.95) * se,
         upper.90 = satt + qnorm(.95) * se) %>% 
  select(-se) %>% #score_submission has rename(se=rmse)
  filter(dataset.num %in% focus_datasets$dataset_num)

#' The curia estimates are now in the same form as the bart estimates, which we
#' can confirm with:
expect_equal(length(unlist(xor.names(curia_ests, bart_ests))), 0)


# Evals for BART and Curia ------------------------------------------------

truths <- estimand_truths %>%
  left_join(distinct(acic_winners_evals, DGP_name, dataset.num)) %>% 
  select(dataset.num, DGP_name, variable, level, year, id.practice, true = SATT)

bart_curia_ests <- rbind(bart_ests, curia_ests)

scores.df <- bart_curia_ests %>% 
  nest_by(version, noise, method) %>% 
  mutate(score = list(score_submission(data, truths=truths))) %>% 
  unnest_wider(score)

saveRDS(scores.df, "Save/bart_curia_full_scores.rds")

bart_curia_evals <- scores.df %>% 
  select(version, noise, method, eval) %>% 
  unnest(eval) %>% 
  left_join(n.patients.combined) %>% 
  left_join(dgp_levels)

saveRDS(bart_curia_evals, "Save/bart_curia_evals.rds")

# Combine w/ acic, summarize ----------------------------------------------

all_evals <- acic_winners_evals %>% 
  left_join(n.patients.combined) %>% 
  left_join(dgp_levels) %>% 
  mutate(version = "fixed") %>% 
  filter(!(method %in% c('t1_rBARTimpT', 'ofpsbart', 't1_H-TMLE'))) %>% 
  bind_rows(bart_curia_evals)

saveRDS(all_evals, "Save/all_evals.rds")

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

message("Completed successfully")
beep(8)
