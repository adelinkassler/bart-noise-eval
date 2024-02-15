#' this script is used to post process the posterior samples for the partial
#' dependence function of an arbitrary BART fit, and turn it into the format of
#' estimates for different subgroups that is used for the ACIC scores. From
#' there the results are able to be used for data, visualization and comparison.
#' 
#' Author: Adelin Kassler, 2024/01/22


# Setup -------------------------------------------------------------------

suppressPackageStartupMessages({
  library(progress)
  library(caret)
  library(glue)
  library(testthat)
  library(BART)
  library(tidyverse)
})

focus_datasets <- read_csv("Data/curia_data/dataset_nums.csv", show_col_types = FALSE)

main_args_df <- map(c("default", "nosdy"), 
                    \(x) {focus_datasets %>% 
  uncount(2) %>% 
  mutate(
    i = 1:n(),
    method = rep(c('BART', 'DART'), times = n()/2),
    suffix = case_when(
      p == '1x' ~ "acic_base",
      p == '10x' ~ "curia_10x"
    ),
    data_file = glue("Data/model_ready/default/patient_wide_preproc_{dataset_num}_{suffix}.rds"),
    pdep_file = glue("Save/{x}/pdep_post/pdep_patient_wide_preproc_{dataset_num}_{suffix}_{method}.rds"),
    out_file = glue("Data/estimates/{x}/ests_{dataset_num}_{suffix}_{method}.rds")
  )
}) %>% bind_rows() %>% 
  mutate(i = 1:n())

# Handle command line parameters ------------------------------------------

#' At the moment, we're not set up to handle command line parameters

# .global.params <- list(
#   weight_cols <- c('n.patients.y3', 'n.patients.y4')
#   # version_name = "default",
#   # dataset_names = c("acic_base", "curia_10x")
# )
# 
# if (interactive()) {
#   .global.params$pdep_files = glue("Save/default/pdep_post/pdep_patient_wide_preproc_0001_acic_base_BART.rds")
#   .global.params$data_files = glue("Data/model_ready/default/patient_wide_preproc_0001_acic_base.rds")
# } else {
#   .args <- commandArgs(TRUE)
#   for (i in 1:length(.args)) {
#     if (.args[i] == "--pdep") {
#       if (dir.exists(.args[i+1])) {
#         .global.params$pdep_files <- dir(.args[i+1])
#       } else if (file.exists(.args[i+1])) {
#         .global.params$pdep_files <- .args[i+1]
#       }
#     } 
#     if (.args[i] == "--data") {
#       if (dir.exists(.args[i+1])) {
#         .global.params$data_files <- dir(.args[i+1])
#       } else if (file.exists(.args[i+1])) {
#         .global.params$data_files <- .args[i+1]
#       }
#     }
#   }
# }
# 
# .global.params$n_files <- length(.global.params$pdep_files)
# expect_equal(length(.global.params$pdep_files), 
#              length(.global.params$data_files))

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


# Subgroup estimates ------------------------------------------------------

calc_weighted_est <- function(pdep_post, weights) {
  weights <- weights / mean(weights)
  est_draws <- apply(pdep_post, 1, \(x) weighted.mean(x, weights))
  est <- mean(est_draws)
  se <- sd(est_draws) / length(est_draws)
  interval.90 <- quantile(est_draws, c(.05, .95))
  return(data.frame(
    satt = est,
    se = se,
    lower.90 = interval.90[[1]],
    upper.90 = interval.90[[2]]
  ))
}

calc_all_subset_ests <- function(pdep_post, data) {
  weights <- data %>% 
    select(n.patients.y3, n.patients.y4) %>% 
    rowSums()
  
  # Overall estimates
  overall_est <- calc_weighted_est(pdep_post, weights) %>% 
    mutate(variable = "Overall", level = NA)
  
  # Individual practice ests
  practice_ests <- map(1:nrow(data), function(i) {
    pdep_row <- pdep_post[, i, drop = FALSE]
    calc_weighted_est(pdep_row, 1) %>% 
      mutate(id.practice = data$id.practice[i])
  }) %>% 
    list_rbind() %>% 
    mutate(variable = NA, level = NA)
  
  # Variable-level ests
  var_level_ests <- map(1:5, function(var_num) {
    # Map over each categorical X variable by number
    map(unique(data[[paste0("X", var_num)]]), function(var_level) {
      # Calc subset weighted est for each level of variable X{i}
      subset_select <- data[[paste0("X", var_num)]] == var_level
      calc_weighted_est(pdep_post[, subset_select, drop = FALSE], 
                        weights[subset_select]) %>% 
        # Level must be character for later rbind
        mutate(level = as.character(var_level))
    }) %>% 
      list_rbind() %>% 
      mutate(variable = paste0("X", var_num))
  })
  
  bind_rows(overall_est, var_level_ests, practice_ests) %>% 
    mutate(year = NA)
}


# Run functions on all arguments ------------------------------------------

n_runs <- nrow(main_args_df)

# Set up progress bar
pb <- progress_bar$new(
  format = ":current/:total in :elapsed [:bar] eta: :eta\n", 
  total = n_runs,
  show_after = 0
)
invisible(pb$tick(0))

# Loop over arguments
pwalk(
  main_args_df,
  function(pdep_file, data_file, out_file, dataset_num, p, method, i, ...) {
    # message(glue("\n***** #{i}/{n_runs} at {format(Sys.time())} *****\n",
    #              "Calculating estimates for dataset={dataset_num}, cols={p}, method={method}",
    #              "\n**********\n"))
    pdep <- readRDS(pdep_file)
    data <- readRDS(data_file)
    ests <- calc_all_subset_ests(pdep, data)
    saveRDS(ests, out_file)
    pb$tick()
    }
)

# Save results to file ----------------------------------------------------


