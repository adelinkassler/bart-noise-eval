#' This script runs BART and DART on ACIC-formatted data (or related files,
#' e.g., from curia). It provides functions to for pre- and post-processing of
#' data for and from BART fits, and a main function and loop at the bottom to do
#' so for the relevant datasets.
#' 
#' Created by Adelin Kassler, Nov 2023
#' 
#' Last updated by Adelin Kassler, Dec 2023

library(caret)
library(glue)
library(testthat)
library(BART)
# library(tidyverse)
library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(tibble)


# Global parameters -------------------------------------------------------

params <- list(
  bart = list(
    sparse = F,
    nskip = 200,
    keepevery = 1
  ),
  dart = list(
    sparse = T,
    nskip = 2000,
    keepevery = 5
  )
)

if (!interactive()) {
  idx_skip <- as.numeric(commandArgs(TRUE)[1])
  idx_offset <- as.numeric(commandArgs(TRUE)[2])
} else {
  idx_skip <- 1
  idx_offset <- 1
}


# Processing curia data ---------------------------------------------------

dedup_cols <- function(x) {
  dup_cols <- duplicated(t(x))
  x[!dup_cols]
}

get_practice_level_vars <- function(x) {
  expect_in(c('id.practice', 'Z', 'year'), names(x))
  
  x %>% 
    group_by(id.practice) %>% 
    summarise(
      n = n(),
      across(!n, n_distinct)
    ) %>% 
    select(-id.practice, -Z) %>% 
    {as.data.frame(t(rbind(
      mean = summarise(., across(everything(), mean)) %>% round(1),
      sd = summarise(., across(everything(), sd)) %>% round(3)
    )))} %>% 
    filter(mean == 1, sd == 0) %>% 
    rownames()
}

select_practice_level <- function(x) {
  expect_in(c('id.practice', 'Z'), names(x))
  practice_level_vars <- get_practice_level_vars(x)
  select(x, id.practice, Z, all_of(practice_level_vars))
}

select_practice_year_level <- function(x) {
  expect_in(c('id.practice', 'Z'), names(x))
  practice_level_vars <- get_practice_level_vars(x)
  select(x, id.practice, Z, !all_of(practice_level_vars))
}



# Format data -------------------------------------------------------------

strings_to_factors <- function(x) {
  mutate(x, across(where(is.character), as.factor))
}

add_sd_y_practice <- function(x, patient_path=NULL) {
  #patient.df <- read_csv(patient_path)
  patient.df <- readRDS("sdY.rds")
  
  left_join(x, patient.df)
}

format_practice_year_wide <- function(practice_year) {
  expect_in(c('id.practice', 'Z', 'year', 'Y', 'n.patients'), 
             names(practice_year))
  
  practice_year %>%
    pivot_longer(!c(id.practice, year, Z)) %>%
    group_by(id.practice, name) %>%
    arrange(id.practice, name, year) %>% 
    mutate(diff1 = value - lag(value),
           diff2 = value - lag(lag(value)),
           diff3 = value - lag(lag(lag(value)))) %>% 
    ungroup() %>%
    pivot_wider(id_cols = c(id.practice, Z),
                names_from = c(name, year),
                values_from = c(value, diff1, diff2, diff3),
                names_glue = "{name}{ifelse(.value=='value','','_diff')}.y{year}{ifelse(.value=='value','',paste0('_y',year-parse_number(.value,na='value')))}") %>% 
    select(!starts_with("post") & 
             !ends_with(c('0', '-1', '-2')) & 
             !starts_with(c("Y_diff.y3", "Y_diff.y4"))) %>% 
    # TODO: outcome variable should be Y_avg.post - Y_avg.pre
    mutate(Y_avg.post = (Y.y3 * n.patients.y3 + Y.y4 * n.patients.y4) /
             (n.patients.y3 + n.patients.y4)) %>%
    select(-Y.y4, -Y.y3)
}

add_pscores <- function(x, method = 'BART') {
  return(x) # Placeholder
  
  data <- as_rhs(x)
  rhs_names <- get_practice_level_vars(x)
  
  # Extensible to other methods
  if (method == 'BART') {
    fit <- gbart(data[rhs_names], data$Z, sparse = F, nskip = params$bart$nskip,
                 keepevery = params$bart$keepevery)
  } else if (method == 'DART') {
    fit <- gbart(data[rhs_names], data$Z, sparse = T, nskip = params$dart$nskip, 
                 keepevery = params$dart$keepevery)
  } else if (method == 'linear') {
    .fmla = as.formula(glue("Z ~ ", glue_collapse(rhs_names, sep = ' + ')))
    fit <- lm(data = data, formula = .fmla)
  }
  
  x %>% 
    mutate(pscore = fit$yhat.train.mean)
}

join_practice_level <- function(practice, practice_year) {
  x <- inner_join(practice, practice_year, by = c('id.practice', 'Z'))
  class(x) <- 'data.frame'
  x
}

get_rhs_names <- function(x) {
  #expect_in(c("Y_avg.post", 'id.practice'), names(x))
  setdiff(names(x), c("Y_avg.post", 'id.practice'))
}

as_rhs <- function(x) {
  if (!(class(x)[1] %in% c('data.frame', 'matrix'))) 
    class(x) <- 'data.frame'
  
  x[get_rhs_names(x)] %>% 
    strings_to_factors()
}


# Performance evals from BART fits ----------------------------------------

pdep <- function(.model, .data) {
  newdata.z1 <- newdata.z0 <- .data %>% 
    select(!c(Y_avg.post, id.practice)) %>%
    as_rhs() %>% 
    bartModelMatrix()
  newdata.z0[,'Z'] <- 0
  newdata.z1[,'Z'] <- 1
  
  pred.z0 <- predict(.model, newdata.z0)
  pred.z1 <- predict(.model, newdata.z1)
  
  pred.diff <- pred.z1 - pred.z0
  return(pred.diff)
  
}

posterior_subgroup_avgs <- function(pred.diff, .data) {
  #pred.diff <- pdep(.model, .data)
  
  # TODO: debug this line
  if (is.null(dim(pred.diff))) return(pred.diff)
  
  # Use total number of patients as weights
  n.patients <- .data %>% 
    # TODO: fix this line
    select(n.patients.y3, n.patients.y4) %>% 
    rowSums()
  
  # Take the weighted average of subgroup
  # Returns a vector of length = # posterior draws
  return(apply(pred.diff, 1, function(.) weighted.mean(., n.patients)))
}

posterior_estimate <- function(pred.diff, .data) {
  pavg <- posterior_subgroup_avgs(pred.diff, .data)
  
  estimate <- mean(pavg)
  # Using a 90% credible interval
  interval <- quantile(pavg, c(.05, .95))
  
  #return(list(estimate = estimate, interval = interval))
  return(data.frame(est_satt = estimate, q05 = interval[[1]], q95 = interval[[2]]))
}

get_all_estimates <- function(pred.diff, .data) {
  # Because .data is masked in dplyr functions
  ..data <- .data
  
  satt.overall <- posterior_estimate(pred.diff, ..data) %>% 
    mutate(variable = "Overall")
  
  # We don't have a satt.yearly because we analyzed a wide dataset with no time
  # groups
  
  satt.practice <- list_rbind(map(unique(..data$id.practice), function(.x) {
    ..data.sub <- ..data %>% mutate(.idx = 1:nrow(..data)) %>% filter(id.practice == .x)
    posterior_estimate(pred.diff[, ..data.sub$.idx], ..data.sub) %>% 
      mutate(id.practice = .x)
  }))
  
  satt.levels <- list_rbind(map(1:5, function(.x)
    list_rbind(map(unique(..data[[paste0("X", .x)]]), function(.y) {
      ..data.sub <- ..data %>% 
        mutate(.idx = 1:nrow(..data)) %>% 
        # filter(variable == paste0("X", .x), level == .y)
        filter(!!sym(paste0("X", .x)) == .y)
      posterior_estimate(pred.diff[, ..data.sub$.idx], ..data.sub) %>% 
        mutate(variable = paste0("X", .x), level = as.character(.y))
    })) 
  ))
  
  return(bind_rows(satt.overall, satt.levels, satt.practice))
  
}


# Main --------------------------------------------------------------------

main.curia_merged <- function(data_path) {
  # merged_data_path <- glue("Data/curia_data/df_merged_raw_practice_level_10x/",
  #                       "merged_dataset_practice_{dataset_num}b.csv")
  merged_data <- quietly(read_csv)(data_path)$result %>% 
    dedup_cols() %>% 
    strings_to_factors() %>% 
    rename(Y = label, n.patients = n_patients_in_practice) %>% 
    select(-mo)
  
  practice <- select_practice_level(merged_data)
  practice_year <- select_practice_year_level(merged_data)
  
  data <- join_practice_level(
    practice_year = practice_year %>% format_practice_year_wide(),
    practice = practice %>% add_pscores() %>% add_sd_y_practice()
  )
  
  # TODO: we need longer burn-in for bart, and especially dart
  # 100->200 for bart
  # for dart, needs iteration
  # skip, don't take more samples
  data.rhs <- as_rhs(data)
  bart.fit <- gbart(data.rhs, data[['Y_avg.post']], sparse = F, 
                    nskip = params$bart$nskip, keepevery = params$bart$keepevery)
  dart.fit <- gbart(data.rhs, data[['Y_avg.post']], sparse = T, 
                    nskip = params$dart$nskip, keepevery = params$dart$keepevery)
  
  pdep_post.bart <- pdep(bart.fit, data)
  pdep_post.dart <- pdep(dart.fit, data)
  
  bart.ests <- get_all_estimates(pdep_post.bart, data)
  dart.ests <- get_all_estimates(pdep_post.dart, data)
  
  dataset_name <- basename(data_path) %>% sub('\\..*$', '', .)
  
  saveRDS(bart.fit, glue("fits/fit_{dataset_name}_BART.rds"))
  saveRDS(dart.fit, glue("fits/fit_{dataset_name}_DART.rds"))
  saveRDS(pdep_post.bart, glue("pdep_post/pdep_{dataset_name}_BART.rds"))
  saveRDS(pdep_post.dart, glue("pdep_post/pdep_{dataset_name}_DART.rds"))
  write_csv(bart.ests, glue("analyzed/{dataset_name}_BART.csv"))
  write_csv(dart.ests, glue("analyzed/{dataset_name}_DART.csv"))
  
  return(NULL)
}

dur <- c()
cap <- list()
idx_select <- seq(from = idx_offset, 
                  to = length(dir("df_merged_raw_practice_level/")),
                  by = idx_skip)
message(glue("Beginning {length(idx_select)} iterations:"))
for (path in dir("df_merged_raw_practice_level/", 
                 full.names = TRUE)[idx_select]) {
  message(glue("(Run #{length(dur)+1}) Working on {path}..."))
  start <- Sys.time()
  try(cap <- append(cap, list(quietly(main.curia_merged)(path))))
  dur <- append(dur, Sys.time() - start)
  message(glue("Completed in {format(dur[length(dur)])}"))
}


# Runtime stats -----------------------------------------------------------

#' Summarize runtimes with key stats and graphs
# Trimming durations avoids outsized impact of, eg, closing lid of laptop during
# a long string of runs
avg_dur <- mean(dur, trim = floor(log10(length(dur))) / length(dur))
dur_quant <- paste(c("Min", '25%', "Med", '75%', "Max"),
                   format(as.numeric(quantile(dur, c(0, .25, .5, .75, 1))), digits = 3), 
                   sep = ' = ', collapse = ', ')
message(glue("Completed {length(dur)} runs\nMean duration/run = {format(avg_dur,",
             "digits = 3)}\nQuantiles: {dur_quant}"))

hist(as.numeric(dur), main = "Histogram of Runtimes", xlab = "Runtime")
abline(v = avg_dur, col = 'red')

#' Captured output:
for (x in cap) {
  cat(x$output)
  if (length(x$warnings)>0) warning(x$warnings)
}

#' Save environment for inspection
# save(file = 'Save/fit_bart_functional.rds')
