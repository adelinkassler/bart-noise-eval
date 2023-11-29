# A script to use on remote nodes to fit BART and DART.
# Adelin Kassler - Nov 2024

# Notes: the original version of the script was meant to be run on data that had
# been formatted into wide format, then had noise added. As of 11/29 it was
# retooled to accept raw data with noise added in ACIC format.


# Setup -------------------------------------------------------------------

library(tidyverse)
library(BART)
library(caret)

if (interactive()) {
  dataset_num <- '0001'
  noise_name <- NULL
} else {
  dataset_num <- commandArgs(TRUE)[1]
  noise_name <- commandArgs(TRUE)[2] %>% {if (is.na(.)) NULL else .}
}

# TODO: if/else to handle base/noise selection from args
dataset_name <- dataset_num

# Paths to load data in ACIC format (maybe +noise). %s = dataset name, i.e. num
# (e.g. "0001"), plus possible noise code (e.g. norm100). Paths may be
# difference for noised data set.
if (is.null(noise_name)) {
  path_p <- "Data/track2_20220404/practice/acic_practice_%s.csv"
  path_py <- "Data/track2_20220404/practice_year/acic_practice_year_%s.csv"
} else {
  path_p <- "Data/noise_gen/practice/acic_practice_%s.csv"
  path_py <- "Data/noise_gen/practice_year/acic_practice_year_%s.csv"
}

# Paths to write analysis output to. %s = BART/DART, %s = dataset name
path_analysis <- "Data/analyzed/estimates_%s_%s.csv"
path_pdep <- "Data/analyzed/pdep_draws_%s_%s.csv"


# Load data ---------------------------------------------------------------

data_py.raw <- read.csv(sprintf(path_py, dataset_name), stringsAsFactors = TRUE)
data_p.raw <- read.csv(sprintf(path_p, dataset_name), stringsAsFactors = TRUE)


# Format data -------------------------------------------------------------

data_py.wide <- data_py.raw %>%
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

data.main <- inner_join(data_py.wide, data_p.raw, by = "id.practice")
class(data.main) <- "data.frame"

rm(data_py.raw, data_p.raw, data_py.wide)


# Fit BART/DART to data+noise ---------------------------------------------

x_names <- setdiff(names(data.main), c("Y_avg.post", 'id.practice'))

# TODO: we need longer burn-in for bart, and especially dart
# 100->200 for bart
# for dart, needs iteration
# skip, don't take more samples
bart.fit <- gbart(data.main[x_names], data.main[['Y_avg.post']], sparse = F, nskip = 200)
dart.fit <- gbart(data.main[x_names], data.main[['Y_avg.post']], sparse = T, nskip = 2000, ndpost = 1000, keepevery = 5)


# Partial dependence calculations -----------------------------------------

pdep <- function(.model, .data) {
  newdata.z1 <- newdata.z0 <- .data %>% 
    select(!c(Y_avg.post, id.practice)) %>%
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


# Run performance functions -----------------------------------------------

# Run the performance functions for BART
pdep_draws.bart <- pdep(bart.fit, data.main)
bart.ests <- get_all_estimates(pdep_draws.bart, data.main)

# Run the performance functions for DART
pdep_draws.dart <- pdep(dart.fit, data.main)
dart.ests <- get_all_estimates(pdep_draws.dart, data.main)


# Output posterior draws and estimates ------------------------------------

write_csv(pdep_draws.bart, sprintf(path_pdep, 'BART', dataset_name))
write_csv(pdep_draws.dart, sprintf(path_pdep, 'DART', dataset_name))

write_csv(bart.ests, sprintf(path_analysis, 'BART', dataset_name))
write_csv(dart.ests, sprintf(path_analysis, 'DART', dataset_name))

# Performance evals in another script/notebook, combined with curia evals &c.
