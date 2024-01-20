#' Run all the loops in one convenient master function

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


# Global parameters -------------------------------------------------------

# Parameters to pass to BART::gbart
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

# Parameters that maybe set at commmand line
idx_skip <- 1
idx_offset <- 1
do_base <- TRUE
do_ps <- TRUE
do_sd <- TRUE
do_2l <- TRUE

if (interactive()) {
  do_curia10x <- TRUE
  do_acic600 <- TRUE
} else {
  args <- commandArgs(TRUE)
  do_curia10x <- 'curia10x' %in% args
  do_acic600 <- 'acic600' %in% args
  do_base <- 'base' %in% args
  do_ps <- 'ps' %in% args
  do_sd <- 'sd' %in% args
  do_2l <- '2l' %in% args
  for (i in seq_along(args)) {
    if (args[i] %in% '--skip') {
      idx_skip <- as.numeric(args[i+1])
    } else if (args[i] %in% '--offset') {
      idx_offset <- as.numeric(args[i+1])
    }
  }
}


# Global Data -------------------------------------------------------------

sdY.df <- readRDS("sdY.rds")


# Convenience functions ---------------------------------------------------

my_timer <- function (expr, gcFirst = TRUE) {
  ppt <- function(y) {
    if (!is.na(y[4L])) 
      y[1L] <- y[1L] + y[4L]
    if (!is.na(y[5L])) 
      y[2L] <- y[2L] + y[5L]
    paste(formatC(y[1L:3L]), collapse = " ")
  }
  if (gcFirst) 
    gc(FALSE)
  time <- Sys.time()
  on.exit(message("Timing stopped at: ", ppt(Sys.time() - 
                                               time)))
  expr
  new.time <- Sys.time()
  on.exit()
  structure(new.time - time, class = "difftime")
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

add_sd_y_practice <- function(x, dsnum) {
  #patient.df <- read_csv(patient_path)
  #patient.df <- readRDS("sdY.rds")
  # Unnecessary because patient.df now loaded globally
  
  left_join(x, sdY.df %>% filter(as.numeric(dataset.num) == dsnum))
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
    # mutate(Y_avg.post = (Y.y3 * n.patients.y3 + Y.y4 * n.patients.y4) /
    #          (n.patients.y3 + n.patients.y4)) %>%
    mutate(Y_outcome = Y.y4 - Y.y3) %>% 
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
  # x <- inner_join(practice, practice_year, by = c('id.practice', 'Z'))
  x <- inner_join(practice, practice_year)
  class(x) <- 'data.frame'
  #TODO: figure out what this join is producing 4x repeated entries
  distinct(x)
}

get_rhs_names <- function(x) {
  #expect_in(c("Y_avg.post", 'id.practice'), names(x))
  setdiff(names(x), c("Y_outcome", 'id.practice', 'dataset.num'))
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
    select(!c(Y_outcome, id.practice)) %>%
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

run_fit_and_pred <- function(data, dataset_name, .long) {
  # Run fit and prediction functions, and save the outputs
  data.rhs <- as_rhs(data)
  if (.long) {
    stretch = 2
  } else {
    stretch = 1
  }
  
  weights <- data %>% 
    select(starts_with('n.patients.')) %>% 
    rowSums()
  weights <- weights / sum(weights)
  
  bart.fit <- gbart(
    data.rhs,
    data[['Y_outcome']],
    sparse = F,
    nskip = stretch * params$bart$nskip,
    keepevery = stretch * params$bart$keepevery,
    w = weights
  )
  dart.fit <- gbart(
    data.rhs,
    data[['Y_outcome']],
    sparse = T,
    nskip = stretch * params$dart$nskip,
    keepevery = stretch * params$dart$keepevery,
    w = weights
  )
  
  saveRDS(bart.fit, glue("Save/fits/fit_{dataset_name}_BART.rds"))
  saveRDS(dart.fit, glue("Save/fits/fit_{dataset_name}_DART.rds"))
  
  pdep_post.bart <- pdep(bart.fit, data)
  pdep_post.dart <- pdep(dart.fit, data)
  saveRDS(pdep_post.bart, glue("Save/pdep_post/pdep_{dataset_name}_BART.rds"))
  saveRDS(pdep_post.dart, glue("Save/pdep_post/pdep_{dataset_name}_DART.rds"))
  
  bart.ests <- get_all_estimates(pdep_post.bart, data)
  dart.ests <- get_all_estimates(pdep_post.dart, data)
  write_csv(bart.ests, glue("Data/analyzed/{dataset_name}_BART.csv"))
  write_csv(dart.ests, glue("Data/analyzed/{dataset_name}_DART.csv"))
  
  return(NULL)
}

main_with_output <- function(main_fn, x, ...) {
  # elapsed <- my_timer(try(cap <- append(cap, list(quietly(main_fn)(x, ...)))))
  elapsed <- my_timer(cap <- append(cap, try(list(quietly(main_fn)(x, ...)))))
  args_list <- as.list(match.call()[-1])
  message(glue("\t...{paste(names(args_list), args_list, sep = '=', collapse = ', ')},",
               " took {format(elapsed)}"))
}

# elapsed <- my_timer(try(cap <- append(cap, list(quietly(main_fn)(path, .pscore=F, .sdy=F)))))
# message(glue("\t...base, took {format(elapsed)}"))
# elapsed <- my_timer(try(cap <- append(cap, list(quietly(main_fn)(path, .sdy=F)))))
# message(glue("\t...w/ p-scores, took {format(elapsed)}"))
# elapsed <- my_timer(try(cap <- append(cap, list(quietly(main_fn)(path)))))
# message(glue("\t...w/ p-scores & sdY, took {format(elapsed)}"))
# elapsed <- my_timer(try(cap <- append(cap, list(quietly(main_fn)(path, .long=T)))))
# message(glue("\t...w/ p-scores, sdY, 2x mcmc runtime, took {format(elapsed)}"))


# Curia Main Loop ---------------------------------------------------------

main.curia_merged <- function(data_path, .pscore = TRUE, .sdy = TRUE, .long = FALSE) {
  # Load/clean curia merged data
  merged_data <- quietly(read_csv)(data_path)$result %>% 
    dedup_cols() %>% 
    strings_to_factors() %>% 
    rename(Y = label, n.patients = n_patients_in_practice) %>% 
    select(-mo)
  
  dsnum <- parse_number(basename(data_path))
  
  # Split, process, and remerge for fit-ready data
  practice <- select_practice_level(merged_data)
  if (.pscore) practice <- practice %>% add_pscores()
  if (.sdy) practice <- practice %>% add_sd_y_practice(dsnum)
  
  practice_year <- select_practice_year_level(merged_data) %>% 
    format_practice_year_wide()
  
  data <- join_practice_level(
    practice_year = practice_year,
    practice = practice
  )
  
  dataset_name <- basename(data_path) %>% sub('\\..*$', '', .)
  #saveRDS(data, glue("Save/wide/wide_{dataset_name}.rds"))
  if (.pscore) dataset_name <- paste0(dataset_name, '_ps')
  if (.sdy) dataset_name <- paste0(dataset_name, '_sd')
  if (.long) dataset_name <- paste0(dataset_name, '_2l')
  
  run_fit_and_pred(data, dataset_name, .long)
}


# ACIC Main Loop ----------------------------------------------------------

dataset_nums <- read_csv("Data/curia_data/dataset_nums.csv", 
                         show_col_types = FALSE) %>% 
  pull(dataset_num)

dataset_paths <- tibble(
  practice = dir("Data/track2_20220404/practice/", full.names = TRUE),
  practice_year = dir("Data/track2_20220404/practice_year/", full.names = TRUE),
  number = 1:3400
) #%>% 
  #filter(number %in% as.numeric(dataset_nums))


main.track2 <- function(dsnum, .pscore = TRUE, .sdy = TRUE, .long = FALSE) {
  dsnum <- as.numeric(dsnum)

  practice_year <- dataset_paths$practice_year[dsnum] %>% 
    read_csv() %>% 
    strings_to_factors() %>% 
    format_practice_year_wide()
  
  practice <- dataset_paths$practice[dsnum] %>% 
    read_csv() %>% 
    strings_to_factors() 
  if (.pscore) practice <- practice %>% add_pscores()
  if (.sdy) practice <- practice %>% add_sd_y_practice(dsnum)

  data <- join_practice_level(
    practice_year = practice_year,
    practice = practice
  )
  
  dataset_name <- sprintf("acic_practice_%.04d", dsnum)
  dataset_name <- paste0(dataset_name, "ydiff")
  #saveRDS(data, glue("Save/wide/wide_{dataset_name}.rds"))
  if (.pscore) dataset_name <- paste0(dataset_name, '_ps')
  if (.sdy) dataset_name <- paste0(dataset_name, '_sd')
  if (.long) dataset_name <- paste0(dataset_name, '_2l')
  
  run_fit_and_pred(data, dataset_name, .long)
}


# Run main loops ----------------------------------------------------------

# Curia 10x noisy
if (do_curia10x) {
dur <- c()
cap <- list()
idx_select <- seq(from = idx_offset, 
                  to = length(dir("Data/curia_data/df_merged_raw_practice_level_10x/")),
                  by = idx_skip)
paths <- dir("Data/curia_data/df_merged_raw_practice_level_10x/", 
             full.names = TRUE)[idx_select]
for (x in paths) {
  message(glue("(Run #{length(dur)+1}) Working on {x}..."))
  start <- Sys.time()
  
  if (do_base) main_with_output(main.curia_merged, x, .pscore = F, .sd = F)
  if (do_ps)   main_with_output(main.curia_merged, x, .pscore = T, .sd = F)
  if (do_sd)   main_with_output(main.curia_merged, x, .pscore = T, .sd = T)
  if (do_2l)   main_with_output(main.curia_merged, x, .pscore = T, .sd = T, .long = T)
  
  dur <- append(dur, Sys.time() - start)
  message(glue("\tCompleted in {format(dur[length(dur)])} at {Sys.time()}"))
}
}

# ACIC default data - curia 600
if (do_acic600) {
dur <- c()
cap <- list()
numbers <- unique(dataset_nums)[seq(from = idx_offset, to = length(unique(dataset_nums)), by = idx_skip)]
for (x in numbers) {
  message(glue("(Run #{length(dur)+1}) Working on ACIC {x}..."))
  start <- Sys.time()
  
  if (do_base) main_with_output(main.track2, x, .pscore = F, .sd = F)
  if (do_ps)   main_with_output(main.track2, x, .pscore = T, .sd = F)
  if (do_sd)   main_with_output(main.track2, x, .pscore = T, .sd = T)
  if (do_2l)   main_with_output(main.track2, x, .pscore = T, .sd = T, .long = T)
  
  dur <- append(dur, Sys.time() - start)
  message(glue("\tCompleted in {format(dur[length(dur)])} at {Sys.time()}"))
}
}

