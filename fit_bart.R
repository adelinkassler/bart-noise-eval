#' This script uses preprocessed datasets and executes the actual fit part of
#' BART/DART runs


# Setup -------------------------------------------------------------------

.version.params <- list(
  name = "default",
  notes = NULL,
  started = Sys.time(),
  finished = NULL,
  exit_status = 1,
  sdy = TRUE,
  pscores = TRUE,
  weight_cols = c('n.patients.y4', 'n.patients.y3'),
  stretch = 1,
  bart = list(
    nskip = 200,
    keepevery = 1
  ),
  dart = list(
    nskip = 2000,
    keepevery = 5
  )
)

.arguments <- commandArgs(TRUE)
for (arg in .arguments) {
  for (i in seq_along(args)) {
    param <- sub('=.*$', '', arg[i])
    value <- sub('^.*=', '', arg[i])
    if (param %in% names(.version.params)) {
      if (grepl("^[0-9]+$", value)) value <- as.numeric(value)
      if (grepl("^(true)|(false)$", value)) value <- as.logical(value)
      .version.params[[param]] <- value
      message(paste("Setting", param, "to", value))
    }
  }
}

suppressPackageStartupMessages({
  # QoL
  library(spsUtils)
  library(progress)
  library(testthat)
  library(glue)
  
  # Computation
  library(caret)
  library(BART)
  
  # library(tidyverse) Tidyverse doesn't load on AWS, individually loading
  # packages for compatibility with cloud
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(readr)
  library(tibble)
})

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

strings_to_factors <- function(x) {
  mutate(x, across(where(is.character), as.factor))
}

# Analysis helper funs ----------------------------------------------------

get_rhs_names <- function(x) {
  setdiff(names(x), c("Y_outcome", 'id.practice', 'dataset.num', 'Y_avg.post', 
                      'Y_avg.pre', 'Y_avg_diff.post_pre'))
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


# Main fit and pred --------------------------------------------------------

run_fit_and_pred <- function(path) {
  dataset_name <- basename(path) %>% sub('\\..*$', '', .)
  ext <- sub('^.*\\.', '', path)
  if (ext == 'rds') {
    data <- readRDS(path)
  } else if (ext == 'csv') {
    data <- read_csv(path) %>% strings_to_factors()
  } else {
    stop(glue("Can't handle filetype .{ext} in {path}"))
  }
  
  if (!.version.params$pscores) data$pscore <- NULL
  if (!.version.params$sdy) data$sdY <- NULL
  
  data$Y_outcome <- data$Y_avg_diff.post_pre

  weights <- data %>% 
    select(one_of(.version.params$weight_cols)) %>% 
    rowSums()
  weights <- weights / sum(weights)
  
  data.rhs <- as_rhs(data)
  
  bart.fit <- quietly(gbart)(
    data.rhs,
    data[['Y_outcome']],
    sparse = F,
    nskip = .version.params$stretch * .version.params$bart$nskip,
    keepevery = .version.params$stretch * .version.params$bart$keepevery,
    w = weights
  )$result
  dart.fit <- quietly(gbart)(
    data.rhs,
    data[['Y_outcome']],
    sparse = T,
    nskip = .version.params$stretch * .version.params$dart$nskip,
    keepevery = .version.params$stretch * .version.params$dart$keepevery,
    w = weights
  )$result
  
  if (!dir.exists(glue("Save/{.version.params$name}"))) dir.create(glue("Save/{.version.params$name}"))
  
  fit_save_dir <- glue("Save/{.version.params$name}/fits/")
  if (!dir.exists(fit_save_dir)) dir.create(fit_save_dir)
  saveRDS(bart.fit, glue("{fit_save_dir}/fit_{dataset_name}_BART.rds"))
  saveRDS(dart.fit, glue("{fit_save_dir}/fit_{dataset_name}_DART.rds"))
  
  pdep_post.bart <- pdep(bart.fit, data)
  pdep_post.dart <- pdep(dart.fit, data)
  pdep_post_save_dir <- glue("Save/{.version.params$name}/pdep_post/")
  if (!dir.exists(pdep_post_save_dir)) dir.create(pdep_post_save_dir)
  saveRDS(pdep_post.bart, glue("{pdep_post_save_dir}/pdep_{dataset_name}_BART.rds"))
  saveRDS(pdep_post.dart, glue("{pdep_post_save_dir}/pdep_{dataset_name}_DART.rds"))
  
  # bart.ests <- get_all_estimates(pdep_post.bart, data)
  # dart.ests <- get_all_estimates(pdep_post.dart, data)
  # saveRDS(bart.ests, glue("Data/estimates/{.version.params$name}/ests_{dataset_name}_BART.rds"))
  # saveRDS(dart.ests, glue("Data/estimates/{.version.params$name}/ests_{dataset_name}_DART.rds"))
  # write_csv(bart.ests, glue("Data/analyzed/{.version.params$name}/{dataset_name}_BART.csv"))
  # write_csv(dart.ests, glue("Data/analyzed/{.version.params$name}/{dataset_name}_DART.csv"))
  
  
  return(NULL)
}


# Run main as a loop --------------------------------------------------------

# input_data_dir <- glue("~/Documents/Consulting/BCBS/Data/model_ready/{.version.params$input_dir}")
# if (!dir.exists(input_data_dir)) stop("Input directory does not exist")

output_ests_dir <- glue("~/Documents/Consulting/BCBS/Data/estimates/{.version.params$name}")
if (!dir.exists(output_ests_dir)) dir.create(output_ests_dir)

on.exit({
  .version.params$finished <- Sys.time()
  write_yaml(.version.params, file = file.path(output_ests_dir, "PARAMETERS.yml"))
})

# if (.version.params$loop) {
#   paths <- dir(input_data_dir)
# } else {
#   paths <- grep("\\.(rds)|(csv)$", .arguments, value = TRUE)
# }
paths <- character(0)
for (arg in .arguments) {
  if (dir.exists(arg)) {
      paths <- c(paths, dir(arg, full.names = TRUE))
  } else if (file.exists(arg)) {
      paths <- c(paths, arg)
  }
}
paths <- sample(paths[grepl('(\\.rds)|(\\.csv)', paths)])

pb <- progress_bar$new(
  format = ":current/:total in :elapsed [:bar] eta: :eta\n", 
  total = length(paths),
  show_after = 0
)
invisible(pb$tick(0))

message(glue("*****Fitting BART/DART models to {length(paths)} datasets.*****"))
iwalk(paths, function(path, i) {
  # elapsed <- my_timer(run_fit_and_pred(path))
  # message(glue("*****({i}/{length(paths)}) Fit BART/DART to {path} in {format(elapsed)}.*****"))
  quiet(run_fit_and_pred(path))
  pb$tick()
})

.version.params$exit_status <- 0
