#' Fit wbcf, analogous to fit_bart.R

# Setup -------------------------------------------------------------------

.version.params <- list(
  name = "bcf_default",
  notes = NULL,
  started = Sys.time(),
  finished = NULL,
  exit_status = 1,
  sdy = TRUE,
  weight_cols = c('n.patients.y4', 'n.patients.y3')
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
  library(testthat)
  library(glue)
  library(stringr)
  library(spsUtil)
  library(progress)
  
  # Computational
  library(caret)
  library(BART)
  library(bcf)
  
  # Tidyverse
  # library(tidyverse) Tidyverse doesn't load on AWS, individually loading
  # packages for compatibility with cloud
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(readr)
  library(tibble)
})

source("acic_bcf.r")
bcf <- bcf::bcf


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

cols_to_remove <- c('id.practice', 'dataset.num', 'Y_avg.post', 
                    'Y_avg.pre', 'Y_avg_diff.post_pre')

get_rhs_names <- function(x) {
  setdiff(names(x), c("Y_outcome", cols_to_remove))
}


# Load datasets -----------------------------------------------------------

focus_datasets <- read_csv("Data/curia_data/dataset_nums.csv", show_col_types = FALSE)


# Modified fns from acic_bcf.r --------------------------------------------

make_bcf_dataset <- function(dataset, pscores, weights=NULL, splityear=FALSE, dir=compdir) {
  # data <- glue('{dir}/track_2_dataset_{dataset}.RDS') %>%
  #   readRDS()
  
  # Modified from original above for new path
  data <- readRDS(glue("Data/t2_joined/acic_ppy_{str_pad(dataset,4,'left','0')}_jd.rds"))
  
  if (splityear) {
    post_data <- data %>%
      filter(post==1) %>%
      #Joined & died are added on these files as time-invariant, so just grab one value
      group_by(id.practice, year, joined, died) %>%
      summarize(Y_post = weighted.mean(Y, n.patients),
                size_post = sum(n.patients),
                across(matches('^V'), ~ weighted.mean(.x, n.patients))) %>%
      ungroup %>%
      rename_all(str_remove,'_avg')
  } else {
    post_data <- data %>%
      filter(post==1) %>%
      group_by(id.practice, joined, died) %>%
      summarize(Y_post = weighted.mean(Y, n.patients),
                size_post = sum(n.patients),
                across(matches('^V'), ~ weighted.mean(.x, n.patients))) %>%
      ungroup %>%
      rename_all(str_remove,'_avg') %>%
      mutate(year = 4)
  }
  
  pre_data <- data %>%
    filter(post==0) %>%
    select(id.practice, year, Z, Y, n.patients, matches('^X')) %>%
    group_by(id.practice) %>%
    mutate(Y_pre = weighted.mean(Y, n.patients),
           size_pre = sum(n.patients),
           across(matches('^V'), ~ weighted.mean(.x, n.patients))) %>%
    ungroup %>%
    pivot_wider(c(id.practice, Z, matches('X'), Y_pre, size_pre), names_from=year, values_from=Y, names_prefix='Y_') %>%
    #In theory we could enter X2/X4 as ordinal, but as participants we can't tell ordinal from categorical
    mutate(Y_trend = Y_1 - Y_2,
           X2A = X2=='A',
           X2B = X2=='B',
           X2C = X2=='C',
           X4A = X4=='A',
           X4B = X4=='B',
           X4C = X4=='C',
           across(matches('X[24][ABC]'), as.numeric))
  
  data <- post_data %>%
    inner_join(pre_data, by='id.practice') %>%
    mutate(w=size_post,
           Y_diff = Y_post - Y_pre) %>%
    inner_join(pscores %>% filter(dataset.num==str_pad(dataset,4,'left','0')) %>% select(-dataset.num), by='id.practice') %>%
    select(id.practice, year, Z, w, Y_diff, Y_post, Y_pre, Y_trend, Y_1, Y_2, size_post, size_pre, pscore, matches('X'), matches('V'), joined, died)
  
  if (!is.null(weights)) {
    data <- data %>%
      inner_join(weights %>% filter(dataset.num==str_pad(dataset,4,'left','0')) %>% select(-dataset.num), by='id.practice') %>%
      mutate(w = w*match_wgt) %>%
      select(-match_wgt)
  }
  
  return(data)
}


# run_bcf (my version) ----------------------------------------------------

run_bcf.data <- function(data, dataset, pscores, ibcf=TRUE, ubcf=FALSE, prac_est=TRUE,
                    n_chains=4, n_bcf_cores=1, nburn=1000, nsim=1000,
                    block_b0_b1=TRUE,
                    ...,
                    dir=compdir,
                    outcome = 'Y_diff',
                    splityear=FALSE,
                    byyear=FALSE,
                    additional_covdata=NULL,
                    weights=NULL) {
  gc()
  #Pretty sure we don't want weights - BCF should be deconfounding on its own
  #data <- make_bcf_dataset(dataset, pscores, weights, splityear, dir=dir)

  #Drop dummied Xs
  #Any reason to include ypre/y_trend AND y1/y2?
  #If I had to pick one pre/trend seems like the more natural transform
  #For OLS it's all collinear, but tree stuff does care I think
  #Worst case including both sets just increases the likelihood of the tree trying to split on what we know a priori is an important var
  xmat <- data %>% 
    # Rename according to mpr naming scheme
    rename(Y_pre = Y_avg.pre, Y_trend = Y_diff.y2_y1) %>% 
    mutate(size_pre = n.patients.y3 + n.patients.y4)
    select(Y_pre, Y_trend, Y_1, Y_2, size_pre, matches('^X'), matches('^V')) %>% 
    select(-X2, -X4) %>% 
    as.matrix
  if (splityear) {
    xmat <- cbind(xmat, data$year)
  }

  if (!is.null(additional_covdata)) {
    new_covars <- setdiff(colnames(additional_covdata), colnames(data))
    if ('dataset' %in% colnames(additional_covdata)) {
      additional_covdata <- additional_covdata[additional_covdata$dataset == dataset]
    } else if ('dataset.num' %in% colnames(additional_covdata)) {
      additional_covdata <- additional_covdata[additional_covdata$dataset.num == str_pad(dataset,4,'left','0')]
    }
    data <- inner_join(data, additional_covdata)
    xmat <- cbind(xmat, data %>% select(all_of(new_covars)) %>% as.matrix)
  }

  #For uBCF, need to feed BCF 1s instead of weights
  #but don't want to change data$w, since we still want weights used to calculate ATTs
  if (ubcf) {
    w_for_bcf <- rep(1, nrow(data))
  } else {
    w_for_bcf <- data$w
  }

  if (n_bcf_cores==1) {
    sink(tempfile(fileext='.log'))
  }
  fit <- bcf(     y           = data[[outcome]],
                  z           = data$Z,
                  x_control   = xmat,
                  x_moderate  = xmat,
                  pihat       = data$pscore,
                  w           = w_for_bcf,
                  verbose     = 0,
                  nburn       = nburn,
                  nsim        = nsim,
                  n_chains    = n_chains,
                  n_cores     = n_bcf_cores,
                  n_threads   = 1,
                  block_b0_b1 = block_b0_b1,
                  random_seed = dataset,
                  include_random_effects = ibcf,
                  save_tree_directory = NULL,
                  log_file   = tempfile(fileext='.log'),
                  simplified_return = TRUE,
                  ...)
  if (n_bcf_cores==1) {
    sink(NULL)
  }

  tau <- lapply(fit$raw_chains, function(x) x$tau) %>% abind::abind(along=3) %>% aperm(c(1,3,2))
  impacts <- bcf_impacts(data, tau, prac_est, splityear, byyear) %>%
    mutate(dataset.num = str_pad(dataset,4,'left','0')) %>%
    select(dataset.num, everything())

  mixing <- get_mixing(fit$raw_chains, data$w)

  ret <- impacts %>%
    mutate(ibcf=ibcf,
           nburn=nburn,
           nsim=nsim,
           n_chains=n_chains) %>%
    bind_cols(mixing)

  if (ibcf) {
    ret$avg_acc_sigv <- mean(lapply(fit$raw_chains, `[`, 'acc_sigv') %>% unlist)
    ret$pct_acc_sigv_lt10 <- mean(lapply(fit$raw_chains, `[`, 'acc_sigv') %>% unlist < .1)
    ret$pct_acc_sigv_lt01 <- mean(lapply(fit$raw_chains, `[`, 'acc_sigv') %>% unlist < .01)
    ret$pct_bad_sigv <- mean(is.na(lapply(fit$raw_chains, `[`, 'sigma_v') %>% unlist))
  }

  return(ret)
}

run_bcf.dsnum <- function(dataset, pscores, ibcf=TRUE, ubcf=FALSE, prac_est=TRUE,
                    n_chains=4, n_bcf_cores=1, nburn=1000, nsim=1000,
                    block_b0_b1=TRUE,
                    ...,
                    dir=compdir,
                    outcome = 'Y_diff',
                    splityear=FALSE,
                    byyear=FALSE,
                    additional_covdata=NULL,
                    weights=NULL) {
  gc()
  #Pretty sure we don't want weights - BCF should be deconfounding on its own
  data <- make_bcf_dataset(dataset, pscores, weights, splityear, dir=dir)
  
  #Drop dummied Xs
  #Any reason to include ypre/y_trend AND y1/y2?
  #If I had to pick one pre/trend seems like the more natural transform
  #For OLS it's all collinear, but tree stuff does care I think
  #Worst case including both sets just increases the likelihood of the tree trying to split on what we know a priori is an important var
  xmat <- data %>% select(Y_pre, Y_trend, Y_1, Y_2, size_pre, matches('^X'), matches('^V')) %>% select(-X2, -X4) %>% as.matrix
  if (splityear) {
    xmat <- cbind(xmat, data$year)
  }
  
  if (!is.null(additional_covdata)) {
    new_covars <- setdiff(colnames(additional_covdata), colnames(data))
    if ('dataset' %in% colnames(additional_covdata)) {
      additional_covdata <- additional_covdata[additional_covdata$dataset == dataset]
    } else if ('dataset.num' %in% colnames(additional_covdata)) {
      additional_covdata <- additional_covdata[additional_covdata$dataset.num == str_pad(dataset,4,'left','0')]
    }
    data <- inner_join(data, additional_covdata)
    xmat <- cbind(xmat, data %>% select(all_of(new_covars)) %>% as.matrix)
  }
  
  #For uBCF, need to feed BCF 1s instead of weights
  #but don't want to change data$w, since we still want weights used to calculate ATTs
  if (ubcf) {
    w_for_bcf <- rep(1, nrow(data))
  } else {
    w_for_bcf <- data$w
  }
  
  if (n_bcf_cores==1) {
    sink(tempfile(fileext='.log'))
  }
  fit <-      bcf(y           = data[[outcome]],
                  z           = data$Z,
                  x_control   = xmat,
                  x_moderate  = xmat,
                  pihat       = data$pscore,
                  w           = w_for_bcf,
                  verbose     = 0,
                  nburn       = nburn,
                  nsim        = nsim,
                  n_chains    = n_chains,
                  # n_cores     = n_bcf_cores,
                  n_threads   = 1,
                  # block_b0_b1 = block_b0_b1,
                  random_seed = dataset,
                  # include_random_effects = ibcf,
                  save_tree_directory = NULL,
                  log_file   = tempfile(fileext='.log'),
                  # simplified_return = TRUE,
                  ...)
  if (n_bcf_cores==1) {
    sink(NULL)
  }
  
  tau <- lapply(fit$raw_chains, function(x) x$tau) %>% abind::abind(along=3) %>% aperm(c(1,3,2))
  impacts <- bcf_impacts(data, tau, prac_est, splityear, byyear) %>%
    mutate(dataset.num = str_pad(dataset,4,'left','0')) %>%
    select(dataset.num, everything())
  
  mixing <- get_mixing(fit$raw_chains, data$w)
  
  ret <- impacts %>%
    mutate(ibcf=ibcf,
           nburn=nburn,
           nsim=nsim,
           n_chains=n_chains) %>%
    bind_cols(mixing)
  
  if (ibcf) {
    ret$avg_acc_sigv <- mean(lapply(fit$raw_chains, `[`, 'acc_sigv') %>% unlist)
    ret$pct_acc_sigv_lt10 <- mean(lapply(fit$raw_chains, `[`, 'acc_sigv') %>% unlist < .1)
    ret$pct_acc_sigv_lt01 <- mean(lapply(fit$raw_chains, `[`, 'acc_sigv') %>% unlist < .01)
    ret$pct_bad_sigv <- mean(is.na(lapply(fit$raw_chains, `[`, 'sigma_v') %>% unlist))
  }
  
  return(ret)
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
  
  if (!.version.params$sdy) data$sdY <- NULL
  
  data$Y_outcome <- data$Y_avg_diff.post_pre
  data[cols_to_remove] <- NULL
  
  weights <- data %>% 
    select(one_of(.version.params$weight_cols)) %>% 
    rowSums()
  weights <- weights / sum(weights)
  
  pscores <- data$pscore 
  data$pscore <- NULL
  
  ## Run BCF fit function
  # run_bcf <- function(dataset, pscores, ibcf=TRUE, ubcf=FALSE, prac_est=TRUE,
  #                     n_chains=4, n_bcf_cores=1, nburn=1000, nsim=1000,
  #                     block_b0_b1=TRUE,
  #                     ...,
  #                     dir=compdir,
  #                     outcome = 'Y_diff',
  #                     splityear=FALSE,
  #                     byyear=FALSE,
  #                     additional_covdata=NULL,
  #                     weights=NULL) {
    
  res <- run_bcf(
    dataset = data,
    pscores = pscores,
    ibcf = FALSE,
    outcome = 'Y_outcome'
  )
  
  return(res)

  if (!dir.exists(glue("Save/{.version.params$name}"))) dir.create(glue("Save/{.version.params$name}"))
  
  fit_save_dir <- glue("Save/{.version.params$name}/fits/")
  if (!dir.exists(fit_save_dir)) dir.create(fit_save_dir)
  saveRDS(fit, glue("{fit_save_dir}/fit_{dataset_name}_bcf.rds"))

  return(NULL)
}

# run_fit_and_pred("Data/model_ready/default/patient_wide_preproc_0001_acic_base.rds") -> res


# Loop over bcf directly --------------------------------------------------

res.l <- pmap(focus_datasets[1,], \(dataset_num, p) {
  preproc_data <- readRDS(glue("Data/model_ready/default/patient_wide_preproc_{dataset_num}_acic_base.rds"))
  
  w <- preproc_data %>% 
    select(one_of(.version.params$weight_cols)) %>% 
    rowSums()
  weights <- tibble(
    id.practice = preproc_data$id.practice,
    dataset.num = dataset_num, 
    match_wgt = 1
  )
  
  pscores <- tibble(
    pscore = preproc_data$pscore,
    id.practice = preproc_data$id.practice,
    dataset.num = dataset_num
  )
  
  res <- run_bcf.dsnum(
    dataset = dataset_num,
    pscores = pscores,
    ibcf = FALSE,
    weights = weights
  )
  return(res)
  
})
