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


# Alternative BCF ---------------------------------------------------------

# bcf <- function (y, z, x_control, x_moderate = x_control, pihat, w = NULL, 
#           random_seed = sample.int(.Machine$integer.max, 1), n_chains = 4, 
#           n_cores = n_chains, n_threads = max((RcppParallel::defaultNumThreads() - 
#                                                  2)/n_cores, 1), nburn, nsim, nthin = 1, update_interval = 100, 
#           ntree_control = 200, sd_control = NULL, base_control = 0.95, 
#           power_control = 2, ntree_moderate = 50, sd_moderate = NULL, 
#           base_moderate = 0.25, power_moderate = 3, save_tree_directory = ".", 
#           log_file = file.path(".", sprintf("bcf_log_%s.txt", format(Sys.time(), 
#                                                                      "%Y%m%d_%H%M%S"))), nu = 3, lambda = NULL, sigq = 0.9, 
#           sighat = NULL, include_pi = "control", use_muscale = TRUE, 
#           use_tauscale = TRUE, include_random_effects = FALSE, batch_size = 100, 
#           block_v_rho = FALSE, block_batch_size = 100, block_b0_b1 = FALSE, 
#           sigu_hyperprior = NULL, ate_prior_sd = NULL, hardcode_sigma_u = FALSE, 
#           hardcode_sigma_v = FALSE, hardcode_rho = FALSE, hardcode_sigma_u_val = 0, 
#           hardcode_sigma_v_val = 0, hardcode_rho_val = 0, rho_beta_prior = TRUE, 
#           rho_beta_a = 2, rho_beta_b = 2, simplified_return = FALSE, 
#           verbose = 1) 
# {
#   if (is.null(w)) {
#     w <- matrix(1, ncol = 1, nrow = length(y))
#   }
#   pihat = as.matrix(pihat)
#   if (!.ident(length(y), length(z), length(w), nrow(x_control), 
#               nrow(x_moderate), nrow(pihat))) {
#     stop("Data size mismatch. The following should all be equal:\n         length(y): ", 
#          length(y), "\n", "length(z): ", length(z), "\n", 
#          "length(w): ", length(w), "\n", "nrow(x_control): ", 
#          nrow(x_control), "\n", "nrow(x_moderate): ", nrow(x_moderate), 
#          "\n", "nrow(pihat): ", nrow(pihat), "\n")
#   }
#   if (any(is.na(y))) 
#     stop("Missing values in y")
#   if (any(is.na(z))) 
#     stop("Missing values in z")
#   if (any(is.na(w))) 
#     stop("Missing values in w")
#   if (any(is.na(x_control))) 
#     stop("Missing values in x_control")
#   if (any(is.na(x_moderate))) 
#     stop("Missing values in x_moderate")
#   if (any(is.na(pihat))) 
#     stop("Missing values in pihat")
#   if (any(!is.finite(y))) 
#     stop("Non-numeric values in y")
#   if (any(!is.finite(z))) 
#     stop("Non-numeric values in z")
#   if (any(!is.finite(w))) 
#     stop("Non-numeric values in w")
#   if (any(!is.finite(x_control))) 
#     stop("Non-numeric values in x_control")
#   if (any(!is.finite(x_moderate))) 
#     stop("Non-numeric values in x_moderate")
#   if (any(!is.finite(pihat))) 
#     stop("Non-numeric values in pihat")
#   if (!all(sort(unique(z)) == c(0, 1))) 
#     stop("z must be a vector of 0's and 1's, with at least one of each")
#   if (!(include_random_effects %in% c(TRUE, FALSE))) 
#     stop("include_random_effects must be TRUE or FALSE")
#   if (round(batch_size) != batch_size | batch_size < 1) 
#     stop("batch_size must be an integer larger than 0")
#   if (!(verbose %in% 0:4)) 
#     stop("verbose must be an integer from 0 to 4")
#   if (length(unique(y)) < 5) 
#     warning("y appears to be discrete")
#   if (nburn < 0) 
#     stop("nburn must be positive")
#   if (nsim < 0) 
#     stop("nsim must be positive")
#   if (nthin < 0) 
#     stop("nthin must be positive")
#   if (nthin > nsim + 1) 
#     stop("nthin must be < nsim")
#   if (nburn < 1000) 
#     warning("A low (<1000) value for nburn was supplied")
#   if (nsim * nburn < 1000) 
#     warning("A low (<1000) value for total iterations after burn-in was supplied")
#   if ((hardcode_rho != hardcode_sigma_v) & block_v_rho) 
#     stop("One of rho and sigma_v is hardcoded, but ablock update was specified")
#   x_c = matrix(x_control, ncol = ncol(x_control))
#   x_m = matrix(x_moderate, ncol = ncol(x_moderate))
#   if (include_pi == "both" | include_pi == "control") {
#     x_c = cbind(x_control, pihat)
#   }
#   if (include_pi == "both" | include_pi == "moderate") {
#     x_m = cbind(x_moderate, pihat)
#   }
#   cutpoint_list_c = lapply(1:ncol(x_c), function(i) .cp_quantile(x_c[, 
#                                                                      i]))
#   cutpoint_list_m = lapply(1:ncol(x_m), function(i) .cp_quantile(x_m[, 
#                                                                      i]))
#   sdy = sqrt(Hmisc::wtd.var(y, w))
#   muy = stats::weighted.mean(y, w)
#   yscale = (y - muy)/sdy
#   if (is.null(lambda)) {
#     if (is.null(sighat)) {
#       lmf = lm(yscale ~ z + as.matrix(x_c), weights = w)
#       sighat = summary(lmf)$sigma
#     }
#     qchi = qchisq(1 - sigq, nu)
#     lambda = (sighat * sighat * qchi)/nu
#   }
#   if (is.null(sd_control)) {
#     con_sd <- 2
#   }
#   else {
#     con_sd = sd_control/sdy
#   }
#   if (is.null(sd_moderate)) {
#     mod_sd <- 1/ifelse(use_tauscale, 0.674, 1)
#   }
#   else {
#     mod_sd = sd_moderate/sdy/ifelse(use_tauscale, 0.674, 
#                                     1)
#   }
#   hardcode_sigma_u_val <- hardcode_sigma_u_val/sdy
#   hardcode_sigma_v_val <- hardcode_sigma_v_val/sdy
#   if (is.null(sigu_hyperprior)) {
#     sigu_hyperprior <- con_sd/3
#   }
#   else {
#     sigu_hyperprior <- sigu_hyperprior/sdy/0.674
#   }
#   if (include_random_effects && is.null(ate_prior_sd)) {
#     stop("a prior SD for the ATE (ate_prior_sd) is required for iBCF")
#   }
#   else if (include_random_effects) {
#     ate_prior_sd <- ate_prior_sd/sdy
#   }
#   else {
#     ate_prior_sd <- 1
#   }
#   dir = tempdir()
#   perm = order(z, decreasing = TRUE)
#   RcppParallel::setThreadOptions(numThreads = n_threads)
#   do_type_config <- .get_do_type(n_cores, log_file)
#   `%doType%` <- do_type_config$doType
#   chain_out <- foreach::foreach(iChain = 1:n_chains) %doType% 
#     {
#       this_seed = random_seed + iChain - 1
#       cat("Calling bcfoverparRcppClean From R\n")
#       set.seed(this_seed)
#       tree_files = .get_chain_tree_files(save_tree_directory, 
#                                          iChain)
#       fitbcf = bcfoverparRcppClean(y_ = yscale[perm], 
#                                    z_ = z[perm], w_ = w[perm], x_con_ = t(x_c[perm, 
#                                                                               , drop = FALSE]), x_mod_ = t(x_m[perm, , drop = FALSE]), 
#                                    x_con_info_list = cutpoint_list_c, x_mod_info_list = cutpoint_list_m, 
#                                    burn = nburn, nd = nsim, thin = nthin, ntree_mod = ntree_moderate, 
#                                    ntree_con = ntree_control, lambda = lambda, 
#                                    nu = nu, con_sd = con_sd, mod_sd = mod_sd, mod_alpha = base_moderate, 
#                                    mod_beta = power_moderate, con_alpha = base_control, 
#                                    con_beta = power_control, treef_con_name_ = tree_files$con_trees, 
#                                    treef_mod_name_ = tree_files$mod_trees, status_interval = update_interval, 
#                                    use_mscale = use_muscale, use_bscale = use_tauscale, 
#                                    b_half_normal = TRUE, randeff = include_random_effects, 
#                                    batch_size = batch_size, acceptance_target = 0.44, 
#                                    verbose = verbose, block_v_rho = block_v_rho, 
#                                    block_batch_size = block_batch_size, block_b0_b1 = block_b0_b1, 
#                                    sigu_hyperprior = sigu_hyperprior, ate_prior_sd = ate_prior_sd, 
#                                    hardcode_sigma_u = hardcode_sigma_u, hardcode_sigma_v = hardcode_sigma_v, 
#                                    hardcode_rho = hardcode_rho, hardcode_sigma_u_val = hardcode_sigma_u_val, 
#                                    hardcode_sigma_v_val = hardcode_sigma_v_val, 
#                                    hardcode_rho_val = hardcode_rho_val, rho_beta_prior = rho_beta_prior, 
#                                    rho_beta_a = rho_beta_a, rho_beta_b = rho_beta_b)
#       cat("bcfoverparRcppClean returned to R\n")
#       ac = fitbcf$m_post[, order(perm)]
#       Tm = fitbcf$b_post[, order(perm)] * (1/(fitbcf$b1 - 
#                                                 fitbcf$b0))
#       Tc = ac * (1/fitbcf$msd)
#       tau_post = sdy * fitbcf$b_post[, order(perm)]
#       mu_post = muy + sdy * (Tc * fitbcf$msd + Tm * fitbcf$b0)
#       yhat_post = muy + sdy * fitbcf$yhat_post[, order(perm)]
#       u_post = sdy * fitbcf$u[, order(perm)]
#       v_post = sdy * fitbcf$v[, order(perm)]
#       if (include_random_effects) {
#         tau_post <- tau_post + v_post
#         mu_post <- mu_post + u_post
#         yhat_post <- yhat_post + u_post + t(t(v_post) * 
#                                               z)
#       }
#       sigma_i = sdy * fitbcf$sigma_i[, order(perm)]
#       names(fitbcf$acceptance) = c("sigma_y", "sigma_u", 
#                                    "sigma_v", "rho")
#       list(sigma_y = sdy * fitbcf$sigma_y, sigma_u = sdy * 
#              fitbcf$sigma_u, sigma_v = sdy * fitbcf$sigma_v, 
#            rho = fitbcf$rho, sigma_i = sigma_i, yhat = yhat_post, 
#            sdy = sdy, con_sd = con_sd, mod_sd = mod_sd, 
#            muy = muy, mu = mu_post, tau = tau_post, u = u_post, 
#            v = v_post, mu_scale = fitbcf$msd, tau_scale = fitbcf$bsd, 
#            b0 = fitbcf$b0, b1 = fitbcf$b1, delta_mu = fitbcf$delta_con, 
#            acceptance = fitbcf$acceptance, perm = perm, 
#            include_pi = include_pi, include_random_effects = include_random_effects, 
#            random_seed = this_seed)
#     }
#   if (!include_random_effects) {
#     chain_out <- lapply(chain_out, function(x) {
#       x$sigma <- x$sigma_y
#       x$sigma_y <- x$sigma_a <- x$sigma_b <- x$sigma_v <- x$rho <- x$sigma_i <- x$u <- x$v <- x$acceptance <- x$acc_sigv <- x$mux <- x$taux <- NULL
#       return(x)
#     })
#   }
#   fitObj <- list(raw_chains = chain_out)
#   if (!simplified_return) {
#     fitObj <- c(fitObj, list(coda_chains = .extract_coda_chains(chain_out)), 
#                 .get_components_from_chains(chain_out))
#   }
#   attr(fitObj, "class") <- "bcf"
#   .cleanup_after_par(do_type_config)
#   return(fitObj)
# }

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
