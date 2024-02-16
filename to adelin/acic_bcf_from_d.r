run_bart <- function(dataset, pscores,
                     n_chains=4, n_bart_cores=1,
                     nburn=1000, nsim=1000,
                     use_pscore=TRUE,
                     ...,
                     dir=compdir,
                     prac_est=TRUE,
                     outcome = 'Y_diff') {
  
  data <- make_bcf_dataset(dataset, pscores, weights=NULL, splityear=FALSE, dir=dir)
  xmat <- data %>% select(Y_pre, Y_trend, Y_1, Y_2, size_pre, matches('^X'), matches('^V')) %>% select(-X2, -X4) %>% as.matrix
  if (use_pscore) {
    xmat <- cbind(data$Z, data$pscore, xmat)
  } else {
    xmat <- cbind(data$Z, xmat)
  }
  xtest <- xmat
  #Could just pass treats and skip controls
  xtest[,'Z'] <- 1 - xtest[,'Z']
  #Bart is very responsive to scale of weights - norm to mean 1
  w <- data$w * nrow(data)/sum(data$w)
  
  fit <- dbarts::bart(x.train=xmat,
                      y.train=data[[outcome]],
                      x.test=xtest,
                      ndpost=nsim,
                      nskip=nburn,
                      nchain=n_chains,
                      nthread=n_bart_cores,
                      verbose=FALSE,
                      seed = dataset,
                      combinechains = FALSE)
  
  tau <- fit$yhat.train - fit$yhat.test
  #For controls need to flip so we get Z=1 - Z=0
  tau[,,data$Z==0] <- -tau[,,data$Z==0]
  #Redim to draws*chains*obs
  tau <- aperm(tau, c(2,1,3))
  
  impacts <- bcf_impacts(data, tau, prac_est=prac_est, splityear=FALSE, byyear=FALSE) %>%
    mutate(dataset.num = str_pad(dataset,4,'left','0')) %>%
    select(dataset.num, everything())
  
  tau_bar <- apply(tau[,,data$Z==1],c(1,2), weighted.mean, w[data$Z==1])
  yhat_bar <- fit$yhat.train %>% aperm(c(2,1,3)) %>% apply(c(1,2), weighted.mean, w)
  stanmix <- abind::abind(tau_bar, yhat_bar, t(fit$sigma), along=3)
  dimnames(stanmix)[[3]] <- c('tau_bar','yhat_bar','sigma')
  mixing <- get_mixing(chains=NULL, w=NULL, stanchains=stanmix)
  
  impacts %>%
    mutate(nburn=nburn,
           nsim=nsim,
           n_chains=n_chains) %>%
    bind_cols(mixing) %>%
    return
}

make_bcf_dataset <- function(dataset, weights=NULL, splityear=FALSE, dir=compdir) {
  
  # data <- glue('{dir}/track_2_dataset_{dataset}.RDS') %>% #mmf change 
  #   readRDS()
  
  dataset.num <- str_pad(dataset,4,'left','0')
  prac_data <- glue('track 2/practice/acic_practice_{dataset.num}.csv') %>% read_csv
  pracyr_data <- glue('track 2/practice_year/acic_practice_year_{dataset.num}.csv') %>% read_csv
  
  data <- 
    pracyr_data %>% 
    left_join(prac_data)
  
  if (splityear) {
    post_data <- data %>%
      filter(post==1) %>%
      #Joined & died are added on these files as time-invariant, so just grab one value
      group_by(id.practice, year) %>%
      summarize(Y_post = weighted.mean(Y, n.patients),
                size_post = sum(n.patients),
                across(matches('^V'), ~ weighted.mean(.x, n.patients))) %>%
      ungroup %>%
      rename_all(str_remove,'_avg')
  } else {
    post_data <- data %>%
      filter(post==1) %>%
      #  group_by(id.practice, joined, died) %>% #mmf change -- joined and died didn't help, right?
      group_by(id.practice) %>%
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
    #inner_join(pscores %>% filter(dataset.num==str_pad(dataset,4,'left','0')) %>% select(-dataset.num), by='id.practice') %>% #mmf change -- we're just doing this for one dataset
    #inner_join(pscores, by='id.practice') %>%
    #select(id.practice, year, Z, w, Y_diff, Y_post, Y_pre, Y_trend, Y_1, Y_2, size_post, size_pre, pscore, matches('X'), matches('V'), joined, died) #mmf change -- no joined/died
    select(id.practice, year, Z, w, Y_diff, Y_post, Y_pre, Y_trend, Y_1, Y_2, size_post, size_pre, matches('X'), matches('V'))
  
  if (!is.null(weights)) {
    data <- data %>%
      inner_join(weights %>% filter(dataset.num==str_pad(dataset,4,'left','0')) %>% select(-dataset.num), by='id.practice') %>%
      mutate(w = w*match_wgt) %>%
      select(-match_wgt)
  }
  
  return(data)
}

bcf_impacts <- function(data, tau, prac_est, splityear, byyear) {
  idx.trt <- data$Z==1
  d.trt <- data[idx.trt,]
  tau <- tau[,,idx.trt]
  w <- d.trt$w
  
  results <- create_sg_df()
  
  if (prac_est) {
    results <- bind_rows(results,
                         tibble(id.practice = sort(unique(d.trt$id.practice))))
  } else {
    results$id.practice = NA_real_
  }
  
  if (byyear) {
    results <- results %>%
      select(-year) %>%
      distinct %>%
      expand_grid(year=c(NA_real_, 3, 4))
  }
  
  draws <- lapply(1:nrow(results), function(i) {
    idx <- rep(TRUE, nrow(d.trt))
    if (!is.na(results$id.practice[i]))   idx <- idx & d.trt$id.practice == results$id.practice[i]
    #If we aren't doing things by year, then only overall will have values of year, and in that case we just duplicate the overall
    if (splityear & !is.na(results$year[i])) idx <- idx & d.trt$year == results$year[i]
    if (!is.na(results$level[i]))         idx <- idx & d.trt[[results$variable[i]]] == results$level[i]
    tau <- tau[,,idx, drop=FALSE]
    w <- w[idx]
    
    return(apply(tau,c(1,2), weighted.mean, w))
  }) %>%
    abind::abind(along=3)
  
  results$satt     <- apply(draws, 3, mean, na.rm=TRUE)
  results$lower.90 <- apply(draws, 3, quantile, .05, na.rm=TRUE)
  results$upper.90 <- apply(draws, 3, quantile, .95, na.rm=TRUE)
  #No clean way to remove NA for these, so they will just return NA
  results$Rhat     <- apply(draws, 3, rstan:::rhat_rfun)
  results$n_eff    <- apply(draws, 3, rstan:::ess_rfun)
  
  return(results %>%
           select(variable, level, year, id.practice, everything()) %>%
           arrange(variable, level, year, id.practice))
  
}

get_mixing <- function(chains, w, pars=NULL, simple=TRUE, wide=TRUE, stanchains=NULL) {
  if (is.null(stanchains)) {
    stanchains <- stanify(chains, w, pars, simple)
  }
  
  mixing <- rstan::monitor(stanchains,warmup=0,print=FALSE,probs=c(.05,.95)) %>%
    as_tibble(rownames='par') %>%
    mutate(width = `95%`-`5%`) %>%
    select(par, mean, sd, width, Rhat, n_eff)
  
  if (wide) {
    mixing <- mixing %>%
      pivot_wider(names_from=par, values_from=c(mean, sd, width, Rhat, n_eff), names_glue='{par}.{.value}')
  }
  
  return(mixing)
}

stanify <- function(chains, w, pars=NULL, simple=TRUE) {
  re   <- 'rho' %in% names(chains[[1]])
  tau0 <- 'tau0' %in% names(chains[[1]])
  
  if (is.null(pars) & re) {
    pars = c('sigma_y','sigma_u','sigma_v','rho','tau_bar','tau_het','mu_bar','mu_het','yhat_bar','tau_scale','mu_scale','b0','b1','delta_mu','real_mu_scale')
  } else if (is.null(pars) & !re) {
    pars = c('sigma','tau_bar','tau_het','mu_bar','mu_het','yhat_bar','tau_scale','mu_scale','b0','b1','delta_mu','real_mu_scale')
  }
  if (tau0) {
    pars = union(pars,'tau0')
  }
  
  if (simple) {
    pars <- setdiff(pars, c('tau_het','mu_het','yhat_bar','mu_scale','b0','b1','delta_mu','real_mu_scale'))
  }
  
  npar = length(pars)
  nchain = length(chains)
  ndraw = nrow(chains[[1]]$mu)
  
  stanobj <- array(0,dim=c(ndraw, nchain, npar))
  dimnames(stanobj)[[3]] <- pars
  for(i in 1:nchain) {
    j <- 0
    for (par in pars) {
      j <- j+1
      if (par=='real_mu_scale') {
        stanobj[,i,j] <- chains[[i]]$mu_scale / sqrt(chains[[i]]$delta_mu)
      } else if (par=='sigv_delta') {
        stanobj[,i,j] <- chains[[i]]$sigma_v^2 + 2*chains[[i]]$rho*chains[[i]]$sigma_v*chains[[i]]$sigma_u
        #NB: heterogeneity here is inclusive of u/v. Maybe want a taux_het and mux_het for ibcf?
      } else if (str_detect(par,'het')) {
        opar <- str_remove(par,'_het')
        stanobj[,i,j] <- sqrt(apply(chains[[i]][[opar]], 1, modi::weighted.var, w))
      } else if (!str_detect(par,'bar')) {
        stanobj[,i,j] <- chains[[i]][[par]]
      } else {
        opar <- str_remove(par,'_bar')
        stanobj[,i,j] <- matrixStats::rowWeightedMeans(chains[[i]][[opar]], w)
      }
      
    }
  }
  return(stanobj)
}

create_sg_df <- function() {
  expand_grid(variable = c('Overall', paste0('X',1:5), NA_character_),
              level=c(NA_character_,'0','1','A','B','C'),
              year=c(NA_real_, 3, 4)) %>%
    filter(is.na(year) | variable=='Overall') %>%
    filter((variable=='Overall' & is.na(level)) |
             (variable %in% c('X1','X3','X5') & level %in% c('0','1')) |
             (variable %in% c('X2','X4') & level %in% c('A','B','C')))
}

run_bcf <- function(dataset, ibcf=FALSE, ubcf=FALSE, prac_est=TRUE, #mmf change -- set ibcf to false
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
  data <- make_bcf_dataset(dataset, weights, splityear, dir=dir)
  
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
  
  pscore_fit <- dbarts::bart(x.train=xmat,
                             y.train=data$Z,
                             ndpost=nsim,
                             nskip=nburn,
                             nchain=n_chains,
                             nthread=n_bcf_cores,
                             verbose=FALSE,
                             seed = dataset,
                             combinechains = TRUE)
  pscore <- colMeans(pscore_fit$yhat.train)
  
  fit <- bcf::bcf(y           = data[[outcome]],
                  z           = data$Z,
                  x_control   = xmat,
                  x_moderate  = xmat,
                  pihat       = pscore,
                  w           = w_for_bcf,
                  verbose     = 0,
                  nburn       = nburn,
                  nsim        = nsim,
                  n_chains    = n_chains,
                  n_cores     = n_bcf_cores,
                  n_threads   = 1,
                  block_b0_b1 = block_b0_b1,
#                  random_seed = dataset,
                  random_seed = 050622,
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