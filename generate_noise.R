# Basic functions to add flexible random noise to a dataset
# Adelin Kassler, 11/16/23

# Basic noise generation: creates a vector
generate_noise <- function(n, distr = "random", ..., verbose = FALSE) {
  
  if(distr == "random") {
    distr <- sample(c("rnorm", "rpois", "rexp", "rgamma", "rbinom"), 1)
  }
  
  args <- list(...)
  
  if(distr == "rnorm") {
    if(is.null(args$mean)) {
      mu <- rnorm(1, mean = 0, sd = 5)
    } else {
      mu <- args$mean
    }
    if(is.null(args$sd)) {
      sigma <- runif(1, min = 0.5, max = 5)
    } else {
      sigma <- args$sd
    }
    noise <- do.call(rnorm, c(list(n = n, mean = mu, sd = sigma)))
    
    if(verbose) {
      message(paste0("Generated noise from normal distribution with mean ", mu, " and sd ", sigma))
    }
    
  } else if(distr == "rpois") {
    if(is.null(args$lambda)) {
      lambda <- sample(1:10, 1) 
    } else {
      lambda <- args$lambda
    }
    noise <- do.call(rpois, c(list(n = n, lambda = lambda)))
    
    if(verbose) {
      message(paste0("Generated noise from Poisson distribution with lambda ", lambda)) 
    }
    
  } else if(distr == "rexp") {
    if(is.null(args$rate)) {
      rate <- runif(1, min = 0.1, max = 1)
    } else {
      rate <- args$rate
    }
    noise <- do.call(rexp, c(list(n = n, rate = rate)))
    
    if(verbose) {
      message(paste0("Generated noise from exponential distribution with rate ", rate))
    }
    
  } else if(distr == "rgamma") {
    if(is.null(args$shape)) {
      shape <- sample(1:5, 1)
    } else {
      shape <- args$shape
    }
    if(is.null(args$scale)) {
      scale <- runif(1, min = 0.5, max = 5)
    } else {
      scale <- args$scale
    }
    noise <- do.call(rgamma, c(list(n = n, shape = shape, scale = scale)))
    
    if(verbose) {
      message(paste0("Generated noise from gamma distribution with shape ", shape, " and scale ", scale))
    }
    
  } else if(distr == "rbinom") {
    if(is.null(args$size)) {
      trials <- sample(5:20, 1)
    } else {
      trials <- args$size
    }
    if(is.null(args$prob)) {
      prob <- runif(1, min = 0.1, max = 0.9)
    } else {
      prob <- args$prob
    }
    noise <- do.call(rbinom, c(list(n = n, size = trials, prob = prob)))
    
    if(verbose) {
      message(paste0("Generated noise from binomial distribution with trials ", trials, " and probability ", prob))
    }
    
  } else {
    stop(paste0("Don't know how to generate from function '", distr, "'. Please set distr to one of: random, rnorm, rpois, rexp, rgamma, rbinom"))
    
  }
  
  return(noise)
  
}

# Generate a bunch of noise variables in one go, and add them on to the dataset. Compatible with dplyr pipelines.
add_noise_cols <- function(.data, num_vars, distr = "random", 
                           params = NULL, var_prefix = "noise_",
                           verbose = FALSE) {
  
  if(!is.data.frame(.data))
    stop(".data must be a data frame")
  
  if(length(num_vars) != 1 || !is.numeric(num_vars) || num_vars < 1)
    stop("num_vars must be a positive integer")
  
  if(!is.character(distr))
    stop("distr must be a character string specifying distribution type")
  
  if(!is.null(params)) {
    if(!is.list(params))
      stop("params must be a list")
    if(length(params) != num_vars)
      stop("Length of params must equal num_vars")
  }
  
  distr_params <- map(1:num_vars, ~list(distr = distr))
  
  if(!is.null(params)) {
    distr_params <- map2(distr_params, params, c)
  }
  
  noise_vars <- map(1:num_vars, function(i) {
    
    noise <- generate_noise(nrow(.data), distr_params[[i]], verbose)
    
    return(noise)
    
  })
  
  # Bind together results
  var_names <- paste0(var_prefix, 1:num_vars)
  names(noise_vars) <- var_names
  noise_df <- do.call(cbind, noise_vars)
  .data <- cbind(.data, noise_df)
  
  # Ensure correct typing to pass to BART function
  class(.data) <- c("data.frame", setdiff(class(.data), "data.frame"))
  
  return(.data)
  
}