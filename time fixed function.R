library(coda)
library(ggplot2)
setwd("~/Downloads/ASG Sampler/code")
source("effective_support_uni_s.R")
source("asg_sampler.R")

## ----------------------------
## 1) TIME-BASED ASG sampler
## ----------------------------
multivariate_gibbs_sample_ASG_time <- function(ker, duration_secs,
                                               theta_init,
                                               burn_in = 0, thin = 1,
                                               tol = 0.01, scale0 = 1,
                                               max_store = Inf) {
  t0 <- Sys.time()
  
  theta_new <- as.numeric(theta_init)
  m <- length(theta_new)
  if (m < 1) stop("theta_init must be a non-empty numeric vector.")
  if (any(!is.finite(theta_new))) stop("theta_init contains NA/Inf.")
  
  ## Precompute bounds at start (like your ASG)
  bounds_list <- vector("list", m)
  for (i in seq_len(m)) {
    g_i <- (function(i) {
      function(xi) {
        tmp <- theta_new
        tmp[i] <- xi
        ker(tmp)
      }
    })(i)
    bnds <- effective.support(g_i, tol = tol, scale0 = scale0)
    bounds_list[[i]] <- c(a = bnds$lower, b = bnds$upper)
  }
  
  ## Store chain rows in a list (fast-ish append)
  chain_list <- vector("list", 10000)
  t <- 0L
  
  ## Run until time budget
  while (as.numeric(difftime(Sys.time(), t0, units = "secs")) < duration_secs) {
    
    ## ONE u for the full sweep (ASG style)
    kcur <- ker(theta_new)
    if (!is.finite(kcur) || kcur <= 0) {
      stop("ker(theta) must be positive and finite at current state; got ", kcur,
           ". Check kernel definition / initial state.")
    }
    u <- runif(1, 0, kcur)
    
    for (i in seq_len(m)) {
      g_i <- (function(i) {
        function(xi) {
          tmp <- theta_new
          tmp[i] <- xi
          ker(tmp)
        }
      })(i)
      
      a_i <- bounds_list[[i]]["a"]; b_i <- bounds_list[[i]]["b"]
      theta_new[i] <- slice1d_fixed_u(g_i, u, a = a_i, b = b_i)
    }
    
    t <- t + 1L
    if (t > length(chain_list)) {
      chain_list <- c(chain_list, vector("list", length(chain_list))) # double capacity
    }
    chain_list[[t]] <- theta_new
    
    if (t >= max_store) break
  }
  
  if (t < 2) stop("Too few iterations collected within duration_secs; increase duration_secs.")
  
  chain_mat <- do.call(rbind, chain_list[1:t])
  colnames(chain_mat) <- paste0("theta_", seq_len(m))
  
  ## burn-in + thinning (on time-based iterations)
  idx <- seq(from = burn_in + 1L, to = t, by = thin)
  idx <- idx[idx >= 1 & idx <= t]
  if (length(idx) < 2) stop("After burn_in/thin, too few samples. Reduce burn_in or thin.")
  
  samples <- chain_mat[idx, , drop = FALSE]
  
  time_taken <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  list(samples = samples, chain = chain_mat, burn_in = burn_in, thin = thin, time_taken = time_taken)
}

## ----------------------------
## 2) ESS/N for prefixes
##    - multivariate: min ESS across dims (conservative)
## ----------------------------
ess_ratio_prefixes <- function(X, N_grid) {
  X <- as.matrix(X)
  nmax <- nrow(X)
  
  N_grid <- sort(unique(pmin(N_grid, nmax)))
  N_grid <- N_grid[N_grid >= 50]
  
  out <- data.frame(N = integer(0), log10N = numeric(0),
                    ess = numeric(0), ess_over_N = numeric(0))
  
  for (N in N_grid) {
    mcm <- coda::mcmc(X[1:N, , drop = FALSE])
    ess_vec <- coda::effectiveSize(mcm)
    ess_val <- if (length(ess_vec) == 1) as.numeric(ess_vec) else min(as.numeric(ess_vec), na.rm = TRUE)
    
    out <- rbind(out, data.frame(
      N = N,
      log10N = log10(N),
      ess = ess_val,
      ess_over_N = ess_val / N
    ))
  }
  out
}

################### MH algorithm ######################

library(mcmc)
library(coda)


mh_time <- function(initial, log_target, scale = 1, run_time ) {
  start_time <- Sys.time()
  batch_size = 1000
  current <- initial
  samples <- list()
  i <- 1
  
  while (as.numeric(difftime(Sys.time(), start_time, units = "secs")) < run_time) {
    # Do a batch at once
    out <- metrop(obj = log_target, initial = current,
                  nbatch = batch_size, blen = 1, scale = scale)
    
    samples[[i]] <- out$batch
    current <- out$final
    i <- i + 1
  }
  
  end_time <- Sys.time()
  total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  list(
    samples = do.call(rbind, samples),   # all samples
    elapsed = total_time                 # actual runtime in seconds
  )
}