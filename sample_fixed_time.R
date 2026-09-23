library(coda)
library(mcmc)
library(mcmcse)
library(ggplot2)

## =====================================================================
## 1) SAMPLE-FIXED ASG sampler
##    Runs exactly burn_in + n_samples * thin sweeps, then drops burn-in
##    and thins, so `samples` has exactly n_samples rows.
## =====================================================================
multivariate_gibbs_sample_ASG_n <- function(ker, n_samples,
                                            theta_init,
                                            burn_in = 200, thin = 1,
                                            tol = 0.01, scale0 = 1) {
  t0 <- Sys.time()
  
  theta_new <- as.numeric(theta_init)
  m <- length(theta_new)
  if (m < 1) stop("theta_init must be a non-empty numeric vector.")
  if (any(!is.finite(theta_new))) stop("theta_init contains NA/Inf.")
  
  ## Precompute bounds at start (same as time-based version)
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
  
  n_iter <- burn_in + n_samples * thin          # total sweeps known up front
  chain_mat <- matrix(NA_real_, nrow = n_iter, ncol = m)
  
  for (t in seq_len(n_iter)) {
    
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
    
    chain_mat[t, ] <- theta_new
  }
  colnames(chain_mat) <- paste0("theta_", seq_len(m))
  
  idx <- seq(from = burn_in + 1L, to = n_iter, by = thin)
  samples <- chain_mat[idx, , drop = FALSE]
  
  time_taken <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  list(samples = samples, chain = chain_mat, burn_in = burn_in, thin = thin,
       n_samples = n_samples, time_taken = time_taken)
}

## =====================================================================
## 2) SAMPLE-FIXED MH (mcmc::metrop), same burn-in
## =====================================================================
mh_n <- function(initial, log_target, n_samples, burn_in = 200, scale = 1) {
  start_time <- Sys.time()
  
  out <- metrop(obj = log_target, initial = initial,
                nbatch = burn_in + n_samples, blen = 1, scale = scale)
  
  samples <- as.matrix(out$batch)[-(1:burn_in), , drop = FALSE]
  
  list(
    samples    = samples,                       # n_samples x d
    accept     = out$accept,
    elapsed    = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  )
}

## =====================================================================
## 3) Helper: run both samplers for one kernel over the N grid
## =====================================================================
ess_of <- function(X) {
  X <- as.matrix(X)
  if (ncol(X) == 1) as.numeric(effectiveSize(X)) else as.numeric(multiESS(X))
}

run_fixed_n <- function(ker, log_ker, theta_init, N_grid, kernel_name,
                        burn_in = 200, seed = 2026) {
  set.seed(seed)
  rows <- list()
  
  for (N in N_grid) {
    ## ASG
    res_p <- multivariate_gibbs_sample_ASG_n(ker = ker, n_samples = N,
                                             theta_init = theta_init,
                                             burn_in = burn_in)
    ess_p <- ess_of(res_p$samples)
    rows[[length(rows) + 1]] <- data.frame(
      Method = "ASG", kernel = kernel_name, N = nrow(res_p$samples),
      ESS = ess_p, ESS_over_N = ess_p / nrow(res_p$samples),
      time = res_p$time_taken, ESS_per_sec = ess_p / res_p$time_taken
    )
    
    ## MH
    res_mh <- mh_n(initial = theta_init, log_target = log_ker,
                   n_samples = N, burn_in = burn_in)
    ess_mh <- ess_of(res_mh$samples)
    rows[[length(rows) + 1]] <- data.frame(
      Method = "MH", kernel = kernel_name, N = nrow(res_mh$samples),
      ESS = ess_mh, ESS_over_N = ess_mh / nrow(res_mh$samples),
      time = res_mh$elapsed, ESS_per_sec = ess_mh / res_mh$elapsed
    )
    
    cat(kernel_name, "| N =", N, "| ASG ESS/N =", round(ess_p / N, 4),
        "| MH ESS/N =", round(ess_mh / N, 4), "\n")
  }
  
  out <- do.call(rbind, rows)
  out$log10N <- log10(out$N)
  out
}

N_grid <- c(1000, 3000, 5000)

## =====================================================================
## 4) Kernels (same as time-based script)
## =====================================================================

## ---- Univariate mixture ----
k1 <- function(theta) {
  0.3 * dbeta((theta[1] + 5) / 2, 0.5, 1) / 2 +
    0.4 * dbeta((theta[1] + 1) / 2, 2, 2) / 2 +
    0.3 * dbeta((theta[1] - 5) / 2, 1, 0.5) / 2
}
log_k1 <- function(theta) {
  dens <- k1(theta[1])
  if (dens <= 0) return(-Inf)
  log(dens)
}

## ---- Banana / Rosenbrock ----
banana_kernel <- function(theta) {
  x1 <- theta[1]; x2 <- theta[2]
  exp(-x1^2 / 10 - x2^2 / 10 - 2 * (x2 - x1^2)^2)
}
log_banana <- function(theta) {
  dens <- banana_kernel(theta)
  if (dens <= 0) return(-Inf)
  log(dens)
}

## ---- Ackley ----
ackley_f <- function(theta, a = 20, b = 0.2, c = 2*pi) {
  if (is.null(dim(theta))) theta <- matrix(theta, nrow = 1)
  d <- ncol(theta)
  ss <- rowSums(theta^2)
  cos_mean <- rowMeans(cos(c * theta))
  as.numeric(-a * exp(-b * sqrt(ss / d)) - exp(cos_mean) + a + exp(1))
}
ackley_log_kernel <- function(theta, temp = 1, box = c(-32.768, 32.768)) {
  if (is.null(dim(theta))) theta <- matrix(theta, nrow = 1)
  in_box <- apply(theta, 1, function(r) all(r >= box[1] & r <= box[2]))
  logk <- -ackley_f(theta) / temp
  logk[!in_box] <- -Inf
  as.numeric(logk)
}
ackley_kernel <- function(theta, temp = 1, box = c(-32.768, 32.768)) {
  exp(ackley_log_kernel(theta, temp = temp, box = box))
}

## =====================================================================
## 5) Run
## =====================================================================
res_uni    <- run_fixed_n(k1, log_k1, theta_init = 0, N_grid = N_grid,
                          kernel_name = "univariate", seed = 2025)
res_banana <- run_fixed_n(banana_kernel, log_banana, theta_init = c(0, 0),
                          N_grid = N_grid, kernel_name = "Rosenbrock Kernel", seed = 2026)
res_ackley <- run_fixed_n(ackley_kernel, ackley_log_kernel, theta_init = c(0, 0),
                          N_grid = N_grid, kernel_name = "Ackley Kernel", seed = 2026)

all_dataset_final_n <- rbind(res_uni, res_banana, res_ackley)
all_dataset_final_n

## =====================================================================
## 6) Plots
## =====================================================================

## ESS/N vs log10(N)  (analogue of the time-based plot)
ggplot(all_dataset_final_n,
       aes(x = log10N, y = ESS_over_N,
           color = Method, shape = factor(N), linetype = kernel)) +
  geom_point(size = 5) +
  geom_line(aes(group = interaction(Method, kernel)), linewidth = 1) +
  labs(x = "Log10(N)", y = "Multi ESS / N",
       shape = "Samples (N)", color = "Method", linetype = "Kernel") +
  theme_minimal(base_size = 14)

## ESS per second (accounts for per-iteration cost at fixed N)
ggplot(all_dataset_final_n,
       aes(x = factor(N), y = ESS_per_sec, fill = Method)) +
  geom_col(position = "dodge") +
  facet_wrap(~ kernel, scales = "free_y") +
  labs(x = "Samples (N, after burn-in of 200)", y = "ESS per second") +
  theme_minimal(base_size = 14)