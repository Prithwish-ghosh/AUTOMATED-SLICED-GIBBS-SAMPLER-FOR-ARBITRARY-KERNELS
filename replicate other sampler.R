## =====================================================================
## Replicated runs (R replications per sampler x kernel) for
## variability estimates. Fully self-contained: every kernel, every
## pseudo-target, and every tuning parameter is defined or defaulted
## inside this file, so nothing depends on what ran before it or in
## what order.
## =====================================================================

library(rhmc)
library(qslice)
library(mcmc)
library(spBayes)
library(coda)
library(mcmcse)

## =====================================================================
## Kernel definitions
## =====================================================================

## ----- Univariate: 3-component Beta mixture on R -----
k1 <- function(theta) {
  x <- theta[1]
  0.3 * dbeta((x + 5) / 2, 0.5, 1) / 2 +
    0.4 * dbeta((x + 1) / 2, 2, 2) / 2 +
    0.3 * dbeta((x - 5) / 2, 1, 0.5) / 2
}

f_k1 <- function(x) {
  dens <- k1(x)
  -log(dens + 1e-300)
}

lf_k1 <- function(x) {
  val <- k1(x)
  log(val + 1e-300)
}

ltd_k1 <- function(theta) {
  val <- k1(theta)
  log(val + 1e-300)
}

## ----- Ackley -----
ackley_f <- function(theta, a = 20, b = 0.2, c = 2 * pi) {
  if (is.null(dim(theta))) theta <- matrix(theta, nrow = 1)
  d <- ncol(theta)
  ss <- rowSums(theta^2)
  cos_mean <- rowMeans(cos(c * theta))
  f <- -a * exp(-b * sqrt(ss / d)) - exp(cos_mean) + a + exp(1)
  as.numeric(f)
}

L_ackley <- function(theta, a = 20, b = 0.2, c = 2 * pi, temp = 1) {
  ackley_f(theta, a = a, b = b, c = c) / temp
}

ackley_log_kernel <- function(theta, a = 20, b = 0.2, c = 2 * pi) {
  theta <- as.numeric(theta)
  -ackley_f(theta, a = a, b = b, c = c)
}

lf_ackley_vec <- function(theta, a = 20, b = 0.2, c = 2 * pi) {
  -ackley_f(theta, a = a, b = b, c = c)
}

make_cond_logtarget <- function(theta_current, j, a = 20, b = 0.2, c = 2 * pi) {
  force(theta_current); force(j)
  function(xj) {
    th <- theta_current
    th[j] <- xj
    lf_ackley_vec(th, a = a, b = b, c = c)
  }
}

## ----- Banana / Rosenbrock -----
banana_kernel <- function(theta) {
  x1 <- theta[1]; x2 <- theta[2]
  exp(-x1^2 / 10 - x2^2 / 10 - 2 * (x2 - x1^2)^2)
}

banana_log_kernel <- function(theta) {
  x1 <- theta[1]; x2 <- theta[2]
  -x1^2 / 10 - x2^2 / 10 - 2 * (x2 - x1^2)^2
}

L_banana <- function(theta) {
  x1 <- theta[1]; x2 <- theta[2]
  (x1^2) / 10 + (x2^2) / 10 + 2 * (x2 - x1^2)^2
}

ltd_banana <- function(theta) banana_log_kernel(theta)

lf_x1_given_x2 <- function(x1, x2_fixed) banana_log_kernel(c(x1, x2_fixed))
lf_x2_given_x1 <- function(x2, x1_fixed) banana_log_kernel(c(x1_fixed, x2))

## =====================================================================
## Pseudo-targets / supporting distributions
## =====================================================================

pseudo <- list(
  ld = function(x) dunif(x, min = -8, max = 8, log = TRUE),
  q  = function(u) qunif(u, min = -8, max = 8)
)

pseudo1d <- list(
  ld = function(x) dunif(x, min = -5, max = 5, log = TRUE),
  q  = function(u) qunif(u, min = -5, max = 5)
)

pseudo_x1 <- list(
  ld = function(x) dunif(x, min = -5, max = 5, log = TRUE),
  q  = function(u) qunif(u, min = -5, max = 5)
)

pseudo_x2 <- list(
  ld = function(x) dunif(x, min = -8, max = 8, log = TRUE),
  q  = function(u) qunif(u, min = -8, max = 8)
)

mu_k1  <- 0
Sig_k1 <- matrix(25, nrow = 1, ncol = 1)

mu_banana  <- c(0, 0)
Sig_banana <- diag(2)

mu_ackley  <- c(0, 0)
Sig_ackley <- diag(2)

## =====================================================================
## Generic replication driver
## run_once_fn(seed) must return list(samples = <n x p matrix>, time = <secs>)
## =====================================================================

R_REPS    <- 10
BASE_SEED <- 5000

run_replicated <- function(run_once_fn, R = R_REPS, base_seed = BASE_SEED,
                           label = "") {
  ess_vals  <- numeric(R)
  time_vals <- numeric(R)
  
  for (r in seq_len(R)) {
    out <- run_once_fn(base_seed + r)
    p <- ncol(out$samples)
    ess_vals[r] <- if (p == 1) {
      as.numeric(effectiveSize(out$samples[, 1]))
    } else {
      as.numeric(multiESS(out$samples))
    }
    time_vals[r] <- out$time
  }
  
  data.frame(
    label     = label,
    R         = R,
    ess_mean  = mean(ess_vals),
    ess_sd    = sd(ess_vals),
    time_mean = mean(time_vals),
    time_sd   = sd(time_vals)
  )
}

## =====================================================================
## HMC -- univariate, Ackley, banana
## All tuning parameters are explicit defaults, not free variables --
## these wrappers do not depend on anything defined outside this file.
## =====================================================================

run_hmc_univariate <- function(seed, numit = 1250, L = 10, eps = 0.05, mass = 1) {
  set.seed(seed)
  t0 <- Sys.time()
  hmc_out <- hmc(f = f_k1, init = 0, numit = numit, L = L, eps = eps, mass = mass)
  list(
    samples = matrix(hmc_out$chain[1, ], ncol = 1),
    time    = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
}

run_hmc_ackley <- function(seed, numit = 1250, L = 10, eps = 0.05, mass = 1) {
  set.seed(seed)
  t0 <- Sys.time()
  hmc_out <- hmc(f = L_ackley, init = c(0, 0), numit = numit, L = L, eps = eps, mass = mass)
  list(
    samples = t(as.matrix(hmc_out$chain)),
    time    = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
}

run_hmc_banana <- function(seed, numit = 1250, L = 10, eps = 0.05, mass = 1) {
  set.seed(seed)
  t0 <- Sys.time()
  hmc_out <- hmc(f = L_banana, init = c(0, 0), numit = numit, L = L, eps = eps, mass = mass)
  list(
    samples = t(as.matrix(hmc_out$chain)),
    time    = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
}

## =====================================================================
## Adaptive Gibbs (adaptMetropGibbs) -- univariate, Ackley, banana
## =====================================================================

run_adapgibbs_univariate <- function(seed, batch = 25, batch.length = 50) {
  set.seed(seed)
  t0 <- Sys.time()
  amg_out <- adaptMetropGibbs(
    ltd = ltd_k1, starting = c(0), tuning = 0.5, accept.rate = 0.5,
    batch = batch, batch.length = batch.length, report = 0, verbose = FALSE
  )
  list(
    samples = matrix(amg_out$p.theta.samples, ncol = 1),
    time    = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
}

run_adapgibbs_ackley <- function(seed, batch = 25, batch.length = 50) {
  set.seed(seed)
  t0 <- Sys.time()
  amg_out <- adaptMetropGibbs(
    ltd = ackley_log_kernel, starting = c(0, 0), tuning = 0.5, accept.rate = 0.5,
    batch = batch, batch.length = batch.length, report = 0, verbose = FALSE
  )
  list(
    samples = amg_out$p.theta.samples,
    time    = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
}

run_adapgibbs_banana <- function(seed, batch = 50, n_keep = 1000, n_burn = 200,
                                 tuning_init = c(0.8, 1.2)) {
  set.seed(seed)
  t0 <- Sys.time()
  n_total <- n_burn + n_keep
  batch.length <- n_total / batch
  if (batch.length != as.integer(batch.length)) {
    stop("batch and n_burn + n_keep must divide evenly for adaptMetropGibbs.")
  }
  amg_out <- adaptMetropGibbs(
    ltd = ltd_banana, starting = c(0, 0), tuning = tuning_init, accept.rate = 0.44,
    batch = batch, batch.length = batch.length, report = 0, verbose = FALSE
  )
  samps_all <- amg_out$p.theta.samples
  if (is.vector(samps_all)) samps_all <- matrix(samps_all, ncol = 2, byrow = TRUE)
  samps <- samps_all[(n_burn + 1):n_total, , drop = FALSE]
  list(
    samples = samps,
    time    = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
}

## =====================================================================
## Qslice (slice_quantile) -- univariate, Ackley, banana
## =====================================================================

run_qslice_univariate <- function(seed, n_draws = 1000) {
  set.seed(seed)
  t0 <- Sys.time()
  draws <- numeric(n_draws)
  draws[1] <- 0
  for (i in 2:n_draws) {
    out <- slice_quantile(x = draws[i - 1], log_target = lf_k1, pseudo = pseudo)
    draws[i] <- out$x
  }
  list(
    samples = matrix(draws, ncol = 1),
    time    = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
}

run_qslice_ackley <- function(seed, n_draws = 1000, d = 2) {
  set.seed(seed)
  t0 <- Sys.time()
  theta <- matrix(NA_real_, nrow = n_draws, ncol = d)
  theta[1, ] <- rep(0, d)
  for (i in 2:n_draws) {
    th <- theta[i - 1, ]
    for (j in 1:d) {
      logt_j <- make_cond_logtarget(th, j)
      out <- slice_quantile(x = th[j], log_target = logt_j, pseudo = pseudo1d)
      th[j] <- out$x
    }
    theta[i, ] <- th
  }
  list(
    samples = theta,
    time    = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
}

run_qslice_banana <- function(seed, n_draws = 1000) {
  set.seed(seed)
  t0 <- Sys.time()
  theta <- matrix(NA_real_, nrow = n_draws, ncol = 2)
  theta[1, ] <- c(0, 0)
  for (i in 2:n_draws) {
    x1_prev <- theta[i - 1, 1]
    x2_prev <- theta[i - 1, 2]
    out1 <- slice_quantile(x = x1_prev, log_target = function(x) lf_x1_given_x2(x, x2_prev), pseudo = pseudo_x1)
    x1_new <- out1$x
    out2 <- slice_quantile(x = x2_prev, log_target = function(x) lf_x2_given_x1(x, x1_new), pseudo = pseudo_x2)
    x2_new <- out2$x
    theta[i, ] <- c(x1_new, x2_new)
  }
  list(
    samples = theta,
    time    = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
}

## =====================================================================
## Elliptical Slice (slice_elliptical_mv) -- univariate, Ackley, banana
## =====================================================================

run_ellip_univariate <- function(seed, n_iter = 1000) {
  set.seed(seed)
  t0 <- Sys.time()
  draws <- numeric(n_iter)
  draws[1] <- 0
  for (i in 2:n_iter) {
    out <- slice_elliptical_mv(x = draws[i - 1], log_target = lf_k1, mu = mu_k1, Sig = Sig_k1, is_chol = FALSE)
    draws[i] <- out$x
  }
  list(
    samples = matrix(draws, ncol = 1),
    time    = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
}

run_ellip_ackley <- function(seed, n_iter = 1000, D = 2) {
  set.seed(seed)
  t0 <- Sys.time()
  draws <- matrix(0, nrow = n_iter, ncol = D)
  draws[1, ] <- rep(0, D)
  for (i in 2:n_iter) {
    out <- slice_elliptical_mv(
      x = draws[i - 1, ],
      log_target = function(x) ackley_log_kernel(x, a = 20, b = 0.2, c = 2 * pi),
      mu = mu_ackley, Sig = Sig_ackley, is_chol = FALSE
    )
    draws[i, ] <- out$x
  }
  list(
    samples = draws,
    time    = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
}

run_ellip_banana <- function(seed, n_iter = 1000) {
  set.seed(seed)
  t0 <- Sys.time()
  draws <- matrix(0, nrow = n_iter, ncol = 2)
  draws[1, ] <- c(0, 0)
  for (i in 2:n_iter) {
    out <- slice_elliptical_mv(
      x = draws[i - 1, ], log_target = banana_log_kernel,
      mu = mu_banana, Sig = Sig_banana, is_chol = FALSE
    )
    draws[i, ] <- out$x
  }
  list(
    samples = draws,
    time    = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
}

## =====================================================================
## Run everything and assemble a single results table
## NOTE: this is compute-heavy (4 methods x 3 kernels x R_REPS runs).
## Do a quick pass with R_REPS <- 2 first to confirm every wrapper runs
## cleanly end to end before committing to the full R_REPS <- 10.
## =====================================================================

results_replicated <- rbind(
  run_replicated(run_hmc_univariate,       label = "HMC -- univariate"),
  run_replicated(run_hmc_ackley,           label = "HMC -- Ackley"),
  run_replicated(run_hmc_banana,           label = "HMC -- Rosenbrock"),
  
  run_replicated(run_adapgibbs_univariate, label = "Adaptive Gibbs -- univariate"),
  run_replicated(run_adapgibbs_ackley,     label = "Adaptive Gibbs -- Ackley"),
  run_replicated(run_adapgibbs_banana,     label = "Adaptive Gibbs -- Rosenbrock"),
  
  run_replicated(run_qslice_univariate,    label = "Qslice -- univariate"),
  run_replicated(run_qslice_ackley,        label = "Qslice -- Ackley"),
  run_replicated(run_qslice_banana,        label = "Qslice -- Rosenbrock"),
  
  run_replicated(run_ellip_univariate,     label = "Elliptical slice -- univariate"),
  run_replicated(run_ellip_ackley,         label = "Elliptical slice -- Ackley"),
  run_replicated(run_ellip_banana,         label = "Elliptical slice -- Rosenbrock")
)

print(results_replicated)

## ---------------------------------------------------------------
## Table-ready strings for direct paste into LaTeX cells
## ---------------------------------------------------------------
results_replicated$ess_cell <- sprintf("%.1f $\\pm$ %.1f",
                                       results_replicated$ess_mean,
                                       results_replicated$ess_sd)
results_replicated$time_cell <- sprintf("%.3f $\\pm$ %.3f",
                                        results_replicated$time_mean,
                                        results_replicated$time_sd)
print(results_replicated[, c("label", "ess_cell", "time_cell")])

## Reminder: apply the same run_replicated() pattern to
## multivariate_gibbs_sample_ASG (ASG rows) and mcmc::metrop() (RW-MH
## rows) so every cell in Tables 2/3/6 -- not just the four baselines
## above -- reflects R_REPS replications rather than a single run.