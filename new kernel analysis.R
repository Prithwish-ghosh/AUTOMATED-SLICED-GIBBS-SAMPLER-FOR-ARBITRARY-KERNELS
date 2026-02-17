rm(list = ls())
library(coda)
library(ggplot2)
setwd("~/Downloads/ASG Sampler/code")
source("effective_support_uni_s.R")
source("asg_sampler.R")
## ============================================================
## Kernels (NOT log-kernels) from "C.5 Target Distributions"
## ============================================================

## ---------- 1) Funnel distribution kernel ----------
## p(x) = N(x_D | 0, sigma^2) * N(x_{1:D-1} | mu, exp(x_D) I )
## so sd for x_{1:D-1} is exp(x_D/2)
funnel_kernel <- function(x, sigma = 3, mu = 0) {
  x <- as.numeric(x)
  D <- 10
  if (D < 2) stop("Funnel kernel needs D >= 2.")
  
  xD <- x[D]
  sd_rest <- exp(xD / 2)
  
  d1 <- dnorm(xD, mean = 0, sd = sigma)
  d2 <- prod(dnorm(x[1:(D - 1)], mean = mu, sd = sd_rest))
  
  d1 * d2
}

funnel_log_kernel <- function(x, sigma = 3, mu = 0) {
  x <- as.numeric(x)
  D <- 10
  if (length(x) != D) stop("x must have length D=10 for this function.")
  if (D < 2) stop("Funnel kernel needs D >= 2.")
  
  xD <- x[D]
  sd_rest <- exp(xD / 2)
  
  log_d1 <- dnorm(xD, mean = 0, sd = sigma, log = TRUE)
  log_d2 <- sum(dnorm(x[1:(D - 1)], mean = mu, sd = sd_rest, log = TRUE))
  
  log_d1 + log_d2
}

funnel_neglog_kernel <- function(x, sigma = 3, mu = 0) {
  -funnel_log_kernel(x, sigma = sigma, mu = mu)
}


## ---------- 2) Hybrid Rosenbrock distribution kernel (2D) ----------
## In the screenshot (2D):
## p(x) = N(x1 | a, 1/2) * N(x2 | x1^2, 1/(2b))
## i.e. sd1 = sqrt(1/2), sd2 = sqrt(1/(2b))
hybrid_rosenbrock_kernel_2d <- function(x, a = 1, b = 100) {
  x <- as.numeric(x)
  if (length(x) != 2) stop("This function is for the 2D hybrid Rosenbrock only.")
  
  x1 <- x[1]
  x2 <- x[2]
  
  sd1 <- sqrt(1/2)
  sd2 <- sqrt(1/(2*b))
  
  dnorm(x1, mean = a, sd = sd1) * dnorm(x2, mean = x1^2, sd = sd2)
}

hybrid_rosenbrock_log_kernel_2d <- function(x, a = 1, b = 100) {
  x <- as.numeric(x)
  if (length(x) != 2) stop("This function is for the 2D hybrid Rosenbrock only.")
  
  x1 <- x[1]
  x2 <- x[2]
  
  sd1 <- sqrt(1/2)
  sd2 <- sqrt(1/(2*b))
  
  dnorm(x1, mean = a,      sd = sd1, log = TRUE) +
    dnorm(x2, mean = x1^2, sd = sd2, log = TRUE)
}

hybrid_rosenbrock_neglog_kernel_2d <- function(x, a = 1, b = 100) {
  -hybrid_rosenbrock_log_kernel_2d(x, a = a, b = b)
}

## ---------- 3) Squiggle distribution kernel ----------
## Transformation in screenshot:
## x1 = z1
## x_{2:D} = z_{2:D} - sin(a z1)
##
## Inverse:
## z1 = x1
## z_{2:D} = x_{2:D} + sin(a x1)
##
## Density:
## p(x) = N(z(x) | mu, Sigma) * |det(dz/dx)|
## and log det(dz/dx)=0 -> determinant = 1
##
## Parameters in screenshot: a = 1.5, mu=0,
## Sigma = diag(5, 1/2, ..., 1/2)
squiggle_kernel <- function(x, a = 1.5) {
  x <- as.numeric(x)
  D <- 3
  if (D < 2) stop("Squiggle kernel needs D >= 2.")
  
  z <- numeric(D)
  z[1] <- x[1]
  z[2:D] <- x[2:D] + sin(a * x[1])
  
  # mu = 0
  # diag variances: var1=5, var2..D = 1/2
  sds <- c(sqrt(5), rep(sqrt(1/2), D - 1))
  
  prod(dnorm(z, mean = 0, sd = sds))
}

squiggle_log_kernel <- function(x, a = 1.5) {
  x <- as.numeric(x)
  D <- 3
  if (length(x) != D) stop("x must have length D=3 for this function.")
  if (D < 2) stop("Squiggle kernel needs D >= 2.")
  
  z <- numeric(D)
  z[1] <- x[1]
  z[2:D] <- x[2:D] + sin(a * x[1])
  
  sds <- c(sqrt(5), rep(sqrt(1/2), D - 1))
  
  sum(dnorm(z, mean = 0, sd = sds, log = TRUE))
}

squiggle_neglog_kernel <- function(x, a = 1.5) {
  -squiggle_log_kernel(x, a = a)
}

## ---------- 4) Allen–Cahn field system kernel (unnormalized) ----------
## Screenshot gives:
## log p(x) = -beta * ( a/(2*Δs) * sum_{i=1}^{D+1} (x_i - x_{i-1})^2
##                    + b*Δs/4   * sum_{i=1}^{D}   (1 - x_i^2)^2 )
##
## with Δs = 1/D, x0 = x_{D+1} = 0, a = 0.1, b = 1/a, D fixed (e.g. 16)
##
## Kernel = exp(log p(x)) (unnormalized)
allen_cahn_kernel <- function(x, beta = 1, a = 0.1) {
  x <- as.numeric(x)
  D <- 10
  
  ds <- 1 / D
  b <- 1 / a
  
  # boundary: x0=0, x_{D+1}=0
  x_ext <- c(0, x, 0)
  
  term1 <- sum((x_ext[2:(D + 2)] - x_ext[1:(D + 1)])^2)      # i = 1..D+1
  term2 <- sum((1 - x^2)^2)                                  # i = 1..D
  
  logp <- -beta * ( (a / (2 * ds)) * term1 + (b * ds / 4) * term2 )
  
  exp(logp)
}

allen_cahn_log_kernel <- function(x, beta = 1, a = 0.1) {
  x <- as.numeric(x)
  D <- 10
  if (length(x) != D) stop("x must have length D=10 for this function.")
  
  ds <- 1 / D
  b <- 1 / a
  
  x_ext <- c(0, x, 0)  # x0=0, x_{D+1}=0
  
  term1 <- sum((x_ext[2:(D + 2)] - x_ext[1:(D + 1)])^2)  # i = 1..D+1
  term2 <- sum((1 - x^2)^2)                               # i = 1..D
  
  -beta * ( (a / (2 * ds)) * term1 + (b * ds / 4) * term2 )
}

allen_cahn_neglog_kernel <- function(x, beta = 1, a = 0.1) {
  -allen_cahn_log_kernel(x, beta = beta, a = a)
}



set.seed(2026)
n_samples = 1000
burn_in = 250
thin = 1

library(coda)

## FUNNEL DISTRIBUTION
res_funnel_kernel = multivariate_gibbs_sample_ASG(funnel_kernel, n_samples, burn_in, theta_init = rep(0,10))
effectiveSize(res_funnel_kernel$samples)
mean(effectiveSize(res_funnel_kernel$samples))
res_funnel_kernel$time_taken

## Hybrid Rosenbrock kernel 
res_hyb_rosen = multivariate_gibbs_sample_ASG(hybrid_rosenbrock_kernel_2d,n_samples, burn_in, theta_init = rep(0,2))
effectiveSize(res_hyb_rosen$samples)
mean(effectiveSize(res_hyb_rosen$samples))
res_hyb_rosen$time_taken

## squiggle_kernel
res_squiggle_kernel = multivariate_gibbs_sample_ASG(squiggle_kernel, n_samples, burn_in, theta_init = rep(0,3))
effectiveSize(res_squiggle_kernel$samples)
mean(effectiveSize(res_squiggle_kernel$samples))
res_squiggle_kernel$time_taken

## allen_cahn_kernel
res_allen_cahn_kernel = multivariate_gibbs_sample_ASG(allen_cahn_kernel, n_samples, burn_in, theta_init = rep(0,10))
effectiveSize(res_allen_cahn_kernel$samples)
mean(effectiveSize(res_allen_cahn_kernel$samples))
res_allen_cahn_kernel$time_taken



###### RW MH algorithm
library(ggplot2)
setwd("~/Downloads/ASG Sampler/code")
source("mh_initials.R") 
# Load required libraries
library(mcmc)
library(mvtnorm)  # For bivariate normal kernel
library(coda)     # For ESS and ACF plots

# Set parameters
set.seed(123)
n_samples <- 1000
burn_in <- 250
thin <- 1
nbatch <- n_samples * thin

# MH funnel
time_k_funn <- system.time({
  out_k_funn <- metrop(obj = funnel_log_kernel, initial = rep(0,10), nbatch = nbatch + burn_in, blen = 1, scale = 1)
})
mean(effectiveSize(out_k_funn$batch))
time_k_funn

# MH Hybrid Rosen 
time_k_hrk <- system.time({
  out_k_hrk <- metrop(obj = hybrid_rosenbrock_log_kernel_2d, initial = rep(0,2), nbatch = nbatch + burn_in, blen = 1, scale = 1)
})
mean(effectiveSize(out_k_hrk$batch))
time_k_hrk

# MH Squiggle
time_k_squig <- system.time({
  out_k_squig <- metrop(obj = squiggle_log_kernel, initial = rep(0,3), nbatch = nbatch + burn_in, blen = 1, scale = 1)
})
mean(effectiveSize(out_k_squig$batch))
time_k_squig

# MH ACK

time_k_ach <- system.time({
  out_k_ach <- metrop(obj = allen_cahn_log_kernel, initial = rep(0,10), nbatch = nbatch + burn_in, blen = 1, scale = 1)
})
mean(effectiveSize(out_k_ach$batch))
time_k_ach



#### HMC Analysis #####

library(rhmc)
library(qslice)
library(mcmc)
set.seed(123)
## HMC parameters
numit <- 1250   # total iterations = burn-in (250) + kept (1000)
L     <- 10     # leapfrog steps
eps   <- 0.05   # step size
mass  <- 1      # mass (scalar for 1D)


## Funnel HMC
time_hmc_funn <- system.time(
  hmc_out_funn <- hmc(
    f     = funnel_neglog_kernel,
    init  = rep(0,10),
    numit = numit,
    L     = L,
    eps   = eps,
    mass  = mass
    # num_grad = 3     # only if your implementation supports it
  )
)

mean(effectiveSize(hmc_out_funn$chain))
time_hmc_funn

# Hybrid Rosenbrock HMC
time_hmc_hrbk <- system.time(
  hmc_out_hrbk <- hmc(
    f     = hybrid_rosenbrock_neglog_kernel_2d,
    init  = rep(0,2),
    numit = numit,
    L     = L,
    eps   = eps,
    mass  = mass
    # num_grad = 3     # only if your implementation supports it
  )
)
mean(effectiveSize(hmc_out_hrbk$chain))
time_hmc_hrbk

# Squiggle HMC
time_hmc_squig <- system.time(
  hmc_out_squig <- hmc(
    f     = squiggle_neglog_kernel,
    init  = rep(0,3),
    numit = numit,
    L     = L,
    eps   = eps,
    mass  = mass
    # num_grad = 3     # only if your implementation supports it
  )
)
mean(effectiveSize(hmc_out_squig$chain))
time_hmc_squig


# Alan C HMC
time_hmc_achk <- system.time(
  hmc_out_achk <- hmc(
    f     = allen_cahn_neglog_kernel,
    init  = rep(0,10),
    numit = numit,
    L     = L,
    eps   = eps,
    mass  = mass
    # num_grad = 3     # only if your implementation supports it
  )
)
mean(effectiveSize(hmc_out_achk$chain))
time_hmc_achk


## Adabtive Gibbs sampling

library(spBayes)
set.seed(123)

n_total       <- 1250      # 250 burn-in + 1000 kept
batch         <- 25
batch.length  <- 50        # 25 * 50 = 1250


# Funnel AGS
time_adap_gibbs_funnel <- system.time(
  amg_out_funnel <- adaptMetropGibbs(
    ltd         = funnel_log_kernel,
    starting    = rep(0,10),      # initial theta
    tuning      = 0.5,       # initial proposal sd (will be adapted)
    accept.rate = 0.5,      # target acceptance
    batch       = batch,
    batch.length = batch.length,
    report      = 5,
    verbose     = TRUE
  ))
mean(effectiveSize(amg_out_funnel$p.theta.samples))
time_adap_gibbs_funnel

# Hybrid Rosenbrock AGS
time_adap_gibbs_hrbrk <- system.time(
  amg_out_hrbrk <- adaptMetropGibbs(
    ltd         = hybrid_rosenbrock_log_kernel_2d,
    starting    = rep(0,2),      # initial theta
    tuning      = 0.5,       # initial proposal sd (will be adapted)
    accept.rate = 0.5,      # target acceptance
    batch       = batch,
    batch.length = batch.length,
    report      = 5,
    verbose     = TRUE
  ))
mean(effectiveSize(amg_out_hrbrk$p.theta.samples))
time_adap_gibbs_hrbrk

# Squiggle AGS
time_adap_gibbs_sqgl <- system.time(
  amg_out_sqgl <- adaptMetropGibbs(
    ltd         = squiggle_log_kernel,
    starting    = rep(0,3),      # initial theta
    tuning      = 0.5,       # initial proposal sd (will be adapted)
    accept.rate = 0.5,      # target acceptance
    batch       = batch,
    batch.length = batch.length,
    report      = 5,
    verbose     = TRUE
  ))
mean(effectiveSize(amg_out_sqgl$p.theta.samples))
time_adap_gibbs_sqgl


# AC kernel AGS
time_adap_gibbs_achk <- system.time(
  amg_out_achk <- adaptMetropGibbs(
    ltd         = allen_cahn_log_kernel,
    starting    = rep(0,10),      # initial theta
    tuning      = 0.5,       # initial proposal sd (will be adapted)
    accept.rate = 0.5,      # target acceptance
    batch       = batch,
    batch.length = batch.length,
    report      = 5,
    verbose     = TRUE
  ))
mean(effectiveSize(amg_out_achk$p.theta.samples))
time_adap_gibbs_achk


# Qslice sampler

library(qslice)

## Funnel

pseudo_funnel <- list(
  ld = function(x) dunif(x, min = -8, max = 8, log = TRUE),
  q  = function(u) qunif(u, min = -8, max = 8)
)

n_draws_funnel <- 1000
draws_funnel <- matrix(0, nrow = n_draws_funnel, ncol = 10)
draws_funnel[1, ] <- rep(0, 10)
nEvaluations_funnel <- 0L

time_qslice_funnel <- system.time(
  for (i in 2:n_draws_funnel) {
    x_curr <- draws_funnel[i - 1, ]
    
    for (j in 1:10) {
      lf_funnel_j <- function(xj) {
        x_tmp <- x_curr
        x_tmp[j] <- xj
        funnel_log_kernel(x_tmp, sigma = 3, mu = 0)   # <-- log-kernel directly
      }
      
      out_qslice_funn <- slice_quantile(
        x          = x_curr[j],
        log_target = lf_funnel_j,
        pseudo     = pseudo_funnel
      )
      
      x_curr[j] <- out_qslice_funn$x
      nEvaluations_funnel <- nEvaluations_funnel + out_qslice_funn$nEvaluations
    }
    
    draws_funnel[i, ] <- x_curr
  }
)

mean(effectiveSize(draws_funnel))
time_qslice_funnel

## ---------- 2) HYBRID ROSENBROCK log-kernel (2D) ----------

pseudo_rosen <- list(
  ld = function(x) dunif(x, min = -8, max = 8, log = TRUE),
  q  = function(u) qunif(u, min = -8, max = 8)
)

n_draws_rosen <- 1000
draws_rosen <- matrix(0, nrow = n_draws_rosen, ncol = 2)
draws_rosen[1, ] <- c(0, 0)
nEvaluations_rosen <- 0L

time_qslice_rosen <- system.time(
  for (i in 2:n_draws_rosen) {
    x_curr <- draws_rosen[i - 1, ]
    
    for (j in 1:2) {
      lf_rosen_j <- function(xj) {
        x_tmp <- x_curr
        x_tmp[j] <- xj
        hybrid_rosenbrock_log_kernel_2d(x_tmp, a = 1, b = 100)  # <-- log-kernel
      }
      
      out <- slice_quantile(
        x          = x_curr[j],
        log_target = lf_rosen_j,
        pseudo     = pseudo_rosen
      )
      
      x_curr[j] <- out$x
      nEvaluations_rosen <- nEvaluations_rosen + out$nEvaluations
    }
    
    draws_rosen[i, ] <- x_curr
  }
)

mean(effectiveSize(draws_rosen))
time_qslice_funnel

## ---------- 3) SQUIGGLE log-kernel (D = 3) ----------

pseudo_squiggle <- list(
  ld = function(x) dunif(x, min = -8, max = 8, log = TRUE),
  q  = function(u) qunif(u, min = -8, max = 8)
)

n_draws_squiggle <- 1000
draws_squiggle <- matrix(0, nrow = n_draws_squiggle, ncol = 3)
draws_squiggle[1, ] <- c(0, 0, 0)
nEvaluations_squiggle <- 0L

time_qslice_squiggle <- system.time(
  for (i in 2:n_draws_squiggle) {
    x_curr <- draws_squiggle[i - 1, ]
    
    for (j in 1:3) {
      lf_squiggle_j <- function(xj) {
        x_tmp <- x_curr
        x_tmp[j] <- xj
        squiggle_log_kernel(x_tmp, a = 1.5)  # <-- log-kernel
      }
      
      out <- slice_quantile(
        x          = x_curr[j],
        log_target = lf_squiggle_j,
        pseudo     = pseudo_squiggle
      )
      
      x_curr[j] <- out$x
      nEvaluations_squiggle <- nEvaluations_squiggle + out$nEvaluations
    }
    
    draws_squiggle[i, ] <- x_curr
  }
)

mean(effectiveSize(draws_squiggle))
time_qslice_squiggle

## ---------- 4) ALLEN–CAHN log-kernel (D = 10) ----------

pseudo_ac <- list(
  ld = function(x) dunif(x, min = -8, max = 8, log = TRUE),
  q  = function(u) qunif(u, min = -8, max = 8)
)

n_draws_ac <- 1000
draws_ac <- matrix(0, nrow = n_draws_ac, ncol = 10)
draws_ac[1, ] <- rep(0, 10)
nEvaluations_ac <- 0L

time_qslice_ac <- system.time(
  for (i in 2:n_draws_ac) {
    x_curr <- draws_ac[i - 1, ]
    
    for (j in 1:10) {
      lf_ac_j <- function(xj) {
        x_tmp <- x_curr
        x_tmp[j] <- xj
        allen_cahn_log_kernel(x_tmp, beta = 1, a = 0.1)  # <-- log-kernel
      }
      
      out <- slice_quantile(
        x          = x_curr[j],
        log_target = lf_ac_j,
        pseudo     = pseudo_ac
      )
      
      x_curr[j] <- out$x
      nEvaluations_ac <- nEvaluations_ac + out$nEvaluations
    }
    
    draws_ac[i, ] <- x_curr
  }
)

mean(effectiveSize(draws_ac))
time_qslice_ac
