rm(list = ls())
source("effective_support_uni_s.R")
source("asg_sampler.R")

## ------------------------------------------------------------------
## Example usage (kept minimal; reuse your kernels)
## ------------
#------------------------------------------------------
library(mcmcse)
# 1) Univariate mixture (as in your code)

k1 <- function(theta) {
  0.3 * dbeta((theta[1] + 5) / 2, 0.5, 1) / 2 +
    0.4 * dbeta((theta[1] + 1) / 2, 2, 2) / 2 +
    0.3 * dbeta((theta[1] - 5) / 2, 1, 0.5) / 2
}

L_k1 <- function(theta) {
  x <- theta[1]
  dens <- 0.3 * dbeta((x + 5) / 2, 0.5, 1) / 2 +
    0.4 * dbeta((x + 1) / 2, 2, 2) / 2 +
    0.3 * dbeta((x - 5) / 2, 1, 0.5) / 2
  if (dens <= 0 || !is.finite(dens)) {
    return(Inf)  # outside support → infinite energy
  } else {
    return(-log(dens))
  }
}


# --- Run examples
n_samples <- 1000
burn_in <- 250
thin <- 1


set.seed(123)
# Univariate
res_uni <- dim1_gibbs_sample_ASG(ker = k1, n_samples = n_samples, burn_in = burn_in, thin = thin,
  theta_init = c(0), tol = 0.01, scale0 = 1
)
res_uni$time_taken
library(coda)
mcmcse::multiESS(res_uni$samples)
mean(res_uni$samples)



# Banana

library(gridExtra)
library(RColorBrewer)

# Set seed for reproducibility
set.seed(24)

# Define the banana-shaped kernel (unnormalized density)
banana_kernel <- function(theta) {
  x1 <- theta[1]
  x2 <- theta[2]
  log_pi <- -x1^2 / 10 - x2^2 / 10 - 2 * (x2 - x1^2)^2
  exp(log_pi)
}

L_banana <- function(theta) {
  x1 <- theta[1]
  x2 <- theta[2]
  # negative log-kernel (Banana / Rosenbrock type)
  (x1^2) / 10 + (x2^2) / 10 + 2 * (x2 - x1^2)^2
}

result_banana = multivariate_gibbs_sample_ASG(banana_kernel, n_samples = n_samples,
                                          burn_in = burn_in, 
                                          thin = thin, theta_init = c(0,0))
#result_banana$samples
multiESS(result_banana$chain)
result_banana$time_taken


# Example: evaluate L and the corresponding kernel at a point
theta <- c(0, 0)
L_value <- L_banana(theta)
kernel_value <- exp(-L_value)   # K(theta) = exp(-L(theta))
cat("L(theta) at", paste(theta, collapse = ", "), "=", L_value, "\n")
cat("Kernel value exp(-L) at", paste(theta, collapse = ", "), "=", kernel_value, "\n")

# Create a grid for x1 and x2
x1_grid <- seq(-4, 4, length.out = 100)
x2_grid <- seq(-2, 6, length.out = 100)

# Compute L and exp(-L) over the grid
zL <- matrix(NA_real_, nrow = length(x1_grid), ncol = length(x2_grid))  # energy
zK <- matrix(NA_real_, nrow = length(x1_grid), ncol = length(x2_grid))  # density

for (i in seq_along(x1_grid)) {
  for (j in seq_along(x2_grid)) {
    Lij <- L_banana(c(x1_grid[i], x2_grid[j]))
    zL[i, j] <- Lij
    zK[i, j] <- exp(-Lij)
  }
}

# Prepare data frames for plotting
grid_data <- expand.grid(x1 = x1_grid, x2 = x2_grid)
grid_data$L  <- as.vector(zL)
grid_data$K  <- as.vector(zK)


# 2) Contour plot of DENSITY (kernel = exp(-L))
p_density <- ggplot(grid_data, aes(x = x1, y = x2, z = K)) +
  geom_contour_filled(bins = 20, aes(fill = after_stat(level))) +
  scale_fill_brewer(palette = "RdBu", direction = -1) +
  geom_point(aes(x = 0, y = 0), color = "green", size = 3) +
  annotate("text", x = 0, y = 0.5, label = "Mode (0, 0)", color = "green", size = 5) +
  labs(title = "Banana Target: Density Contours (K = exp(-L))", x = "x1", y = "x2") +
  theme_minimal() +
  theme(text = element_text(size = 20))

# Print plots
#p_energy
p_density

## Assuming you already have:
## - grid_data with columns x1, x2, K
## - p_density defined as in your snippet
## - sampler output in `res$samples` (n x 2 matrix)

# 0) Put samples into a data.frame and name columns
samps <- as.data.frame(result_banana$samples)
if (ncol(samps) < 2) stop("Need at least two coordinates to overlay.")
names(samps)[1:2] <- c("theta_1", "theta_2")

# (Optional) show a short recent trajectory to visualize mixing
path_n <- min(400, nrow(samps))
path_df <- samps[(nrow(samps) - path_n + 1):nrow(samps), , drop = FALSE]

# 1) Overlay samples on the existing density plot
# Overlay only sample points (no connecting lines)
p_density_overlay <- p_density +
  geom_point(
    data = samps,
    aes(x = theta_1, y = theta_2),
    inherit.aes = FALSE,
    alpha = 1,   # transparency
    color = "purple",
    size = 0.8,     # point size
    color = "black" # you can change to "white" or another color
  ) +
  coord_cartesian(
    xlim = range(grid_data$x1, na.rm = TRUE),
    ylim = range(grid_data$x2, na.rm = TRUE)
  )

# Plot it
p_density_overlay

## ------------------------------------------------------------------
## Repeated-run wrapper: runs any sampler function R times with
## different seeds, collects Multi ESS (and optionally timing), and
## returns mean/sd -- ready to drop into any table.
## ------------------------------------------------------------------
run_replicated <- function(sampler_fn, kernel, theta_init, n_samples, burn_in,
                            tol = 0.01, scale0 = 1, R = 10, base_seed = 1000,
                            ess_fun = multiESS) {
  ess_vals  <- numeric(R)
  time_vals <- numeric(R)

  for (r in seq_len(R)) {
    set.seed(base_seed + r)
    t0 <- Sys.time()
    out <- sampler_fn(kernel, theta_init = theta_init,
                       n_samples = n_samples, burn_in = burn_in,
                       tol = tol, scale0 = scale0, make_plots = FALSE)
    time_vals[r] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    ess_vals[r]  <- ess_fun(out$samples)
  }

  list(
    ess_mean  = mean(ess_vals),  ess_sd  = sd(ess_vals),
    time_mean = mean(time_vals), time_sd = sd(time_vals),
    ess_raw   = ess_vals,        time_raw = time_vals
  )
}

## Example: Rosenbrock, R = 10 replications
res_rosenbrock_rep <- run_replicated(
  sampler_fn = multivariate_gibbs_sample_ASG,
  kernel     = banana_kernel,
  theta_init = c(0, 0),
  n_samples  = 1000,
  burn_in    = 250,
  R          = 10
)

sprintf("Multi ESS: %.1f (+/- %.1f)", res_rosenbrock_rep$ess_mean, res_rosenbrock_rep$ess_sd)



# Ackley_function

ackley_f <- function(theta, a = 20, b = 0.2, c = 2*pi) {
  if (is.null(dim(theta))) {
    # vector case -> make 1-row matrix
    theta <- matrix(theta, nrow = 1)
  }
  d <- ncol(theta)
  # sum of squares per row
  ss <- rowSums(theta^2)
  # mean of cos(c * x_i) per row
  cos_mean <- rowMeans(cos(c * theta))
  
  # standard d-dim Ackley:
  # f(x) = -a * exp(-b * sqrt( (1/d) * sum x_i^2 )) - exp( mean cos(c*x_i) ) + a + e
  f <- -a * exp(-b * sqrt(ss / d)) - exp(cos_mean) + a + exp(1)
  as.numeric(f)
}

# Negative log-kernel for the Ackley target
L_ackley <- function(theta, a = 20, b = 0.2, c = 2*pi, temp = 1) {
  f <- ackley_f(theta, a = a, b = b, c = c)  # your function above
  f / temp
}


L_ackley <- function(theta, a = 20, b = 0.2, c = 2*pi, temp = 1) {
  if (is.null(dim(theta))) theta <- matrix(theta, nrow = 1)
  d  <- ncol(theta)
  ss <- rowSums(theta^2)
  cos_mean <- rowMeans(cos(c * theta))
  f <- -a * exp(-b * sqrt(ss / d)) - exp(cos_mean) + a + exp(1)  # Ackley value ≥ 0
  as.numeric(f / temp)
}


# Log-kernel: log π(θ) = - f(θ) / temp, with optional domain box constraint
ackley_log_kernel <- function(theta, temp = 1, box = c(-32.768, 32.768)) {
  if (is.null(dim(theta))) {
    theta <- matrix(theta, nrow = 1)
  }
  # domain check (hypercube)
  in_box <- apply(theta, 1, function(r) all(r >= box[1] & r <= box[2]))
  logk <- -ackley_f(theta) / temp
  logk[!in_box] <- -Inf
  as.numeric(logk)
}

# Kernel (unnormalized density): π(θ) ∝ exp( - f(θ) / temp )
# For vector theta returns scalar; for matrix returns vector (rowwise).
# Prefer using ackley_log_kernel in samplers to avoid overflow/underflow.
ackley_kernel <- function(theta, temp = 1, box = c(-32.768, 32.768)) {
  lk <- ackley_log_kernel(theta, temp = temp, box = box)
  exp(lk)
}

ackley_kernel <- compiler::cmpfun(ackley_kernel)
ackley_log_kernel <- compiler::cmpfun(ackley_log_kernel)
ackley_f <- compiler::cmpfun(ackley_f)


# -------------------------------
# Quick examples
# -------------------------------

# 1) Single point in d = 2 (global minimum at 0)
ackley_f(c(0, 0))                # should be 0
ackley_kernel(c(0, 0), temp = 1) # equals 1

# 2) Batch of points (matrix: rows are points)
set.seed(1)
X <- matrix(runif(10 * 2, -5, 5), ncol = 2)  # 10 points in R^2
ackley_f(X)
ackley_log_kernel(X, temp = 1)

# 3) Higher dimension (d = 10)
theta10 <- rnorm(10)
ackley_f(theta10)
ackley_log_kernel(theta10, temp = 2)

result_ackley_2d = multivariate_gibbs_sample_ASG(ackley_kernel ,theta_init = rep(0,2), n_samples = 1000, 
                                             burn_in = 250, thin = 1, tol = 0.1,make_plots = T)

result_ackley_2d$time_taken
multiESS(result_ackley_2d$samples)
multiESS(result_banana$chain)

result_ackley_2d$time_taken
library(ggplot2)

# --- Background: heatmap of -f (your code) -----------------------------------
xs <- seq(-5, 5, length.out = 400)
ys <- seq(-5, 5, length.out = 400)
grid <- expand.grid(x1 = xs, x2 = ys)
fvals <- ackley_f(as.matrix(grid))
grid$neg_f <- -fvals  # proxy for log-kernel

p_bg <- ggplot(grid, aes(x1, x2)) +
  geom_raster(aes(fill = neg_f), interpolate = TRUE) +
  geom_contour(aes(z = neg_f), color = "white", linewidth = 0.25, alpha = 0.6) +
  scale_fill_viridis_c(name = "-f(x)") +
  coord_fixed() +
  labs(title = "Ackley (d=2): -f(x) background with samples",
       x = "x1", y = "x2") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank())

# --- Overlay: your 1000 samples ---------------------------------------------
# If your object is named result_ackley_2d, just change the name below.
samps <- as.data.frame(result_ackley_2d$samples)
names(samps) <- c("x1", "x2")  # ensure column names

p <- p_bg +
  geom_point(data = samps, aes(x = x1, y = x2),
             inherit.aes = FALSE,
             color = "red", alpha = 0.55, size = 1)
p
## mode
d11 <- density(result_ackley_2d$samples[,1])
d22 <- density(result_ackley_2d$samples[,2])
c(d11$x[which.max(d11$y)], d22$x[which.max(d22$y)])

## =====================================================================
## Replicated results for ALL kernels (append after everything already
## defined above: k1, banana_kernel, ackley_kernel, lasso_kernel,
## bridge_kernel, euclid_norm_kernel, dim1_gibbs_sample_ASG,
## multivariate_gibbs_sample_ASG, and run_replicated() itself).
##
## R is scaled by cost: R = 10 for the cheap, low-dimensional kernels
## (univariate, Rosenbrock, Ackley, Euclidean-norm -- all <=1000
## samples, m <= 10); R = 3 for LASSO/Bridge (m = 21, 10,000 samples,
## 2,500 burn-in -- roughly 10x the cost per run of everything else).
## Raise LASSO/Bridge's R if your compute budget allows; state
## whatever R you actually used explicitly in the paper.
## =====================================================================
## ------------------------------------------------------------------
## Repeated-run wrapper: runs any sampler function R times with
## different seeds, collects Multi ESS (and optionally timing), and
## returns mean/sd -- ready to drop into any table.
## ------------------------------------------------------------------
run_replicated <- function(sampler_fn, kernel, theta_init, n_samples, burn_in,
                           tol = 0.01, scale0 = 1, R = 10, base_seed = 1000,
                           ess_fun = multiESS) {
  ess_vals  <- numeric(R)
  time_vals <- numeric(R)
 
  for (r in seq_len(R)) {
    set.seed(base_seed + r)
    t0 <- Sys.time()
    out <- sampler_fn(kernel, theta_init = theta_init,
                      n_samples = n_samples, burn_in = burn_in,
                      tol = tol, scale0 = scale0, make_plots = FALSE)
    time_vals[r] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    ess_vals[r]  <- ess_fun(out$samples)
  }
 
  list(
    ess_mean  = mean(ess_vals),  ess_sd  = sd(ess_vals),
    time_mean = mean(time_vals), time_sd = sd(time_vals),
    ess_raw   = ess_vals,        time_raw = time_vals
  )
}

# 1) Univariate mixture (as in your code)

k1 <- function(theta) {
  0.3 * dbeta((theta[1] + 5) / 2, 0.5, 1) / 2 +
    0.4 * dbeta((theta[1] + 1) / 2, 2, 2) / 2 +
    0.3 * dbeta((theta[1] - 5) / 2, 1, 0.5) / 2
}

L_k1 <- function(theta) {
  x <- theta[1]
  dens <- 0.3 * dbeta((x + 5) / 2, 0.5, 1) / 2 +
    0.4 * dbeta((x + 1) / 2, 2, 2) / 2 +
    0.3 * dbeta((x - 5) / 2, 1, 0.5) / 2
  if (dens <= 0 || !is.finite(dens)) {
    return(Inf)  # outside support → infinite energy
  } else {
    return(-log(dens))
  }
}
## Univariate ESS needs effectiveSize, not multiESS (p = 1)
ess_univariate <- function(samples) as.numeric(effectiveSize(samples[, 1]))

## ---------------------------------------------------------------
## 1) Univariate (Beta mixture, m = 1)
## ---------------------------------------------------------------
res_uni_rep <- run_replicated(
  sampler_fn = dim1_gibbs_sample_ASG,
  kernel     = k1,
  theta_init = c(0),
  n_samples  = 1000,
  burn_in    = 250,
  tol        = 0.01,
  scale0     = 1,
  R          = 10,
  ess_fun    = ess_univariate
)
sprintf("Univariate  -- ESS: %.1f (+/- %.1f)  |  time: %.4f (+/- %.4f)",
        res_uni_rep$ess_mean, res_uni_rep$ess_sd,
        res_uni_rep$time_mean, res_uni_rep$time_sd)


## ---------------------------------------------------------------
## 3) Ackley (m = 2)
## ---------------------------------------------------------------


# Negative log-kernel for the Ackley target
L_ackley <- function(theta, a = 20, b = 0.2, c = 2*pi, temp = 1) {
  f <- ackley_f(theta, a = a, b = b, c = c)  # your function above
  f / temp
}


L_ackley <- function(theta, a = 20, b = 0.2, c = 2*pi, temp = 1) {
  if (is.null(dim(theta))) theta <- matrix(theta, nrow = 1)
  d  <- ncol(theta)
  ss <- rowSums(theta^2)
  cos_mean <- rowMeans(cos(c * theta))
  f <- -a * exp(-b * sqrt(ss / d)) - exp(cos_mean) + a + exp(1)  # Ackley value ≥ 0
  as.numeric(f / temp)
}


# Log-kernel: log π(θ) = - f(θ) / temp, with optional domain box constraint
ackley_log_kernel <- function(theta, temp = 1, box = c(-32.768, 32.768)) {
  if (is.null(dim(theta))) {
    theta <- matrix(theta, nrow = 1)
  }
  # domain check (hypercube)
  in_box <- apply(theta, 1, function(r) all(r >= box[1] & r <= box[2]))
  logk <- -ackley_f(theta) / temp
  logk[!in_box] <- -Inf
  as.numeric(logk)
}

ackley_kernel <- function(theta, temp = 1, box = c(-32.768, 32.768)) {
  lk <- ackley_log_kernel(theta, temp = temp, box = box)
  exp(lk)
}

ackley_kernel <- compiler::cmpfun(ackley_kernel)
ackley_log_kernel <- compiler::cmpfun(ackley_log_kernel)
ackley_f <- compiler::cmpfun(ackley_f)

res_ackley_rep <- run_replicated(
  sampler_fn = multivariate_gibbs_sample_ASG,
  kernel     = ackley_kernel,
  theta_init = rep(0, 2),
  n_samples  = 1000,
  burn_in    = 250,
  tol        = 0.1,
  scale0     = 1,
  R          = 10
)
sprintf("Ackley      -- ESS: %.1f (+/- %.1f)  |  time: %.4f (+/- %.4f)",
        res_ackley_rep$ess_mean, res_ackley_rep$ess_sd,
        res_ackley_rep$time_mean, res_ackley_rep$time_sd)


## ===============================================================
## LASSO Kernel Function with Coefficient Path Plot 
## ===============================================================

# Required packages
library(glmnet)
library(reshape2)   # replaces pivot_longer

# ==========================================================
# LASSO Kernel for d = 8  (β0, β1, β2)
# ==========================================================

data("QuickStartExample")
QuickStartExample_data = as.data.frame(QuickStartExample)
make_lasso_kernel <- function(y, X, lambda = 1, sigma2 = 1.0, alpha = 0.5) {
  N <- length(y)
  function(theta) {
    beta0 <- theta[1]
    beta  <- theta[-1]
    N   <- length(y)
    mu <- beta0 + as.vector(X %*% beta)
    rss <- sum((y - mu)^2)
    
    # Negative log-likelihood (proportional)
    neg_log_lik <- rss / (2 * sigma2)
    
    # Negative log-prior: LASSO = Laplace prior on beta (not beta0)
    neg_log_prior <- N*lambda * sum(abs(beta)^alpha)
    
    # Unnormalized posterior ∝ exp(- (neg_log_lik + neg_log_prior))
    exp(-(neg_log_lik + neg_log_prior))
  }
}



dat <- QuickStartExample_data
y   <- QuickStartExample$y
X   <- QuickStartExample$x
#X = X[,-20]

lambda_lasso_kernel <- function(theta, sigma2 = 0.1, alpha = 1) {
  dat <- QuickStartExample_data
  y   <- QuickStartExample$y
  X   <- QuickStartExample$x
  N   <- length(y)
  
  lambda <- exp(theta[1])          # <-- simple: enforce λ ≥ 0
  beta0  <- theta[2]
  beta   <- theta[-c(1,2)]
  
  mu  <- beta0 + as.vector(X %*% beta)
  rss <- mean((y - mu)^2)
  
  neg_log_lik   <- rss / (2 * sigma2 )  # mild scaling
  neg_log_prior <- N*lambda * mean(abs(beta^alpha))
  
   exp(-(neg_log_lik + neg_log_prior))

}

lambda_lasso_L <- function(theta, sigma2 = 1.0) {
  dat <- QuickStartExample_data
  y   <- QuickStartExample$y
  X   <- QuickStartExample$x
  N   <- length(y)
  
  lambda <- exp(theta[1])          # enforce λ ≥ 0
  beta0  <- theta[2]
  beta   <- theta[-c(1,2)]
  
  mu  <- beta0 + as.vector(X %*% beta)
  rss <- mean((y - mu)^2)
  
  neg_log_lik   <- rss / (2 * sigma2)     # scaled negative log-likelihood
  neg_log_prior <- N * lambda * mean(abs(beta))  # L1 penalty term
  
  L_value <- neg_log_lik + neg_log_prior  # total energy (−log kernel)
  return(as.numeric(L_value))
}


lambda_lasso_L(rep(1,22))




lasso_kernel <- make_lasso_kernel(y, X, lambda = 0.001, alpha = 0.1)
lasso_kernel(rep(1,21))

# ==========================================================
# Example: test kernel call
# ==========================================================
set.seed(123)
n_samples <- 10000
burn_in <- 2500
thin <- 1

res_lasso = multivariate_gibbs_sample_ASG(lasso_kernel ,theta_init = rep(0,21), n_samples, burn_in, tol = 0.1, scale0 = 1)
dim(res_lasso$samples)
hist((exp(res_lasso$samples[,1])))


mean(effectiveSize(res_lasso$samples))
multiESS(res_lasso$samples)
## --- Extract posterior samples ---
S <- res_lasso$samples
if (!is.data.frame(S)) S <- as.data.frame(S)

# Ensure proper names
colnames(S) <- paste0("beta_", 0:(ncol(S) - 1))

## --- Posterior mode estimation using kernel density mode ---
posterior_mode <- function(x) {
  d <- density(x)     # smooth density
  d$x[which.max(d$y)]              # mode = x at maximum density
}

modes <- sapply(S, posterior_mode)
modes
ci_bridge <- sapply(S, quantile, probs = c(0.025, 0.975))
 ci_bridge




library(ggplot2)
library(tidyr)
library(dplyr)

res_lasso$time_taken
# Convert samples to data frame
S <- as.data.frame(res_lasso$samples)

colnames(S) <- paste0("beta_", 0:(ncol(S) - 1))


# Convert to long format (serial order)
longS <- S %>%
  mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "parameter", values_to = "value") %>%
  mutate(parameter = factor(parameter, levels = paste0("beta_", 0:(ncol(S) - 1))))

hist(res_lasso$samples[,3])

# Plot posterior densities: 7 columns × 3 rows = 21 panels
ggplot(longS, aes(x = value)) +
  geom_density(fill = "skyblue", alpha = 0.5, linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = 2, color = "red") +
  facet_wrap(~ parameter, scales = "free", ncol = 7, nrow = 3) +
  labs(
    title = "Posterior Densities of LASSO Coefficients (β0–β20)",
    x = "Coefficient value",
    y = "Density",
    subtitle = "Arranged serially (7×3 grid); dashed red line = 0"
  ) +
  theme_bw(base_size = 12)+
  coord_cartesian(xlim = c(-2, 2))

y_pred_asg <- X %*% posterior_df$PosteriorMode
y_pred_lasso <- X %*% posterior_df$Parameter



###### Bridge ######
## ===============================================================
## Bridge Regression Kernel (generalized LASSO, alpha < 1)
## Section 6.2: applies the generalized version of the lasso kernel
## (make_lasso_kernel with alpha != 1) to the same QuickStartExample
## dataset, with lambda = 0.001 and alpha = 0.1
## ===============================================================

# ==========================================================
# Bridge Regression Kernel for d = 21  (beta0, beta1, ..., beta20)
# ==========================================================

data("QuickStartExample")
QuickStartExample_data <- as.data.frame(QuickStartExample)

dat <- QuickStartExample_data
y   <- QuickStartExample$y
X   <- QuickStartExample$x

## Same kernel constructor as the LASSO case, alpha < 1 gives the
## Bridge Regression penalty sum |beta_j|^alpha instead of the L1 norm
bridge_kernel <- make_lasso_kernel(y, X, lambda = 0.001, alpha = 0.1)
bridge_kernel(rep(1, 21))

# ==========================================================
# Sampling: same setup as Section 6.1, alpha = 0.1
# ==========================================================
set.seed(123)
n_samples <- 10000
burn_in   <- 2500
thin      <- 1

res_bridge <- multivariate_gibbs_sample_ASG(
  bridge_kernel,
  theta_init = rep(0, 21),
  n_samples  = n_samples,
  burn_in    = burn_in,
  tol        = 0.1,
  scale0     = 1
)

dim(res_bridge$samples)
hist(res_bridge$samples[, 1])

mean(effectiveSize(res_bridge$samples))
multiESS(res_bridge$samples)

## --- Extract posterior samples ---
S_bridge <- res_bridge$samples
if (!is.data.frame(S_bridge)) S_bridge <- as.data.frame(S_bridge)

colnames(S_bridge) <- paste0("beta_", 0:(ncol(S_bridge) - 1))

## --- Posterior mode estimation using kernel density mode ---
posterior_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

modes_bridge <- sapply(S_bridge, posterior_mode)
modes_bridge

## --- 95% credible intervals, for the paper's Table 9 format ---
ci_bridge <- sapply(S_bridge, quantile, probs = c(0.025, 0.975))
ci_bridge



####### f(x) = ||x||^2 function #######

## Kernel: K(x) ∝ exp(-||x||), x ∈ R^m
## Section 5.3 setup: m = 10, n_samples = 1000, burn_in = 200

euclid_norm_kernel <- function(theta) {
  exp(-sqrt(sum(theta^2)))
}

set.seed(123)

m <- 10
theta_init <- rep(0.1, m)   # avoid exact 0 (||x|| = 0 is fine here, but keep away from any edge cases)

time_asg_euclid <- system.time(
  asg_out_euclid <- multivariate_gibbs_sample_ASG(
    ker        = euclid_norm_kernel,
    n_samples  = 1000,
    burn_in    = 200,
    thin       = 1,
    theta_init = theta_init,
    tol        = 0.01,
    scale0     = 1,
    make_plots = FALSE
  )
)

time_asg_euclid

## ---------------------------------------------------------
## Diagnostics matching the paper's Section 5.3 checks
## ---------------------------------------------------------
library(coda)

theta_samples <- asg_out_euclid$samples   # n_samples x m matrix

## Effective sample size per dimension (Table 7)
ess_by_dim <- apply(theta_samples, 2, effectiveSize)
print(ess_by_dim)
multiESS(theta_samples)
## Multivariate ESS (single joint summary)
mESS_euclid <- multiESS(theta_samples)
mESS_euclid

## Radial component r = ||x|| should follow Gamma(m, 1)
r_samples <- sqrt(rowSums(theta_samples^2))

## ACF of ||x|| (Figure 6a)
acf(r_samples, lag.max = 100, main = "Autocorrelation of ||x||")

## Empirical vs. theoretical Gamma(m, 1) density (Figure 6b)
library(ggplot2)
df_r <- data.frame(r = r_samples)
xs <- seq(0, max(r_samples) * 1.1, length.out = 500)
df_gamma <- data.frame(x = xs, dens = dgamma(xs, shape = m, rate = 1))

ggplot() +
  geom_histogram(data = df_r, aes(x = r, y = after_stat(density)),
                 bins = 40, fill = "skyblue", color = "black", alpha = 0.7) +
  geom_line(data = df_gamma, aes(x = x, y = dens), color = "red", linewidth = 1.2) +
  labs(title = "Distribution of ||x|| vs Theoretical Gamma(m, 1)",
       x = "r = ||x||", y = "Density") +
  theme_minimal(base_size = 12)
