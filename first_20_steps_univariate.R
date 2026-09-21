rm(list = ls())
library(nortest)
library(ggplot2)
library(dplyr)
library(mvtnorm)
library(doParallel)
library(gridExtra)
library(coda)

source("effective_support_uni_s.R") 
source("asg sampler.R")

# Assumes these already exist in your session:
#   - effective.support()
#   - univariate_slice_sample()
#   - multivariate_gibbs_sample()
#   - k1()

library(ggplot2)
library(gridExtra)
library(coda)

set.seed(2025)
n_samples <- 1000
burn_in   <- 250
thin      <- 1

k1 <- function(theta) {
  0.3 * dbeta((theta[1] + 5) / 2, 0.5, 1) / 2 + 
    0.4 * dbeta((theta[1] + 1) / 2, 2, 2) / 2 + 
    0.3 * dbeta((theta[1] - 5) / 2, 1, 0.5) / 2
}

# Run YOUR sampler (univariate case: theta_init = c(0))
res <- multivariate_gibbs_sample_ASG(
  ker        = k1,
  n_samples  = n_samples,
  burn_in    = burn_in,
  thin       = thin,
  theta_init = c(0),
  tol        = 0.1,
  scale0     = 1
)

# Effective support + normalized target
es <- effective.support(function(x) k1(c(x)))
a  <- es$lower
b  <- es$upper
f  <- es$f

# How many first steps to overlay
N_show <- 10  # switch to 10 for first 10 steps

stopifnot(length(res$samples) >= N_show)
steps_df <- data.frame(
  Iter  = 1:N_show,
  Theta = res$samples[1:N_show]
)
steps_df$dens <- f(steps_df$Theta)   # normalized target at those steps

# ---------------- Panel A: Trace of first N steps ----------------
p_trace_steps <- ggplot(steps_df, aes(Iter, Theta)) +
  geom_path(linewidth = 0.9, color = "steelblue") +
  geom_point(size = 1.8, color = "steelblue") +
  labs(
    title = sprintf("First %d Steps (post burn-in): Trace", N_show),
    x = "Iteration (post burn-in)", y = expression(theta)
  ) +
  theme_minimal(base_size = 13)
p_trace_steps
# ---------------- Panel B: Histogram backdrop + first N steps path -------------
# Histogram uses density on y; we overlay (theta, f(theta)) so scales match
all_df <- data.frame(Theta = res$samples)

p_hist_path <- ggplot() +
  geom_histogram(data = all_df,
                 aes(x = theta_1, y = after_stat(density)),
                 bins = 60, fill = "grey85", color = "grey40") +
  # effective-support bounds
  geom_vline(xintercept = c(a, b), linetype = "dashed", color = "darkgreen") +
  # first N steps path over the histogram (uses normalized density for y)
  geom_path(data = steps_df, aes(x = Theta, y = dens),
            color = "black", linewidth = 0.9, alpha = 0.9) +
  geom_point(data = steps_df, aes(x = Theta, y = dens),
             color = "black", size = 1.6) +
  geom_rug(data = steps_df, aes(x = Theta, y = NULL), sides = "b", alpha = 0.35) +
  coord_cartesian(xlim = c(a, b)) +
  labs(
    title = sprintf("First %d Steps Over Histogram of Samples", N_show),
    x = expression(theta), y = "Density"
  ) +
  theme_minimal(base_size = 13)

p_hist_path

# Quick checks
cat("ESS (effectiveSize):", round(effectiveSize(res$samples), 1), "\n")
cat("Bounds [a,b]:", sprintf("[%.4f, %.4f]\n", a, b))
cat(sprintf("Total time recorded by your sampler: %.2f sec\n", res$time_taken))

