rm(list = ls())
library(nortest)
library(ggplot2)
library(dplyr)
library(mvtnorm)
library(doParallel)
library(gridExtra)
library(coda)

source("effective_support_uni_s.R") 


# Multivariate Gibbs sampling function with density plots
multivariate_gibbs_sample <- function(ker, n_samples, burn_in, thin, tol = 0.1, scale0 = 1, theta_init ) {
  start_time <- Sys.time()
  
  m <- length(theta_init)
  
  total_iter <- burn_in + n_samples * thin
  theta_chain <- matrix(NA, nrow = total_iter, ncol = m)
  theta_new <- theta_init
  
  # Lists to store bounds and normalized densities for each dimension
  bounds_list <- vector("list", m)
  f_list <- vector("list", m)
  
  # Gibbs sampling loop
  for (t in 1:total_iter) {
    for (i in 1:m) {
      K_i <- function(theta_i) {
        theta_temp <- theta_new
        theta_temp[i] <- theta_i
        return(ker(theta_temp))
      }
      
      init_i <- theta_new[i]
      result <- univariate_slice_sample(ker = K_i, n_samples = 1, burn_in = 0, thin = 1, init = init_i, tol = tol, scale0 = scale0)
      theta_new[i] <- result$samples[1]
      
      # Store bounds and normalized density (only need to do this once per dimension, so we do it at t=1)
      if (t == 1) {
        bounds_list[[i]] <- list(lower = result$lower, upper = result$upper)
        f_list[[i]] <- result$f
      }
    }
    theta_chain[t, ] <- theta_new
  }
  
  sample_indices <- seq(burn_in + 1, total_iter, by = thin)
  theta_samples <- theta_chain[sample_indices, ]
  
  end_time <- Sys.time()
  time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Diagnostics and density plots for each dimension
  for (i in 1:m) {
    chain_i <- theta_chain[, i]
    total_iter <- length(chain_i)
    # Extract samples for this dimension (works for both m=1 and m>1)
    samples_i <- if (m == 1) theta_samples else theta_samples[, i]
    
    # 1. Density Plot with Bounds
    a <- bounds_list[[i]]$lower
    b <- bounds_list[[i]]$upper
    f <- f_list[[i]]
    
    x_vals <- seq(a, b, length.out = 10000)
    target_density <- f(x_vals)
    df_target <- data.frame(x = x_vals, density = target_density)
    
    p_density <- ggplot() +
      geom_histogram(aes(x = samples_i, y = after_stat(density)), bins = 50, fill = "skyblue", color = "black", alpha = 0.7) +
      geom_line(data = df_target, aes(x = x, y = density), color = "red", linewidth = 1) +
      geom_vline(xintercept = a, color = "green4", linetype = "dashed", linewidth = 1) +
      geom_vline(xintercept = b, color = "green4", linetype = "dashed", linewidth = 1) +
      labs(title = paste("Density Plot for Dimension", i, "with Bounds"),
           x = paste("theta_", i, sep = ""), y = "Density") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    # Add the time taken on the left side (fixed sprintf)
    p_density <- p_density + annotate("text", x = -Inf, y = Inf, label = sprintf("Time Taken: %.2f sec", time_taken),
                                      hjust = 0.005, vjust = 1.1, size = 5, color = "black")
    
    # 2. Trace Plot
    df_chain <- data.frame(Iteration = 1:total_iter, X = chain_i)
    p1 <- ggplot(df_chain, aes(x = Iteration, y = X)) +
      geom_line(color = "steelblue", alpha = 0.7) +
      geom_vline(xintercept = burn_in, color = "red", linetype = "dashed", linewidth = 0.5) +
      labs(title = paste("Trace Plot of MCMC Chain for Dimension", i),
           x = "Iteration", y = paste("theta_", i, sep = "")) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            panel.grid.minor = element_blank()) +
      annotate("text", x = burn_in, y = max(chain_i), label = "Burn-in", hjust = -0.1, color = "red")
    
    # 3. Autocorrelation Plot
    acf_values <- acf(chain_i, plot = FALSE, lag.max = 10)
    df_acf <- data.frame(Lag = acf_values$lag, ACF = acf_values$acf)
    p2 <- ggplot(df_acf, aes(x = Lag, y = ACF)) +
      geom_hline(yintercept = 0, color = "black") +
      geom_segment(aes(x = Lag, xend = Lag, y = 0, yend = ACF), color = "steelblue", linewidth = 0.5) +
      labs(title = paste("Autocorrelation Plot for Dimension", i),
           x = "Lag", y = "Autocorrelation") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            panel.grid.minor = element_blank())
    
    # 4. Running Mean Plot
    running_mean <- cumsum(chain_i) / (1:total_iter)
    df_running_mean <- data.frame(Iteration = 1:total_iter, RunningMean = running_mean)
    p3 <- ggplot(df_running_mean, aes(x = Iteration, y = RunningMean)) +
      geom_line(color = "darkgreen", linewidth = 1) +
      geom_vline(xintercept = burn_in, color = "red", linetype = "dashed", linewidth = 1) +
      labs(title = paste("Running Mean for Dimension", i),
           x = "Iteration", y = "Running Mean") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            panel.grid.minor = element_blank()) +
      annotate("text", x = burn_in, y = max(running_mean), label = "Burn-in", hjust = -0.1, color = "red")
    
    # Print the plots
    print(p_density)
    print(p1)
    print(p2)
    print(p3)
  }
  
  return(list(samples = theta_samples, chain = theta_chain, burn_in = burn_in, thin = thin, time_taken = time_taken))
}

# Example multivariate kernel (2D for demonstration)




## Univariate
k1 <- function(theta) {
  0.3 * dbeta((theta[1] + 5) / 2, 0.5, 1) / 2 + 
    0.4 * dbeta((theta[1] + 1) / 2, 2, 2) / 2 + 
    0.3 * dbeta((theta[1] - 5) / 2, 1, 0.5) / 2
}

bivariate_normal_ker =
  function(theta, mean = c(0, 0), sigma = matrix(c(1, 0.5, 0.5, 1), 2, 2)) {
    dmvnorm(theta, mean = mean, sigma = sigma)
  }

### Triivariate

# 3. Trivariate kernel (for testing)
K_trivariate <- function(theta) {
  0.5 * dnorm(theta[1], mean = 0, sd = 1) * dnorm(theta[2], mean = 0, sd = 1) * dnorm(theta[3], mean = 0, sd = 1) +
    0.5 * dnorm(theta[1], mean = 1, sd = 1) * dnorm(theta[2], mean = 1, sd = 1) * dnorm(theta[3], mean = 1, sd = 1)
}


k4 <- function(x) {
  0.5 * dunif(x, -30, -20) + 0.5 * dunif(x, 20, 30)
}


# Required package for bivariate normal kernel
library(mvtnorm)

# 1. Univariate Kernel: k1 (Mixture of Beta Distributions)
k1_mean <- function() {
  # Weights for the mixture
  w1 <- 0.3
  w2 <- 0.4
  w3 <- 0.3
  
  # First component: Beta(0.5, 1), transformed as theta = 2U - 5
  # E[U] = 0.5 / (0.5 + 1) = 1/3
  # E[theta_1] = 2 * E[U] - 5 = 2 * (1/3) - 5 = -13/3
  E_theta1 <- 2 * (0.5 / (0.5 + 1)) - 5
  
  # Second component: Beta(2, 2), transformed as theta = 2U - 1
  # E[U] = 2 / (2 + 2) = 0.5
  # E[theta_2] = 2 * E[U] - 1 = 2 * 0.5 - 1 = 0
  E_theta2 <- 2 * (2 / (2 + 2)) - 1
  
  # Third component: Beta(1, 0.5), transformed as theta = 2U + 5
  # E[U] = 1 / (1 + 0.5) = 2/3
  # E[theta_3] = 2 * E[U] + 5 = 2 * (2/3) + 5 = 13/3
  E_theta3 <- 2 * (1 / (1 + 0.5)) + 5
  
  # Mixture mean: E[theta] = w1 * E[theta_1] + w2 * E[theta_2] + w3 * E[theta_3]
  mean_k1 <- w1 * E_theta1 + w2 * E_theta2 + w3 * E_theta3
  
  return(mean_k1)
}

# 2. Bivariate Normal Kernel
bivariate_normal_mean <- function() {
  # The mean is directly the mean vector of the bivariate normal
  mean_vec <- c(0, 0)  # As specified in the kernel
  return(mean_vec)
}

# 3. Trivariate Kernel: K_trivariate (Mixture of Normals)
K_trivariate_mean <- function() {
  # Weights for the mixture
  w1 <- 0.5
  w2 <- 0.5
  
  # First component: Product of N(0,1) for each dimension
  mean1 <- c(0, 0, 0)
  
  # Second component: Product of N(1,1) for each dimension
  mean2 <- c(1, 1, 1)
  
  # Mixture mean: w1 * mean1 + w2 * mean2
  mean_trivariate <- w1 * mean1 + w2 * mean2
  
  return(mean_trivariate)
}

# 4. Univariate Kernel: k4 (Mixture of Uniforms)
k4_mean <- function() {
  # Weights for the mixture
  w1 <- 0.5
  w2 <- 0.5
  
  # First component: Uniform(-30, -20)
  # E[X] = (a + b) / 2 = (-30 + -20) / 2 = -25
  E_X1 <- (-30 + -20) / 2
  
  # Second component: Uniform(20, 30)
  # E[X] = (20 + 30) / 2 = 25
  E_X2 <- (20 + 30) / 2
  
  # Mixture mean: w1 * E[X1] + w2 * E[X2]
  mean_k4 <- w1 * E_X1 + w2 * E_X2
  
  return(mean_k4)
}

# Compute and print the means
cat("Mean of k1 (univariate beta mixture):", k1_mean(), "\n")
cat("Mean of bivariate_normal_ker:", bivariate_normal_mean(), "\n")
cat("Mean of K_trivariate (trivariate mixture):", K_trivariate_mean(), "\n")
cat("Mean of k4 (univariate uniform mixture):", k4_mean(), "\n")


# Run the multivariate Gibbs sampler
set.seed(123)
n_samples <- 1000
burn_in <- 250
thin <- 1
library(coda)

result1 <- multivariate_gibbs_sample(ker = k1, n_samples = n_samples,
                                     burn_in = burn_in, thin = thin, theta_init = c(0))
effectiveSize(result1$samples)
result1$time_taken
mean(result1$samples)
var(result1$samples)

#result1
result2 <- multivariate_gibbs_sample(ker = K_trivariate, n_samples = n_samples, 
                                     burn_in = burn_in, thin = thin, theta_init =  c(0,0,0))
effectiveSize(result2$samples)
result2$time_taken
c(mean(result2$samples[,1]), mean(result2$samples[,2]), mean(result2$samples[,3]))
c(var(result2$samples[,1]), var(result2$samples[,2]), var(result2$samples[,3]))


#result2
result3 <- multivariate_gibbs_sample(ker = bivariate_normal_ker, n_samples = n_samples, 
                                     burn_in = burn_in, thin = thin,  theta_init =  c(0,0))
effectiveSize(result3$samples)
result3$time_taken
#result3

result4 <- multivariate_gibbs_sample(ker =k4, n_samples = n_samples, 
                                     burn_in = burn_in, thin = thin, theta_init = c(0))
effectiveSize(result4$samples)
result4$time_taken


# 1. Univariate kernel k1
k_simple <-  function(x) {
  0.5 * dnorm(x, 0, 1) +
    0.5 * dnorm(x, 0, 1)
}

result_simp <- multivariate_gibbs_sample(ker =k_simple, n_samples = n_samples, 
                                     burn_in = burn_in, thin = thin)
effectiveSize(result_simp$samples)
result_simp$time_taken



# Load required libraries
library(coda)
library(ggplot2)
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

# Example: Evaluate the kernel at a point
theta <- c(0, 0)
kernel_value <- banana_kernel(theta)
cat("Kernel value at theta =", theta, ":", kernel_value, "\n")

# Create a grid for x1 and x2
x1_grid <- seq(-4, 4, length.out = 100)
x2_grid <- seq(-2, 6, length.out = 100)

# Compute the kernel values over the grid
z <- matrix(NA, nrow = length(x1_grid), ncol = length(x2_grid))
for (i in 1:length(x1_grid)) {
  for (j in 1:length(x2_grid)) {
    z[i, j] <- banana_kernel(c(x1_grid[i], x2_grid[j]))
  }
}

z

# Create data frame for plotting
grid_data <- expand.grid(x1 = x1_grid, x2 = x2_grid)
grid_data$z <- as.vector(z)

# Create contour plot
p <- ggplot(grid_data, aes(x = x1, y = x2, z = z)) +
  geom_contour_filled(bins = 20, aes(fill = after_stat(level))) +
  scale_fill_brewer(palette = "RdBu", direction = -1) +
  geom_point(aes(x = 0, y = 0), color = "green", size = 3) +
  annotate("text", x = 0, y = 0.5, label = "Mode (0, 0)", color = "green", size = 5) +
  labs(title = "Contour Plot of Banana-Shaped Kernel", x = "x1", y = "x2") +
  theme_minimal() +
  theme(text = element_text(size = 20))
p

result_banana = multivariate_gibbs_sample(banana_kernel, n_samples = n_samples,
                                          burn_in = burn_in, 
                                          thin = thin, theta_init = c(0,0))
result_banana$samples
effectiveSize(result_banana$samples)
result_banana$time_taken
acf(result_banana$samples[,2])
acf(result_banana$samples[,1])
c(mean(result_banana$samples[,1]),mean(result_banana$samples[,2]))

## mode
d1 <- density(result_banana$samples[,1])
d2 <- density(result_banana$samples[,2])
c(d1$x[which.max(d1$y)], d2$x[which.max(d2$y)])



# --- Libraries ---------------------------------------------------------------
library(ggplot2)

# -------------------------------
# d-dimensional Ackley + kernels
# -------------------------------

# Ackley objective (returns f(theta))
# theta: numeric vector of length d, or matrix with rows = points, cols = dims
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

# -------------------------------
# (Optional) 2D visualization
# -------------------------------
# Pretty heatmap of -f (so brighter == higher kernel mass)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  xs <- seq(-10, 10, length.out = 400)
  ys <- seq(-10, 10, length.out = 400)
  grid <- expand.grid(x1 = xs, x2 = ys)
  fvals <- ackley_f(as.matrix(grid))
  grid$neg_f <- -fvals
  
  ggplot(grid, aes(x1, x2)) +
    geom_raster(aes(fill = neg_f), interpolate = TRUE) +
    geom_contour(aes(z = neg_f), color = "white", linewidth = 0.25, alpha = 0.6) +
    scale_fill_viridis_c(name = "-f(x)") +
    coord_fixed() +
    labs(title = "Ackley (d=2): heatmap of -f (proxy for log-kernel)",
         x = "x1", y = "x2") +
    theme_minimal(base_size = 14) +
    theme(panel.grid = element_blank())
}

result_achley_2d = multivariate_gibbs_sample(ackley_kernel ,theta_init = c(0,0), n_samples = n_samples, 
                                          burn_in = burn_in, thin = thin)
effectiveSize(result_achley_2d$samples)

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
samps <- as.data.frame(result_achley_2d$samples)
names(samps) <- c("x1", "x2")  # ensure column names

p <- p_bg +
  geom_point(data = samps, aes(x = x1, y = x2),
             inherit.aes = FALSE,
             color = "red", alpha = 0.55, size = 1)

p

## mode
d1_ach <- density(result_achley_2d$samples[,1])
d2_ach <- density(result_achley_2d$samples[,2])
c(d1_ach$x[which.max(d1_ach$y)], d2_ach$x[which.max(d2_ach$y)])
