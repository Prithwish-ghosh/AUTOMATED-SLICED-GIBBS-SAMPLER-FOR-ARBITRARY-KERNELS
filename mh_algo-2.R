rm(list = ls())
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

# Run MH for k1 with timing
time_k_ackley <- system.time({
  out_k_ackley <- metrop(obj = ackley_log_kernel, initial = c(0,0), nbatch = nbatch + burn_in, blen = 1, scale = 1)
})

samples_k_ackley <- out_k_ackley$batch
samples_k_ackley
# Compute ESS and ACF
ess_k_ackley <- effectiveSize(samples_k_ackley)
ess_k_ackley
acf_k_ackley <- acf(samples_k_ackley, plot = F)

# Plotting
#par(mfrow = c(2, 2))

# --- Background: heatmap of -f (your code) -----------------------------------
xs <- seq(-5, 5, length.out = 400)
ys <- seq(-5, 5, length.out = 400)
grid <- expand.grid(x1 = xs, x2 = ys)
fvals <- ackley_log_kernel(as.matrix(grid))
grid$neg_f <- fvals  # proxy for log-kernel

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
samps <- as.data.frame(samples_k_ackley)
names(samps) <- c("x1", "x2")  # ensure column names

p <- p_bg +
  geom_point(data = samps, aes(x = x1, y = x2),
             inherit.aes = FALSE,
             color = "red", alpha = 0.55, size = 1)

p

plot(samples_k_ackley[,1], type = "l", main = "Trace plot for k_ackley dim 1", xlab = "Iteration", ylab = "Theta")
plot(samples_k_ackley[,2], type = "l", main = "Trace plot for k_ackley dim 2", xlab = "Iteration", ylab = "Theta")
acf(samples_k_ackley[,1], main = "ACF for k_ackley dim 1")

acf_vals <- acf(samples_k_ackley[, 1], plot = FALSE)

plot(
  acf_vals$lag,
  acf_vals$acf,
  type = "h",       # vertical lines
  lwd = 3,          # <-- increase thickness here
  main = "ACF for k_ackley",
  xlab = "Lag",
  ylab = "ACF"
)

acf(samples_k_ackley[,2], main = "ACF for k_ackley dim 2")
plot(running_mean(samples_k_ackley[,1]),  lwd = 3, type = "l", main = "Running Mean for k_ackley function", xlab = "Iteration", ylab = "Mean")
plot(running_mean(samples_k_ackley[,2]), type = "l", main = "Running Mean for k_ackley dim2", xlab = "Iteration", ylab = "Mean")
# Summary

cat("Time taken (seconds):", time_k_ackley[3], "\n")
cat("Effective Sample Size:", ess_k_ackley, "\n")


d1_ach_mh <- density(samples_k_ackley[,1])
d2_ach_mh <- density(samples_k_ackley[,2])
c(d1_ach_mh$x[which.max(d1_ach_mh$y)], d2_ach_mh$x[which.max(d2_ach_mh$y)])





##### Second Univariate
k1 <- function(x) {
  0.3 * dbeta((x + 5) / 2, 0.5, 1) / 2 +
    0.4 * dbeta((x + 1) / 2, 2, 2) / 2 +
    0.3 * dbeta((x - 5) / 2, 1, 0.5) / 2
}

log_k1 <- function(theta) {
  dens <- k1(theta[1])
  if (dens <= 0) return(-Inf)
  log(dens)
}

# Run MH for k1 with timing
time_k1 <- system.time({
  out_k1 <- metrop(obj = log_k1, initial = 0, nbatch = nbatch + burn_in, blen = 1, scale = 1)
})

effectiveSize(out_k1$batch)
out_k1$final
samples_k1 <- out_k1$batch

# Compute ESS and ACF
ess_k1 <- effectiveSize(samples_k1)
acf_k1 <- acf(samples_k1, plot = F)

# Plotting
#par(mfrow = c(2, 2))
plot(samples_k1, type = "l", main = "Trace plot for k1", xlab = "Iteration", ylab = "Theta")
acf(samples_k1, main = "ACF for k1")

acf_vals <- acf(samples_k1, plot = FALSE)

plot(
  acf_vals$lag,
  acf_vals$acf,
  type = "h",       # vertical lines
  lwd = 3,          # <-- increase thickness here
  main = "ACF for Univariate Kernel",
  xlab = "Lag",
  ylab = "ACF"
)

plot(running_mean(samples_k1),lwd = 2, type = "l", main = "Running Mean for k1", xlab = "Iteration", ylab = "Mean")
mean(samples_k1)


# Common x-range for both hist and curve
xrange <- c(-5, 7)
xs <- seq(xrange[1], xrange[2], length.out = 2000)

# Evaluate target curve safely over the full range
kvals <- k1(xs)
kvals[!is.finite(kvals)] <- NA

# Pre-compute histogram (no plotting) to get densities for ylim
h_tmp <- hist(samples_k1, breaks = 50, plot = FALSE)
ylim_max <- max(h_tmp$density, kvals, na.rm = TRUE)
if (!is.finite(ylim_max) || ylim_max <= 0) ylim_max <- max(h_tmp$density, na.rm = TRUE)
ylim_max <- ylim_max * 1.05


hist(samples_k1, breaks = 50, freq = FALSE,
     xlim = xrange, ylim = c(0, ylim_max),
     main = "Density for k1", xlab = "Theta",
     col = "lightblue", border = "white")

lines(xs, kvals, col = "red", lwd = 2)  # overlay curve across full [-7,7]


# Summary
cat("k1 Results:\n")
cat("Acceptance rate:", out_k1$accept, "\n")
cat("Time taken (seconds):", time_k1[3], "\n")
cat("Effective Sample Size:", ess_k1, "\n")
mean(exp(samples_k1))
var(samples_k1)
#assess_mixing(out_k1$accept, ess_k1, acf_k1, samples_k1, is_multimodal = TRUE)

# 2. Bivariate normal kernel
bivariate_normal_ker <- function(theta, mean = c(0, 0), sigma = matrix(c(1, 0.5, 0.5, 1), 2, 2)) {
  dmvnorm(theta, mean = mean, sigma = sigma)
}

log_bivariate <- function(theta) {
  dens <- bivariate_normal_ker(theta)
  if (dens <= 0) return(-Inf)
  log(dens)
}

# Run MH for bivariate
time_bivariate <- system.time({
  out_bivariate <- metrop(obj = log_bivariate, initial = c(0, 0), nbatch = nbatch + burn_in, blen = 1, scale = 0.5)
})


samples_bivariate <- out_bivariate$batch
samples_bivariate <- samples_bivariate

# Compute ESS and ACF
ess_bivariate <- effectiveSize(samples_bivariate)
acf_bivariate <- acf(samples_bivariate, plot = FALSE)

# Plotting
#par(mfrow = c(2, 3))
plot(samples_bivariate[, 1], samples_bivariate[, 2], pch = 20, main = "Bivariate Samples", xlab = "Theta1", ylab = "Theta2")
acf(samples_bivariate[, 1], main = "ACF for Theta1 ")
acf(samples_bivariate[, 2], main = "ACF for Theta2 ")
plot(samples_bivariate[, 1], type = "l", main = "Trace for Theta1 ", xlab = "Iteration", ylab = "Theta1")
plot(samples_bivariate[, 2], type = "l", main = "Trace for Theta2 ", xlab = "Iteration", ylab = "Theta2")
plot(running_mean(samples_bivariate[, 1]), type = "l", main = "Running Mean for Theta1 ", xlab = "Iteration", ylab = "Mean")
plot(running_mean(samples_bivariate[, 2]), type = "l", main = "Running Mean for Theta2", xlab = "Iteration", ylab = "Mean")

# Summary
cat("Bivariate Normal Results:\n")
cat("Acceptance rate:", out_bivariate$accept, "\n")
cat("Time taken (seconds):", time_bivariate[3], "\n")
cat("Effective Sample Size (Theta1, Theta2):", ess_bivariate, "\n")
assess_mixing(out_bivariate$accept, ess_bivariate, acf_bivariate, samples_bivariate)

# 3. Trivariate kernel
K_trivariate <- function(theta) {
  0.5 * dnorm(theta[1], mean = 0, sd = 1) * dnorm(theta[2], mean = 0, sd = 1) * dnorm(theta[3], mean = 0, sd = 1) +
    0.5 * dnorm(theta[1], mean = 1, sd = 1) * dnorm(theta[2], mean = 1, sd = 1) * dnorm(theta[3], mean = 1, sd = 1)
}

log_K_trivariate <- function(theta) {
  dens <- K_trivariate(theta)
  if (dens <= 0) return(-Inf)
  log(dens)
}

# Run MH for trivariate
time_trivariate <- system.time({
  out_trivariate <- metrop(obj = log_K_trivariate, initial = c(0, 0, 0), nbatch = nbatch + burn_in, blen = 1, scale = 0.5)
})


samples_trivariate <- out_trivariate$batch
samples_trivariate <- samples_trivariate

# Compute ESS and ACF
ess_trivariate <- effectiveSize(samples_trivariate)
acf_trivariate <- acf(samples_trivariate, plot = FALSE)

# Plotting
#par(mfrow = c(2, 3))
plot(samples_trivariate[, 1], samples_trivariate[, 2], pch = 20, main = "Theta1 vs Theta2 ", xlab = "Theta1", ylab = "Theta2")
acf(samples_trivariate[, 1], main = "ACF for Theta1 ")
acf(samples_trivariate[, 2], main = "ACF for Theta2 ")
acf(samples_trivariate[, 3], main = "ACF for Theta3 ")
plot(samples_trivariate[, 1], type = "l", main = "Trace for Theta1 ", xlab = "Iteration", ylab = "Theta1")
plot(samples_trivariate[, 2], type = "l", main = "Trace for Theta2 ", xlab = "Iteration", ylab = "Theta2")
plot(samples_trivariate[, 3], type = "l", main = "Trace for Theta3 ", xlab = "Iteration", ylab = "Theta3")
plot(running_mean(samples_trivariate[, 1]), type = "l", main = "Running Mean for Theta1 ", xlab = "Iteration", ylab = "Mean")
plot(running_mean(samples_trivariate[, 2]), type = "l", main = "Running Mean for Theta2 ", xlab = "Iteration", ylab = "Mean")
plot(running_mean(samples_trivariate[, 3]), type = "l", main = "Running Mean for Theta3 ", xlab = "Iteration", ylab = "Mean")

# Summary
cat("Trivariate Results:\n")
cat("Acceptance rate:", out_trivariate$accept, "\n")
cat("Time taken (seconds):", time_trivariate[3], "\n")
cat("Effective Sample Size (Theta1, Theta2, Theta3):", ess_trivariate, "\n")
c(mean(samples_trivariate[,1]), mean(samples_trivariate[,2]), mean(samples_trivariate[,3]))
assess_mixing(out_trivariate$accept, ess_trivariate, acf_trivariate, samples_trivariate, is_multimodal = TRUE)

# 4. Univariate kernel k4
k4 <- function(x) {
  0.5 * dunif(x, -30, -20) + 0.5 * dunif(x, 20, 30)
}

log_k4 <- function(x) {
  dens <- k4(x[1])
  if (dens <= 0) return(-Inf)
  log(dens)
}

# Run MH for k4
time_k4 <- system.time({
  out_k4 <- metrop(obj = log_k4, initial = 20, nbatch = nbatch + burn_in, blen = 1, scale = 5)
})

samples_k4 <- out_k4$batch
samples_k4 <- samples_k4

# Compute ESS and ACF
ess_k4 <- effectiveSize(samples_k4)
acf_k4 <- acf(samples_k4, plot = FALSE)

# Plotting
#par(mfrow = c(2, 2))
plot(samples_k4, type = "l", main = "Trace plot for k4 ", xlab = "Iteration", ylab = "Theta")
hist(samples_k4, breaks = 50, main = "Density for k4 ", xlab = "Theta", freq = FALSE)
curve(k4, from = -30, to = 30, add = TRUE, col = "red", lwd = 2)
acf(samples_k4, main = "ACF for k4 ")
plot(running_mean(samples_k4), type = "l", main = "Running Mean for k4 ", xlab = "Iteration", ylab = "Mean")

# Summary
cat("k4 Results:\n")
cat("Acceptance rate:", out_k4$accept, "\n")
cat("Time taken (seconds):", time_k4[3], "\n")
cat("Effective Sample Size:", ess_k4, "\n")
assess_mixing(out_k4$accept, ess_k4, acf_k4, samples_k4, is_multimodal = TRUE)



### Banana Shaped Kernel

# Load required libraries
library(mcmc)
library(coda)
library(ggplot2)
library(RColorBrewer)

# Set parameters
set.seed(123)
n_samples <- 1000
burn_in <- 250
thin <- 1
nbatch <- n_samples * thin

# Function to compute running mean
running_mean <- function(x) {
  cumsum(x) / seq_along(x)
}

# Banana-shaped kernel (unnormalized density)
banana_kernel <- function(theta) {
  x1 <- theta[1]
  x2 <- theta[2]
  log_pi <- -x1^2 / 10 - x2^2 / 10 - 2 * (x2 - x1^2)^2
  exp(log_pi)
}

# Log-density for MH sampling
log_banana <- function(theta) {
  dens <- banana_kernel(theta)
  if (dens <= 0) return(-Inf)
  log(dens)
}

# Run MH for banana kernel with timing
time_banana <- system.time({
  out_banana <- metrop(obj = log_banana, initial = c(0, 0), nbatch = nbatch + burn_in, blen = 1, scale = 0.5)
})
samples_banana <- out_banana$batch

# Compute ESS and ACF
ess_banana <- effectiveSize(samples_banana)
acf_banana <- acf(samples_banana, plot = FALSE)

# Create grid for contour plot
x1_grid <- seq(-4, 4, length.out = 100)
x2_grid <- seq(-2, 6, length.out = 100)
z <- matrix(NA, nrow = length(x1_grid), ncol = length(x2_grid))
for (i in 1:length(x1_grid)) {
  for (j in 1:length(x2_grid)) {
    z[i, j] <- banana_kernel(c(x1_grid[i], x2_grid[j]))
  }
}
grid_data <- expand.grid(x1 = x1_grid, x2 = x2_grid)
grid_data$z <- as.vector(z)

# Plotting
par(mfrow = c(1, 1))

# 2D scatter plot of samples
plot(samples_banana[, 1], samples_banana[, 2], pch = 20, main = "Banana Samples", xlab = "Theta1", ylab = "Theta2")

# Contour plot with samples
plot(NULL, xlim = c(-4, 4), ylim = c(-2, 6), main = "Banana Kernel with Samples", xlab = "Theta1", ylab = "Theta2")
contour(x1_grid, x2_grid, z, levels = exp(seq(-10, 0, by = 0.5)), col = "blue", add = TRUE)
curve(x^2, from = -4, to = 4, col = "red", lty = 2, lwd = 2, add = TRUE)
points(0, 0, pch = 19, col = "green", cex = 1.5)
text(0, 0.5, "Mode (0, 0)", col = "green", pos = 4)
points(samples_banana[, 1], samples_banana[, 2], pch = 20, col = "black", cex = 0.8)

# ACF plots
acf(samples_banana[, 1], main = "ACF for Theta1")
acf(samples_banana[, 2], main = "ACF for Theta2")

acf_vals <- acf(samples_banana[, 1], plot = FALSE)

plot(
  acf_vals$lag,
  acf_vals$acf,
  type = "h",       # vertical lines
  lwd = 3,          # <-- increase thickness here
  main = "ACF for Rosenbrock Kernel",
  xlab = "Lag",
  ylab = "ACF"
)

# Trace plots
plot(samples_banana[, 1], type = "l", main = "Trace for Theta1", xlab = "Iteration", ylab = "Theta1")
plot(samples_banana[, 2], type = "l", main = "Trace for Theta2", xlab = "Iteration", ylab = "Theta2")

# Running mean plots
plot(running_mean(samples_banana[, 1]), type = "l",lwd = 2, main = "Running Mean for Rosenbrock Kernel", xlab = "Iteration", ylab = "Mean")
plot(running_mean(samples_banana[, 2]), type = "l", main = "Running Mean for Theta2", xlab = "Iteration", ylab = "Mean")

# Reset par
par(mfrow = c(1, 1))

# ggplot2 contour plot with samples (for consistency with previous style)
p_contour <- ggplot(grid_data, aes(x = x1, y = x2)) +
  geom_contour_filled(aes(z = z, fill = after_stat(level)), bins = 20) +
  scale_fill_brewer(palette = "RdBu", direction = -1) +
  geom_point(aes(x = 0, y = 0), color = "green", size = 3) +
  geom_point(data = data.frame(x1 = samples_banana[, 1], x2 = samples_banana[, 2]), 
             aes(x = x1, y = x2), color = "purple", size = 1, alpha = 0.8) +
  labs(title = "Contour Plot of Banana-Shaped Kernel with MH Samples", x = "Theta1", y = "Theta2") +
  theme_minimal() +
  theme(text = element_text(size = 20))
p_contour


set.seed(1144)

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

grid.arrange(p_contour, p, ncol = 2)
c(mean(samples_banana[,1]),mean(samples_banana[,2]))

## Mode
d1 <- density(samples_banana[,1])
d2 <- density(samples_banana[,2])
c(d1$x[which.max(d1$y)], d2$x[which.max(d2$y)])


# Summary
cat("Banana Kernel Results:\n")
cat("Acceptance rate:", out_banana$accept, "\n")
cat("Time taken (seconds):", time_banana[3], "\n")
cat("Effective Sample Size (Theta1, Theta2):", ess_banana, "\n")



#### Lasso #####

log_lasso_kernel <- (function() {
  dat <- mtcars
  y   <- dat$mpg
  X   <- cbind(hp = dat$hp, wt = dat$wt)
  X   <- scale(X, center = TRUE, scale = TRUE)
  N   <- length(y)
  lambda <- 0.25
  
  function(theta) {
    beta0 <- theta[1]
    beta  <- theta[-1]                 # (β1, β2)
    rss   <- sum((y - beta0 - as.vector(X %*% beta))^2)
    penalty <- lambda * sum(abs(beta))
    # log of unnormalized target:
    -(rss / (2 * N) + penalty)
  }
})()


set.seed(123)
n_samples <- 10000
burn_in <- 2500
thin <- 1
nbatch <- n_samples * thin

time_k_lasso <- system.time({
  out_k_lasso <- metrop(obj = log_lasso_kernel, initial = c(0,0,0), nbatch = nbatch + burn_in, blen = 1, scale = 1)
})

#out_k_lasso$batch

estimated_samples_mh = c(mean(out_k_lasso$batch[,1]), mean(out_k_lasso$batch[,2]), mean(out_k_lasso$batch[,3]))
estimated_samples_mh
effectiveSize(out_k_lasso$batch)

acf_lasso_dim1 <- acf(out_k_lasso$batch[,1], lag.max = 50)
acf_lasso_dim2 <- acf(out_k_lasso$batch[,2], lag.max = 50)
acf_lasso_dim3 <- acf(out_k_lasso$batch[,3], lag.max = 50)
# Plotting
#par(mfrow = c(2, 2))
plot(out_k_lasso$batch[,1], type = "l", main = "Trace plot for lasso beta_0", xlab = "Iteration", ylab = "Theta")
plot(out_k_lasso$batch[,2], type = "l", main = "Trace plot for lasso beta_1", xlab = "Iteration", ylab = "Theta")
plot(out_k_lasso$batch[,3], type = "l", main = "Trace plot for lasso beta_2", xlab = "Iteration", ylab = "Theta")

plot(running_mean(out_k_lasso$batch[,2]), type = "l", main = "Running Mean for k1", xlab = "Iteration", ylab = "Mean")
mean(samples_k1)
