rm(list = ls())
source("effective_support_uni_s.R")
source("asg_sampler.R")

## ------------------------------------------------------------------
## Example usage (kept minimal; reuse your kernels)
## ------------------------------------------------------------------

# 1) Univariate mixture (as in your code)
k1 <- function(theta) {
  0.3 * dbeta((theta[1] + 5) / 2, 0.5, 1) / 2 +
    0.4 * dbeta((theta[1] + 1) / 2, 2, 2) / 2 +
    0.3 * dbeta((theta[1] - 5) / 2, 1, 0.5) / 2
}

# 2) Bivariate normal kernel (example)
library(mvtnorm)
bivariate_normal_ker <- function(theta, mean = c(0, 0),
                                 sigma = matrix(c(1, 0.6, 0.6, 1), 2, 2)) {
  dmvnorm(theta, mean = mean, sigma = sigma)
}

# --- Run examples

n_samples <- 1000
burn_in <- 250
thin <- 1


set.seed(123)
# Univariate
res_uni <- dim1_gibbs_sample_ASG(
  ker = k1, n_samples = n_samples, burn_in = burn_in, thin = thin,
  theta_init = c(0), tol = 0.01, scale0 = 1
)
res_uni$time_taken
library(coda)
effectiveSize(res_uni$samples)
mean(res_uni$samples)
# Bivariate normal
res_bvn <- multivariate_gibbs_sample_ASG(
  ker = bivariate_normal_ker, n_samples = n_samples, burn_in = burn_in, thin = thin,
  theta_init = c(0, 0), tol = 0.01, scale0 = 1
)
effectiveSize(res_bvn$samples)


## ---- Parameters ------------------------------------------------------------
mu    <- c(0, 0)
Sigma <- matrix(c(1, 0.6, 0.6, 1), 2, 2)

## ---- Grid for contours -----------------------------------------------------
# Choose plotting window using marginal sds (wide enough to cover most mass)
sd_x <- sqrt(Sigma[1, 1])
sd_y <- sqrt(Sigma[2, 2])

x_seq <- seq(mu[1] - 4 * sd_x, mu[1] + 4 * sd_x, length.out = 200)
y_seq <- seq(mu[2] - 4 * sd_y, mu[2] + 4 * sd_y, length.out = 200)

grid <- expand.grid(x = x_seq, y = y_seq)
grid$density <- bivariate_normal_ker(cbind(grid$x, grid$y), mean = mu, sigma = Sigma)



p <- ggplot() +
  # Contour lines (True Kernel)
  geom_contour(data = grid, aes(x = x, y = y, z = density, colour = "True Kernel"), bins = 12) +
  # Samples (points)
  geom_point(data = res_bvn$samples, aes(x = theta_1, y = theta_2, colour = "Samples"), size = 0.8, alpha = 0.5) +
  coord_equal() +
  labs(
    title = "Bivariate Normal Kernel: Contours with Sample Overlay",
    subtitle = sprintf("mean = (%.1f, %.1f),  Sigma = [[%.1f, %.1f],[%.1f, %.1f]]",
                       mu[1], mu[2], Sigma[1,1], Sigma[1,2], Sigma[2,1], Sigma[2,2]),
    x = expression(theta[1]),
    y = expression(theta[2]),
    colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = c(0.85, 0.85),       # top-right corner
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  ) +
  scale_colour_manual(
    values = c("True Kernel" = "darkblue", "Samples" = "red")
  )

print(p)

# Trivariate
# 3. Trivariate kernel (for testing)
K_trivariate <- function(theta) {
  0.5 * dnorm(theta[1], mean = 0, sd = 1) * dnorm(theta[2], mean = 0, sd = 1) * dnorm(theta[3], mean = 0, sd = 1) +
    0.5 * dnorm(theta[1], mean = 1, sd = 1) * dnorm(theta[2], mean = 1, sd = 1) * dnorm(theta[3], mean = 1, sd = 1)
}

res_tri = multivariate_gibbs_sample_ASG(ker = K_trivariate, n_samples = n_samples, burn_in = burn_in, 
                                        theta_init = c(0,0,0), tol = 0.01, scale0 = 1)

res_tri$time_taken
effectiveSize(res_tri$samples)
mean(res_tri$samples[,1])
mean(res_tri$samples[,2])
mean(res_tri$samples[,3])

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

#z

# Create data frame for plotting
grid_data <- expand.grid(x1 = x1_grid, x2 = x2_grid)
grid_data$z <- as.vector(z)

# Create contour plot
p_banana <- ggplot(grid_data, aes(x = x1, y = x2, z = z)) +
  geom_contour_filled(bins = 20, aes(fill = after_stat(level))) +
  scale_fill_brewer(palette = "RdBu", direction = -1) +
  geom_point(aes(x = 0, y = 0), color = "green", size = 3) +
  annotate("text", x = 0, y = 0.5, label = "Mode (0, 0)", color = "green", size = 5) +
  labs(title = "Contour Plot of Banana-Shaped Kernel", x = "x1", y = "x2") +
  theme_minimal() +
  theme(text = element_text(size = 20))
p_banana

result_banana = multivariate_gibbs_sample_ASG(banana_kernel, n_samples = n_samples,
                                          burn_in = burn_in, 
                                          thin = thin, theta_init = c(0,0))
#result_banana$samples
effectiveSize(result_banana$samples)
result_banana$time_taken

## mode
d1 <- density(result_banana$samples[,1])
d2 <- density(result_banana$samples[,2])
c(d1$x[which.max(d1$y)], d2$x[which.max(d2$y)])

samps_banana <- as.data.frame(result_banana$samples)
names(samps_banana) <- c("x1", "x2")  # ensure column names

p_banana_final <- p_banana +
  geom_point(data = samps_banana, aes(x = x1, y = x2),
             inherit.aes = FALSE,
             color = "purple", alpha = 0.55, size = 0.8)

p_banana_final






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

result_ackley_2d = multivariate_gibbs_sample_ASG(ackley_kernel ,theta_init = c(0,0), n_samples = n_samples, 
                                             burn_in = burn_in, thin = thin)
effectiveSize(result_ackley_2d$samples)
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




## ===============================================================
## LASSO Kernel Function with Coefficient Path Plot (Fixed Version)
## ===============================================================

# Required packages
library(glmnet)
library(reshape2)   # replaces pivot_longer

# ==========================================================
# LASSO Kernel for d = 2  (β0, β1, β2)
# ==========================================================

lasso_kernel <- (function() {
  # Prepare data once (global inside closure)
  dat <- mtcars
  y   <- dat$mpg
  X   <- cbind(hp = dat$hp, wt = dat$wt)
  X   <- scale(X, center = TRUE, scale = TRUE)
  N   <- length(y)
  lambda <- 0.25
  
  # Return the actual kernel function: depends only on theta
  function(theta) {
    beta0 <- theta[1]
    beta  <- theta[-1]
    rss   <- sum((y - beta0 - as.vector(X %*% beta))^2)
    penalty <- lambda * sum(abs(beta))
    exp(-(rss / (2 * N) + penalty))
  }
})()

# ==========================================================
# Example: test kernel call
# ==========================================================
set.seed(123)
n_samples <- 10000
burn_in <- 2500
thin <- 1
theta_test <- c(25,-2, -3)
print(lasso_kernel(theta_test))  # should return a positive scalar

res_lasso = multivariate_gibbs_sample_ASG(lasso_kernel ,n_samples = n_samples, burn_in = burn_in, 
                                          theta_init = c(0,0,0), tol = 0.01, scale0 = 1)

res_lasso$time_taken
# Recreate the data exactly as in your lasso_kernel closure
dat <- mtcars
y   <- dat$mpg
X   <- cbind(hp = dat$hp, wt = dat$wt)
X   <- scale(X, center = TRUE, scale = TRUE)

lambda <- 0.25  # same penalty used in your kernel

# Fit glmnet at fixed lambda (same objective form: (1/2N)*RSS + lambda*||beta||_1)
fit_fixed <- glmnet(X, y, alpha = 1, intercept = TRUE, standardize = FALSE, lambda = lambda)
coef_fixed <- as.numeric(coef(fit_fixed, s = lambda))
beta0_glm  <- coef_fixed[1]
beta1_glm  <- coef_fixed[2]
beta2_glm  <- coef_fixed[3]

# (Optional) cross-validated choice, if you want it:
# cvfit <- cv.glmnet(X, y, alpha = 1, intercept = TRUE, standardize = FALSE)
# coef_cv <- as.numeric(coef(cvfit, s = "lambda.min"))

cat("GLMNET (fixed lambda) coefficients:\n",
    sprintf("beta0 = %.3f, beta1 = %.3f, beta2 = %.3f\n", beta0_glm, beta1_glm, beta2_glm))


estimated_samples = c(mean(res_lasso$samples[,1]), mean(res_lasso$samples[,2]), mean(res_lasso$samples[,3]))
estimated_samples
effectiveSize(res_lasso$samples)




#### Hinge Loss function ######


library(ggplot2)
set.seed(2)

# --- symmetric toy data ---
n <- 40
x <- seq(-2, 2, length.out = n)
y <- ifelse(x >= 0, 1, -1)
x <- scale(x, center = TRUE, scale = TRUE)

# Hinge kernel (θ = (β0, β1))
hinge_kernel_2d <- function(theta) {
  beta0 <- theta[1]
  beta1 <- theta[2]
  margins <- y * (beta0 + beta1 * as.numeric(x))
  loss <- sum(pmax(0, 1 - margins))
  exp(-loss)
}

set.seed(123)
n_samples <- 10000
burn_in <- 2500
thin <- 1

res_hinge = multivariate_gibbs_sample_ASG(hinge_kernel_2d ,n_samples = n_samples, burn_in = burn_in, 
                                          theta_init = c(0,0), tol = 0.01, scale0 = 1)
res_hinge$time_taken
effectiveSize(res_hinge$samples)
