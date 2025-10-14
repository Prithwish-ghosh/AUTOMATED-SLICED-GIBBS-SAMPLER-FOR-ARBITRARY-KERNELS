
library(ggplot2)
library(RColorBrewer)

banana_kernel <- function(theta) {
  x1 <- theta[1]
  x2 <- theta[2]
  log_pi <- -x1^2 / 10 - x2^2 / 10 - 2 * (x2 - x1^2)^2
  exp(log_pi)
}

theta <- c(0, 0)
kernel_value <- banana_kernel(theta)
cat("Kernel value at theta =", theta, ":", kernel_value, "\n")

x1_grid <- seq(-4, 4, length.out = 100)
x2_grid <- seq(-2, 6, length.out = 100)

z <- matrix(NA_real_, nrow = length(x1_grid), ncol = length(x2_grid))
for (i in seq_along(x1_grid)) {
  for (j in seq_along(x2_grid)) {
    z[i, j] <- banana_kernel(c(x1_grid[i], x2_grid[j]))
  }
}

grid_data <- expand.grid(x1 = x1_grid, x2 = x2_grid)
grid_data$z <- as.vector(z)

p <- ggplot(grid_data, aes(x = x1, y = x2, z = z)) +
  geom_contour_filled(bins = 20, aes(fill = after_stat(level))) +
  scale_fill_brewer(palette = "RdBu", direction = -1) +
  geom_point(aes(x = 0, y = 0), color = "green", size = 3) +
  labs(title = "Contour Plot of Banana-Shaped Kernel", x = "x1", y = "x2") +
  theme_minimal()

source("effective_support_uni_s.R") 

kernel_x2_given_x1_factory <- function(x1_fixed) {
  function(x2) banana_kernel(c(x1_fixed, x2))
}
kernel_x1_given_x2_factory <- function(x2_fixed) {
  function(x1) banana_kernel(c(x1, x2_fixed))
}

support_x2_given_x1 <- function(x1_fixed, tol = 0.01, scale0 = 1) {
  eff <- effective.support(kernel_x2_given_x1_factory(x1_fixed), tol = tol, scale0 = scale0)
  c(lower = eff$lower, upper = eff$upper)
}
support_x1_given_x2 <- function(x2_fixed, tol = 0.01, scale0 = 1) {
  eff <- effective.support(kernel_x1_given_x2_factory(x2_fixed), tol = tol, scale0 = scale0)
  c(lower = eff$lower, upper = eff$upper)
}

## ===============================
## Choose the common slice set c
## c = {coordinate at global mode, -2, 1}
## ===============================

## Find global mode (argmax) numerically (should be ~ (0,0))
neg_log <- function(th) {  # minimize negative of log kernel
  x1 <- th[1]; x2 <- th[2]
  # -log kernel up to additive constant
  (x1^2)/10 + (x2^2)/10 + 2*(x2 - x1^2)^2
}
opt <- optim(c(0,0), neg_log, method = "BFGS")
x1_mode <- opt$par[1]
x2_mode <- opt$par[2]
cat(sprintf("Estimated mode at (x1, x2) ≈ (%.4f, %.4f)\n", x1_mode, x2_mode))

## Use the same set for both slice directions:
c_values <- unique(c(signif(x1_mode, 6), -2, 1))
## (x1_mode ≈ x2_mode here, so using x1_mode is fine; it equals 0)

## ===============================
## Build overlay data (segments + mode curve)
## ===============================
x1_min <- min(x1_grid); x1_max <- max(x1_grid)
x2_min <- min(x2_grid); x2_max <- max(x2_grid)

## Vertical segments: bounds of X2 | X1 = c, for c in c_values
vertical_df <- do.call(rbind, lapply(c_values, function(x1c) {
  bds <- support_x2_given_x1(x1c, tol = 0.01, scale0 = 1)
  data.frame(x1 = x1c, y0 = bds[1], y1 = bds[2])
}))
## Clip to plotting window
vertical_df$y0 <- pmax(vertical_df$y0, x2_min)
vertical_df$y1 <- pmin(vertical_df$y1, x2_max)
vertical_df <- vertical_df[vertical_df$y1 > vertical_df$y0, , drop = FALSE]

## Horizontal segments: bounds of X1 | X2 = c, for c in c_values
horizontal_df <- do.call(rbind, lapply(c_values, function(x2c) {
  bds <- support_x1_given_x2(x2c, tol = 0.01, scale0 = 1)
  data.frame(y2 = x2c, x0 = bds[1], x1 = bds[2])
}))
## Clip to plotting window
horizontal_df$x0 <- pmax(horizontal_df$x0, x1_min)
horizontal_df$x1 <- pmin(horizontal_df$x1, x1_max)
horizontal_df <- horizontal_df[horizontal_df$x1 > horizontal_df$x0, , drop = FALSE]

## Conditional mode curve x2 = x1^2 (reference)
mode_df <- data.frame(x1 = seq(x1_min, x1_max, length.out = 400))
mode_df$x2 <- mode_df$x1^2

## ===============================
## Add overlays to your contour
## ===============================
p <- p +
  ## vertical effective-support segments (blue) for X2 | X1 = c
  geom_segment(data = vertical_df,
               aes(x = x1, xend = x1, y = y0, yend = y1),
               inherit.aes = FALSE, linewidth = 0.8, color = "blue") +
  ## horizontal effective-support segments (red) for X1 | X2 = c
  geom_segment(data = horizontal_df,
               aes(x = x0, xend = x1, y = y2, yend = y2),
               inherit.aes = FALSE, linewidth = 0.8, color = "red") +
  ## labels
  annotate("text", x = x1_min + 0.3, y = x2_max - 0.5,
           label = "x2 | x1 = c", color = "blue", size = 5, hjust = 0) +
  annotate("text", x = x1_min + 0.3, y = x2_max - 1.1,
           label = "x1 | x2 = c", color = "red",  size = 5, hjust = 0)

## Show plot
print(p)

## (Optional) Inspect the bounds printed:
cat("\nVertical bounds: X2 | X1 = c\n")
print(vertical_df)
cat("\nHorizontal bounds: X1 | X2 = c\n")
print(horizontal_df)
