# Define the banana kernel (bivariate)
banana_kernel <- function(theta) {
  x1 <- theta[1]
  x2 <- theta[2]
  log_pi <- -x1^2 / 10 - x2^2 / 10 - 2 * (x2 - x1^2)^2
  exp(log_pi)
}

# Provided effective.support function
effective.support <- function(ker, tol = 0.01, scale0 = 1, subdivisions = 20) {
  k <- Vectorize(ker)
  
  k.star <- function(u) {
    x <- qcauchy(u, scale = scale0)
    dens <- k(x)
    cauchy_dens <- dcauchy(x, scale = scale0)
    ifelse(cauchy_dens > 1e-10 & dens >= 0 & is.finite(dens), dens / cauchy_dens, 0)
  }
  
  norm.const <- tryCatch(
    integrate(k.star, 0, 1, subdivisions = subdivisions, rel.tol = 1e-2, abs.tol = 1e-4)$value,
    error = function(e) NA
  )
  
  if (is.finite(norm.const) && norm.const > 0) {
    f.star <- function(u) { k.star(u) / norm.const }
    cdf <- function(u) {
      if (u <= 0) return(0)
      if (u >= 1) return(1)
      val <- tryCatch(
        integrate(f.star, lower = 0, upper = u, subdivisions = subdivisions,
                  rel.tol = 1e-2, abs.tol = 1e-4)$value,
        error = function(e) NA
      )
      ifelse(is.finite(val), val, NA)
    }
    F.star <- Vectorize(cdf)
    f <- function(x) { k(x) / norm.const }
    F <- function(x) {
      u <- pcauchy(x, scale = scale0)
      val <- F.star(u)
      ifelse(is.finite(val), val, ifelse(x <= qcauchy(0, scale0), 0, 1))
    }
    
    g.lower <- function(u) { F.star(u) - 0.5 * tol }
    g.upper <- function(u) { 1 - F.star(u) - 0.5 * tol }
    a.star <- tryCatch(
      optimize(function(u) abs(g.lower(u)), interval = c(0, 0.5), tol = 1e-4)$minimum,
      error = function(e) 0.001
    )
    b.star <- tryCatch(
      optimize(function(u) abs(g.upper(u)), interval = c(a.star, 1), tol = 1e-4)$minimum,
      error = function(e) 0.999
    )
    a <- qcauchy(a.star, scale = scale0)
    b <- qcauchy(b.star, scale = scale0)
  } else {
    x_test <- seq(-100, 100, length.out = 1000)
    k_vals <- k(x_test)
    support_idx <- which(k_vals > 1e-10)
    if (length(support_idx) == 0) stop("No non-zero support detected for kernel.")
    
    a_approx <- min(x_test[support_idx])
    b_approx <- max(x_test[support_idx])
    a_test <- a_approx
    b_test <- b_approx
    
    norm.const <- integrate(k, a_test, b_test, subdivisions = subdivisions,
                            rel.tol = 1e-2, abs.tol = 1e-4)$value
    if (!is.finite(norm.const) || norm.const <= 0) {
      stop("Normalizing constant in fallback method is non-finite or zero.")
    }
    
    f <- function(x) { k(x) / norm.const }
    cdf <- function(x) {
      if (x <= a_test) return(0)
      if (x >= b_test) return(1)
      val <- integrate(f, a_test, x, subdivisions = subdivisions,
                       rel.tol = 1e-2, abs.tol = 1e-4)$value
      ifelse(is.finite(val), val, ifelse(x < a_test, 0, 1))
    }
    F <- Vectorize(cdf)
    
    g.lower <- function(x) { F(x) - 0.5 * tol }
    g.upper <- function(x) { 1 - F(x) - 0.5 * tol }
    a <- tryCatch(
      optimize(function(x) abs(g.lower(x)), interval = c(a_test, b_test), tol = 1e-2)$minimum,
      error = function(e) a_approx
    )
    b <- tryCatch(
      optimize(function(x) abs(g.upper(x)), interval = c(a, b_test), tol = 1e-2)$minimum,
      error = function(e) b_approx
    )
  }
  
  attr(f, "norm.const") <- norm.const
  attr(f, "lower") <- a
  attr(f, "upper") <- b
  
  return(list(lower = a, upper = b, f = f, norm.const = norm.const))
}

# Profile the kernel for x1 (fix x2 at conditional mode or a reasonable value)
# Conditional mode of x2 given x1: x2 = x1^2 (from minimizing (x2 - x1^2)^2)
kernel_x1 <- function(x) {
  x2_opt <- x^2  # Mode of x2 given x1
  banana_kernel(c(x, x2_opt))
}

# Profile the kernel for x2 (fix x1 at 0, near the mode)
kernel_x2 <- function(x) {
  banana_kernel(c(0, x))
}

# Compute effective support for x1 and x2
x1_support <- effective.support(kernel_x1, tol = 0.01, scale0 = 1, subdivisions = 100)
x2_support <- effective.support(kernel_x2, tol = 0.01, scale0 = 1, subdivisions = 100)

x1_support
x2_support



# Function to find conditional max/min of x2 given x1
find_x2_bounds <- function(x1, lower = -10, upper = 10) {
  # Maximize kernel for x2 given x1
  max_x2 <- optimize(
    f = function(x2) banana_kernel(c(x1, x2)),
    interval = c(lower, upper),
    maximum = TRUE
  )$maximum
  
  # Find min by checking boundaries (where density is low but non-zero)
  dens_lower <- banana_kernel(c(x1, lower))
  dens_upper <- banana_kernel(c(x1, upper))
  min_x2 <- if (dens_lower > dens_upper) lower else upper
  
  list(max = max_x2, min = min_x2)
}

# Function to find conditional max/min of x1 given x2
find_x1_bounds <- function(x2, lower = -10, upper = 10) {
  # Maximize kernel for x1 given x2
  max_x1 <- optimize(
    f = function(x1) banana_kernel(c(x1, x2)),
    interval = c(lower, upper),
    maximum = TRUE
  )$maximum
  
  # Find min by checking boundaries
  dens_lower <- banana_kernel(c(lower, x2))
  dens_upper <- banana_kernel(c(upper, x2))
  min_x1 <- if (dens_lower > dens_upper) lower else upper
  
  list(max = max_x1, min = min_x1)
}

# Profile kernels for x1 and x2 using max and min values
# Case 1: Fix at max (mode-like)
kernel_x1_max <- function(x) {
  x2_opt <- find_x2_bounds(x)$max  # Max of x2 given x1
  banana_kernel(c(x, x2_opt))
}
kernel_x2_max <- function(x) {
  x1_opt <- find_x1_bounds(x)$max  # Max of x1 given x2
  banana_kernel(c(x1_opt, x))
}


# Compute effective support for x1 and x2
x1_support_max <- effective.support(kernel_x1_max, tol = 0.01, scale0 = 1, subdivisions = 100)
x2_support_max <- effective.support(kernel_x2_max, tol = 0.01, scale0 = 1, subdivisions = 100)

# Define the banana kernel
banana_kernel <- function(theta) {
  x1 <- theta[1]
  x2 <- theta[2]
  log_pi <- -x1^2 / 10 - x2^2 / 10 - 2 * (x2 - x1^2)^2
  exp(log_pi)
}

# Vectorized kernel functions
kernel_x1 <- Vectorize(function(x) {
  x2_opt <- x^2  # Mode of x2 given x1
  banana_kernel(c(x, x2_opt))
})

kernel_x2 <- Vectorize(function(x) {
  banana_kernel(c(0, x))
})

find_x2_bounds <- function(x1, lower = -10, upper = 10) {
  max_x2 <- optimize(
    f = function(x2) banana_kernel(c(x1, x2)),
    interval = c(lower, upper),
    maximum = TRUE
  )$maximum
  dens_lower <- banana_kernel(c(x1, lower))
  dens_upper <- banana_kernel(c(x1, upper))
  min_x2 <- if (dens_lower > dens_upper) lower else upper
  list(max = max_x2, min = min_x2)
}

find_x1_bounds <- function(x2, lower = -10, upper = 10) {
  max_x1 <- optimize(
    f = function(x1) banana_kernel(c(x1, x2)),
    interval = c(lower, upper),
    maximum = TRUE
  )$maximum
  dens_lower <- banana_kernel(c(lower, x2))
  dens_upper <- banana_kernel(c(upper, x2))
  min_x1 <- if (dens_lower > dens_upper) lower else upper
  list(max = max_x1, min = min_x1)
}

kernel_x1_max <- Vectorize(function(x) {
  x2_opt <- find_x2_bounds(x)$max
  banana_kernel(c(x, x2_opt))
})

kernel_x2_max <- Vectorize(function(x) {
  x1_opt <- find_x1_bounds(x)$max
  banana_kernel(c(x1_opt, x))
})

# Compute effective support
x1_support <- effective.support(kernel_x1, tol = 0.01, scale0 = 1, subdivisions = 100)
x2_support <- effective.support(kernel_x2, tol = 0.01, scale0 = 1, subdivisions = 100)
x1_support_max <- effective.support(kernel_x1_max, tol = 0.01, scale0 = 1, subdivisions = 100)
x2_support_max <- effective.support(kernel_x2_max, tol = 0.01, scale0 = 1, subdivisions = 100)


# Additional cases with random fixed values
random_fixed_values <- c(1, -2)

# For x1 profiling (fix x2 at random values)
kernel_x1_random1 <- Vectorize(function(x) {
  banana_kernel(c(x, random_fixed_values[1]))
})
kernel_x1_random2 <- Vectorize(function(x) {
  banana_kernel(c(x, random_fixed_values[2]))
})

# For x2 profiling (fix x1 at random values)
kernel_x2_random1 <- Vectorize(function(x) {
  banana_kernel(c(random_fixed_values[1], x))
})
kernel_x2_random2 <- Vectorize(function(x) {
  banana_kernel(c(random_fixed_values[2], x))
})

# Compute effective support for random cases
x1_support_random1 <- effective.support(kernel_x1_random1, tol = 0.01, scale0 = 1, subdivisions = 100)
x1_support_random2 <- effective.support(kernel_x1_random2, tol = 0.01, scale0 = 1, subdivisions = 100)
x2_support_random1 <- effective.support(kernel_x2_random1, tol = 0.01, scale0 = 1, subdivisions = 100)
x2_support_random2 <- effective.support(kernel_x2_random2, tol = 0.01, scale0 = 1, subdivisions = 100)

# Evaluate density for histograms (use wider grid to cover all supports)
x1_grid <- seq(min(c(x1_support$lower, x1_support_max$lower, x1_support_random1$lower, x1_support_random2$lower)) - 0.5,
               max(c(x1_support$upper, x1_support_max$upper, x1_support_random1$upper, x1_support_random2$upper)) + 0.5, length.out = 100)
x2_grid <- seq(min(c(x2_support$lower, x2_support_max$lower, x2_support_random1$lower, x2_support_random2$lower)) - 0.5,
               max(c(x2_support$upper, x2_support_max$upper, x2_support_random1$upper, x2_support_random2$upper)) + 0.5, length.out = 100)

# Normalize kernels
norm_x1 <- integrate(kernel_x1, x1_support$lower, x1_support$upper, subdivisions = 100)$value
norm_x1_max <- integrate(kernel_x1_max, x1_support_max$lower, x1_support_max$upper, subdivisions = 100)$value
norm_x1_random1 <- integrate(kernel_x1_random1, x1_support_random1$lower, x1_support_random1$upper, subdivisions = 100)$value
norm_x1_random2 <- integrate(kernel_x1_random2, x1_support_random2$lower, x1_support_random2$upper, subdivisions = 100)$value

norm_x2 <- integrate(kernel_x2, x2_support$lower, x2_support$upper, subdivisions = 100)$value
norm_x2_max <- integrate(kernel_x2_max, x2_support_max$lower, x2_support_max$upper, subdivisions = 100)$value
norm_x2_random1 <- integrate(kernel_x2_random1, x2_support_random1$lower, x2_support_random1$upper, subdivisions = 100)$value
norm_x2_random2 <- integrate(kernel_x2_random2, x2_support_random2$lower, x2_support_random2$upper, subdivisions = 100)$value

# Compute densities
x1_density <- kernel_x1(x1_grid) / norm_x1
x1_density_max <- kernel_x1_max(x1_grid) / norm_x1_max
x1_density_random1 <- kernel_x1_random1(x1_grid) / norm_x1_random1
x1_density_random2 <- kernel_x1_random2(x1_grid) / norm_x1_random2

x2_density <- kernel_x2(x2_grid) / norm_x2
x2_density_max <- kernel_x2_max(x2_grid) / norm_x2_max
x2_density_random1 <- kernel_x2_random1(x2_grid) / norm_x2_random1
x2_density_random2 <- kernel_x2_random2(x2_grid) / norm_x2_random2


# Compute true marginal densities
library(cubature)  # For numerical integration

# Marginal density for x1 (integrate over x2)
marginal_x1 <- Vectorize(function(x1) {
  integrate(
    Vectorize(function(x2) banana_kernel(c(x1, x2))),
    lower = -10, upper = 10
  )$value
})

# Marginal density for x2 (integrate over x1)
marginal_x2 <- Vectorize(function(x2) {
  integrate(
    Vectorize(function(x1) banana_kernel(c(x1, x2))),
    lower = -10, upper = 10
  )$value
})

# Evaluate marginal density over grid
x1_grid <- seq(min(c(x1_support$lower, x1_support_max$lower, x1_support_random1$lower, x1_support_random2$lower)) - 0.5,
               max(c(x1_support$upper, x1_support_max$upper, x1_support_random1$upper, x1_support_random2$upper)) + 0.5, length.out = 100)
x2_grid <- seq(min(c(x2_support$lower, x2_support_max$lower, x2_support_random1$lower, x2_support_random2$lower)) - 0.5,
               max(c(x2_support$upper, x2_support_max$upper, x2_support_random1$upper, x2_support_random2$upper)) + 0.5, length.out = 100)

# Normalize marginal densities
norm_x1 <- integrate(marginal_x1, min(x1_grid), max(x1_grid), subdivisions = 100)$value
norm_x2 <- integrate(marginal_x2, min(x2_grid), max(x2_grid), subdivisions = 100)$value

x1_density <- marginal_x1(x1_grid) / norm_x1
x2_density <- marginal_x2(x2_grid) / norm_x2

# Plot histograms
par(mfrow = c(1, 1))

# x1 density plot
plot(x1_grid, x1_density, type = "l", col = "black", lwd = 2,
     main = "Marginal Density of x1 with Support Bounds", xlab = "x1", ylab = "Density",
     ylim = c(0, max(x1_density) * 1.1))
abline(v = x1_support$lower, col = "blue", lty = 2, lwd = 2)
abline(v = x1_support$upper, col = "blue", lty = 2, lwd = 2)
abline(v = x1_support_max$lower, col = "red", lty = 1, lwd = 2)
abline(v = x1_support_max$upper, col = "red", lty = 1, lwd = 2)
abline(v = x1_support_random1$lower, col = "green", lty = 4, lwd = 2)
abline(v = x1_support_random1$upper, col = "green", lty = 4, lwd = 2)
abline(v = x1_support_random2$lower, col = "purple", lty = 5, lwd = 2)
abline(v = x1_support_random2$upper, col = "purple", lty = 5, lwd = 2)
legend("topright", legend = c("Mode Bounds", "Max Bounds", "Random1 Bounds", "Random2 Bounds"),
       col = c("blue", "red", "green", "purple"), lty = c(2, 1, 4, 5), lwd = 2,
       cex = 0.7, box.lwd = 1, bg = "white")

# x2 density plot
plot(x2_grid, x2_density, type = "l", col = "black", lwd = 2,
     main = "Marginal Density of x2 with Support Bounds", xlab = "x2", ylab = "Density",
     ylim = c(0, max(x2_density) * 1.1))
abline(v = x2_support$lower, col = "blue", lty = 2, lwd = 2)
abline(v = x2_support$upper, col = "blue", lty = 2, lwd = 2)
abline(v = x2_support_max$lower, col = "red", lty = 1, lwd = 2)
abline(v = x2_support_max$upper, col = "red", lty = 1, lwd = 2)
abline(v = x2_support_random1$lower, col = "green", lty = 4, lwd = 2)
abline(v = x2_support_random1$upper, col = "green", lty = 4, lwd = 2)
abline(v = x2_support_random2$lower, col = "purple", lty = 5, lwd = 2)
abline(v = x2_support_random2$upper, col = "purple", lty = 5, lwd = 2)
legend("topright", legend = c("Mode Bounds", "Max Bounds", "Random1 Bounds", "Random2 Bounds"),
       col = c("blue", "red", "green", "purple"), lty = c(2, 1, 4, 5), lwd = 2,
       cex = 0.7, box.lwd = 1, bg = "white")


