library(nortest)
library(ggplot2)
library(dplyr)
library(mvtnorm)
library(doParallel)
library(gridExtra)
library(coda)

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



# Univariate slice sampling function (modified to return bounds and normalized density)
univariate_slice_sample <- function(ker, n_samples, burn_in, thin, init, tol = 0.1, scale0 = 1) {
  
  bounds <- effective.support(ker)
  
  a <- bounds$lower
  b <- bounds$upper
  f <- bounds$f
  k <- function(x) f(x) * bounds$norm.const
  
  x <- init
  total_iter <- burn_in + n_samples * thin
  chain <- numeric(total_iter)
  
  for (i in 1:total_iter) {
    y <- runif(1, 0, k(x))
    accept <- FALSE
    while (!accept) {
      x_new <- runif(1, a, b)
      if (k(x_new) > y) {
        x <- x_new
        accept <- TRUE
      }
    }
    chain[i] <- x
  }
  
  sample_indices <- seq(burn_in + 1, total_iter, by = thin)
  samples <- chain[sample_indices]
  
  return(list(samples = samples, chain = chain, lower = a, upper = b, f = f))
}

