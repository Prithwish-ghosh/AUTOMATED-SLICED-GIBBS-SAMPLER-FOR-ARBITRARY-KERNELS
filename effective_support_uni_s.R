library(nortest)
library(ggplot2)
library(dplyr)
library(mvtnorm)
library(doParallel)
library(gridExtra)
library(coda)

## ------------------------------------------------------------------
## Effective support estimator (speed-optimized: grid + interpolation
## replaces the nested optimize()/integrate() root-finding, without
## changing what is being solved for). Identical to the baseline file
## -- this fix is independent of the refresh-schedule question, it's
## what makes per-sweep recomputation affordable at all.
## ------------------------------------------------------------------
effective.support <- function(ker, tol = 0.01, scale0 = 1, subdivisions = 20,
                              grid_n = 400) {
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
    f.star <- function(u) k.star(u) / norm.const
    
    u_grid <- seq(0, 1, length.out = grid_n + 1)
    dens_grid <- f.star(u_grid)
    trapz <- (dens_grid[-1] + dens_grid[-length(dens_grid)]) / 2 * diff(u_grid)
    cdf_grid <- c(0, cumsum(trapz))
    cdf_grid <- cdf_grid / cdf_grid[length(cdf_grid)]
    
    F.star <- approxfun(u_grid, cdf_grid, rule = 2)
    Finv   <- approxfun(cdf_grid, u_grid, rule = 2)
    
    f <- function(x) k(x) / norm.const
    F <- function(x) F.star(pcauchy(x, scale = scale0))
    
    a.star <- Finv(0.5 * tol)
    b.star <- Finv(1 - 0.5 * tol)
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
    
    f <- function(x) k(x) / norm.const
    cdf <- function(x) {
      if (x <= a_test) return(0)
      if (x >= b_test) return(1)
      val <- integrate(f, a_test, x, subdivisions = subdivisions,
                       rel.tol = 1e-2, abs.tol = 1e-4)$value
      ifelse(is.finite(val), val, ifelse(x < a_test, 0, 1))
    }
    F <- Vectorize(cdf)
    
    g.lower <- function(x) F(x) - 0.5 * tol
    g.upper <- function(x) 1 - F(x) - 0.5 * tol
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

## ------------------------------------------------------------------
## Fixed-u 1D slice step
## ------------------------------------------------------------------
slice1d_fixed_u <- function(g, u, a, b, max_tries = 100000L) {
  for (t in seq_len(max_tries)) {
    x_new <- runif(1, a, b)
    gx <- g(x_new)
    if (is.finite(gx) && gx > u) return(x_new)
  }
  stop(sprintf(
    "slice1d_fixed_u: no acceptance in %d tries -- bracket [%.4g, %.4g] may not intersect the slice.",
    max_tries, a, b
  ))
}

## ------------------------------------------------------------------
## Univariate ASG (m = 1): unchanged from baseline. Bracket computed
## once, before sampling -- exact for m = 1 regardless of refresh
## policy, since g_i never depends on "other" coordinates.
## ------------------------------------------------------------------
dim1_gibbs_sample_ASG <- function(ker, n_samples = 10, burn_in = 2, thin = 1,
                                  theta_init, tol = 0.01, scale0 = 1,
                                  make_plots = TRUE) {
  start_time <- Sys.time()
  m <- 1
  total_iter <- burn_in + n_samples * thin
  theta_chain <- matrix(NA_real_, nrow = total_iter, ncol = m)
  colnames(theta_chain) <- paste0("theta_", seq_len(m))
  theta_new <- as.numeric(theta_init)
  
  bounds_list <- vector("list", m)
  f_list <- vector("list", m)
  for (i in seq_len(m)) {
    g_i <- (function(i) {
      function(xi) {
        tmp <- theta_new
        tmp[i] <- xi
        ker(tmp)
      }
    })(i)
    bnds <- effective.support(g_i, tol = tol / m, scale0 = scale0)
    bounds_list[[i]] <- c(a = bnds$lower, b = bnds$upper)
    f_list[[i]] <- bnds$f
  }
  
  for (t in seq_len(total_iter)) {
    u <- runif(1, 0, ker(theta_new))
    for (i in seq_len(m)) {
      g_i <- (function(i) {
        function(xi) {
          tmp <- theta_new
          tmp[i] <- xi
          ker(tmp)
        }
      })(i)
      a_i <- bounds_list[[i]]["a"]; b_i <- bounds_list[[i]]["b"]
      theta_new[i] <- slice1d_fixed_u(g_i, u, a = a_i, b = b_i)
    }
    theta_chain[t, ] <- theta_new
  }
  
  sample_idx <- seq(from = burn_in + 1L, to = total_iter, by = thin)
  theta_samples <- theta_chain[sample_idx, , drop = FALSE]
  end_time <- Sys.time()
  time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  if (make_plots) {
    for (i in seq_len(m)) {
      chain_i <- theta_chain[, i]
      samples_i <- theta_samples[, i]
      a <- bounds_list[[i]]["a"]; b <- bounds_list[[i]]["b"]
      f <- f_list[[i]]
      
      xs <- seq(a, b, length.out = 2000)
      df_target <- data.frame(x = xs, dens = f(xs))
      p_density <- ggplot() +
        geom_histogram(aes(x = samples_i, y = after_stat(density)), bins = 50, fill = "skyblue", color = "black", alpha = 0.7) +
        geom_line(data = df_target, aes(x = x, y = dens), linewidth = 2, color = "red") +
        geom_vline(xintercept = c(a, b), linetype = "dashed") +
        labs(title = paste("Dim", i, ": Samples Generated by ASG Desity plot"),
             x = paste0("theta_", i), y = "Density") +
        theme_minimal(base_size = 12) +
        annotate("text", x = min(samples_i, na.rm = TRUE), y = Inf,
                 label = sprintf("Time: %.2fs", time_taken),
                 hjust = 0, vjust = 1.2, size = 4)
      
      df_trace <- data.frame(iter = seq_along(chain_i), val = chain_i)
      p_trace <- ggplot(df_trace, aes(iter, val)) +
        geom_line(alpha = 0.8, size = 1.5) +
        geom_vline(xintercept = burn_in, linetype = "dashed", color = "red") +
        labs(title = paste("Dim", i, ": Trace"), x = "Iteration", y = paste0("theta_", i)) +
        theme_minimal(base_size = 12)
      
      acf_vals <- acf(chain_i, plot = FALSE, lag.max = 50)
      df_acf <- data.frame(lag = as.numeric(acf_vals$lag), acf = as.numeric(acf_vals$acf))
      p_acf <- ggplot(df_acf, aes(lag, acf)) +
        geom_hline(yintercept = 0) +
        geom_segment(aes(xend = lag, yend = 0)) +
        labs(title = paste("Dim", i, ": ACF"), x = "Lag", y = "ACF") +
        theme_minimal(base_size = 12)
      
      rmn <- cumsum(chain_i) / seq_along(chain_i)
      df_rm <- data.frame(iter = seq_along(chain_i), mean = rmn)
      p_rm <- ggplot(df_rm, aes(iter, mean)) +
        geom_line(size = 2) +
        geom_vline(xintercept = burn_in, linetype = "dashed", color = "red") +
        labs(title = paste("Dim", i, ": Running mean"), x = "Iteration", y = "Mean") +
        theme_minimal(base_size = 12)
      
      print(p_density); print(p_trace); print(p_acf); print(p_rm)
    }
  }
  list(samples = theta_samples, chain = theta_chain,
       burn_in = burn_in, thin = thin, time_taken = time_taken)
}

## ------------------------------------------------------------------
## Multivariate ASG -- PER-SWEEP variant
## Bracket [a_k, b_k] is recomputed at EVERY sweep t, for EVERY
## coordinate k -- matches the original per-Gibbs-cycle theory:
##   [a_k, b_k] <- EFFECTIVESUPPORT1D(g_k(.), eps/m, s0; warm-start)
## invoked inside BOTH the "for t" and "for k" loops. tol/m applies
## Remark 2.1's Bonferroni correction. No repeat/Corollary-2.2 branch
## needed: this matches the "Approximate draws" claim via Proposition
## 2.1's per-sweep coverage guarantee (P(S_u subset [a,b]) >= 1-eps/m
## at every t), not via an exactness correction.
## ------------------------------------------------------------------
multivariate_gibbs_sample_ASG <- function(ker, n_samples = 10, burn_in = 2, thin = 1,
                                          theta_init, tol = 0.01, scale0 = 1,
                                          make_plots = TRUE) {
  start_time <- Sys.time()
  m <- length(theta_init)
  total_iter <- burn_in + n_samples * thin
  theta_chain <- matrix(NA_real_, nrow = total_iter, ncol = m)
  colnames(theta_chain) <- paste0("theta_", seq_len(m))
  theta_new <- as.numeric(theta_init)
  
  bounds_list <- vector("list", m)   # kept for diagnostics only; updated every sweep below
  
  ## Sampling
  for (t in seq_len(total_iter)) {
    u <- runif(1, 0, ker(theta_new))
    for (i in seq_len(m)) {
      g_i <- (function(i) {
        function(xi) {
          tmp <- theta_new
          tmp[i] <- xi
          ker(tmp)
        }
      })(i)
      
      ## Recompute the bracket every sweep, every coordinate -- matches
      ## the per-Gibbs-cycle theory. tol/m applies the Bonferroni
      ## correction of Remark 2.1.
      bnds_i <- effective.support(g_i, tol = tol / m, scale0 = scale0)
      a_i <- bnds_i$lower; b_i <- bnds_i$upper
      bounds_list[[i]] <- c(a = a_i, b = b_i)
      
      theta_new[i] <- slice1d_fixed_u(g_i, u, a = a_i, b = b_i)
    }
    theta_chain[t, ] <- theta_new
  }
  
  sample_idx <- seq(from = burn_in + 1L, to = total_iter, by = thin)
  theta_samples <- theta_chain[sample_idx, , drop = FALSE]
  end_time <- Sys.time()
  time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  if (make_plots) {
    diag_series <- theta_chain[, 1]
    series_name <- colnames(theta_chain)[1]
    
    df_trace <- data.frame(iter = seq_along(diag_series), val = diag_series)
    p_trace <- ggplot(df_trace, aes(x = iter, y = val)) +
      geom_line(alpha = 0.8, size = 1.5) +
      geom_vline(xintercept = burn_in, linetype = "dashed", color = "red") +
      labs(title = "Trace for theta (first coordinate)", x = "Iteration", y = series_name) +
      theme_minimal(base_size = 12) +
      coord_cartesian(ylim = c(-8, 8))
    
    acf_vals <- acf(diag_series, plot = FALSE, lag.max = 50)
    df_acf <- data.frame(lag = as.numeric(acf_vals$lag), acf = as.numeric(acf_vals$acf))
    p_acf <- ggplot(df_acf, aes(x = lag, y = acf)) +
      geom_hline(yintercept = 0, size = 1) +
      geom_segment(aes(xend = lag, yend = 0)) +
      labs(title = "ACF for theta (first coordinate)", x = "Lag", y = "ACF") +
      theme_minimal(base_size = 12)
    
    running_mean <- cumsum(diag_series) / seq_along(diag_series)
    df_rm <- data.frame(iter = seq_along(diag_series), mean = running_mean)
    p_rm <- ggplot(df_rm, aes(x = iter, y = mean)) +
      geom_line(size = 2) +
      geom_vline(xintercept = burn_in, linetype = "dashed", color = "red") +
      labs(title = "Running mean for theta (first coordinate)", x = "Iteration", y = "Running mean") +
      theme_minimal(base_size = 12)
    
    print(p_trace); print(p_acf); print(p_rm)
    
    tiny <- 1e-300
    logK_series <- apply(theta_chain, 1, function(th) {
      val <- ker(as.numeric(th))
      log(pmax(val, tiny))
    })
    
    df_logk_trace <- data.frame(iter = seq_along(logK_series), logK = logK_series)
    p_logk_trace <- ggplot(df_logk_trace, aes(x = iter, y = logK)) +
      geom_line(alpha = 0.8, size = 2) +
      geom_vline(xintercept = burn_in, linetype = "dashed", color = "red") +
      labs(title = "Trace of log-kernel: log K(theta^{(t)})", x = "Iteration t", y = "log K(theta)") +
      theme_minimal(base_size = 12)
    
    acf_logk <- acf(logK_series[(burn_in + 1L):length(logK_series)], plot = FALSE, lag.max = 50)
    df_logk_acf <- data.frame(lag = as.numeric(acf_logk$lag), acf = as.numeric(acf_logk$acf))
    p_logk_acf <- ggplot(df_logk_acf, aes(x = lag, y = acf)) +
      geom_hline(yintercept = 0, size = 2) +
      geom_segment(aes(xend = lag, yend = 0)) +
      labs(title = "ACF of log-kernel: log K(theta^{(t)})", x = "Lag", y = "ACF") +
      theme_minimal(base_size = 12)
    
    print(p_logk_trace); print(p_logk_acf)
  }
  
  list(
    samples    = theta_samples,
    chain      = theta_chain,
    burn_in    = burn_in,
    thin       = thin,
    time_taken = time_taken
  )
}
