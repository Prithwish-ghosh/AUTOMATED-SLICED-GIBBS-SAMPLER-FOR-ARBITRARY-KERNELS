## ------------------------------------------------------------------
## Multivariate ASG  (kernel K on its natural scale -- no log kernel)
##
## Requires effective.support() from effective_support_uni_s.R.
## slice1d_fixed_u() is NOT needed from that file any more: the
## version used here is slice1d_fixed_u_x0() below.
##
## [MODIFIED] Optional supervisor check condition: check_bracket
##   check_bracket = TRUE  -> bracket [a, b] = smallest and largest
##       roots of K = u, found on the Cauchy scale v = F_C(x; scale0)
##       by a fixed grid scan + bisection. Each endpoint satisfies
##         K(a_out) <= u < K(a_in),  |a_in - a_out| < e   (same for b)
##       and every grid point outside [a, b] has K <= u, so [a, b]
##       contains ALL components of the slice S_u (under assumption
##       A3: every slice component is wider than the grid spacing on
##       the Cauchy scale). Outside ends are returned, and x0 is used
##       only for the check a < x0 < b, never to anchor the search.
##   check_bracket = FALSE -> original behaviour: bracket = effective-
##       support box from effective.support() (no check).
## In both cases x0 is RETAINED if it falls outside [a, b] or the
## rejection step runs out of tries (Tierney retain rule), which keeps
## the update exactly invariant.
## ------------------------------------------------------------------


## [NEW] Supervisor's check condition: extreme roots of K = u.
## g : conditional kernel (natural scale), u : slice height,
## x0: current value (only checked, never used to search).
extreme_roots_bracket <- function(g, u, x0, scale0 = 1, n_grid = 1000L,
                                  e = 1e-8, max_refine = 40L) {
  inside <- function(v) { val <- g(qcauchy(v, scale = scale0)); isTRUE(val > u) }
  vs <- seq_len(n_grid) / (n_grid + 1)          # fixed grid on (0, 1)
  
  ## bisection keeping inside(v_in) = TRUE, inside(v_out) = FALSE;
  ## returns the OUTSIDE end, so no slice mass is cut off
  bisect_out <- function(v_out, v_in) {
    while (abs(v_in - v_out) >= e) {
      mid <- (v_in + v_out) / 2
      if (inside(mid)) v_in <- mid else v_out <- mid
    }
    v_out
  }
  ## slice reaches beyond the outermost grid point: refine toward edge
  refine_edge <- function(v_in, edge) {
    for (j in seq_len(max_refine)) {
      v_try <- edge + (v_in - edge) / 2
      if (!inside(v_try)) return(list(v = bisect_out(v_try, v_in), capped = FALSE))
      v_in <- v_try
    }
    list(v = v_in, capped = TRUE)
  }
  
  kL <- NA_integer_                              # first grid point in slice
  for (k in seq_along(vs)) if (inside(vs[k])) { kL <- k; break }
  if (is.na(kL)) return(list(lower = NA_real_, upper = NA_real_, ok = FALSE))
  kR <- kL                                       # last grid point in slice
  for (k in rev(seq_along(vs))) if (inside(vs[k])) { kR <- k; break }
  
  capped <- FALSE
  if (kL == 1L) { r <- refine_edge(vs[1], 0); a_v <- r$v; capped <- r$capped
  } else a_v <- bisect_out(vs[kL - 1L], vs[kL])
  if (kR == n_grid) { r <- refine_edge(vs[n_grid], 1); b_v <- r$v; capped <- capped || r$capped
  } else b_v <- bisect_out(vs[kR + 1L], vs[kR])
  
  a <- qcauchy(a_v, scale = scale0); b <- qcauchy(b_v, scale = scale0)
  list(lower = a, upper = b,
       ok = is.finite(a) && is.finite(b) && (a < x0) && (x0 < b) && !capped)
}


## [NEW] Fixed-u 1D slice step with the retain rule (replaces
## slice1d_fixed_u from effective_support_uni_s.R). Draws uniformly
## in x on [a, b] until K > u. Keeps x0 if x0 is outside [a, b] or if
## max_tries is exhausted (previously returned NULL and crashed).
slice1d_fixed_u_x0 <- function(g, u, a, b, x0, max_tries = 100000L) {
  if (!is.finite(a) || !is.finite(b) || x0 < a || x0 > b)
    return(list(x = x0, retained = TRUE))
  for (t in seq_len(max_tries)) {
    x_new <- runif(1, a, b)
    gx <- g(x_new)
    if (is.finite(gx) && gx > u) return(list(x = x_new, retained = FALSE))
  }
  list(x = x0, retained = TRUE)
}


## [NEW] One coordinate update, with or without the check.
asg_coordinate_update <- function(g, u, x0, check_bracket, box = NULL,
                                  tol = 0.01, scale0 = 1,
                                  n_grid = 1000L, bisect_e = 1e-8) {
  if (check_bracket) {
    br <- extreme_roots_bracket(g, u, x0, scale0 = scale0,
                                n_grid = n_grid, e = bisect_e)
    if (!br$ok) return(list(x = x0, retained = TRUE))
    a <- br$lower; b <- br$upper
  } else {
    if (is.null(box)) box <- effective.support(g, tol = tol, scale0 = scale0)
    a <- box$lower; b <- box$upper
  }
  slice1d_fixed_u_x0(g, u, a, b, x0)
}



dim1_gibbs_sample_ASG <- function(ker, n_samples = 10 , burn_in = 2, thin = 1,
                                  theta_init, tol = 0.01, scale0 = 1,
                                  bisect_e = 1e-8,
                                  check_bracket = TRUE,   # [MODIFIED] option
                                  n_grid = 1000L,         # [MODIFIED] grid size for check
                                  make_plots = TRUE) {
  start_time <- Sys.time()
  m <- 1
  n_retained <- 0L                                        # [MODIFIED]
  total_iter <- burn_in + n_samples * thin
  theta_chain <- matrix(NA_real_, nrow = total_iter, ncol = m)
  colnames(theta_chain) <- paste0("theta_", seq_len(m))
  theta_new <- as.numeric(theta_init)
  
  ## [MODIFIED comment] effective.support() is called ONCE here (m = 1):
  ## its box is the sampling bracket when check_bracket = FALSE, and its
  ## density f is the overlay in the diagnostic plots in both cases.
  bounds_list <- vector("list", m)   # box (used when check_bracket = FALSE)
  f_list <- vector("list", m)        # density overlay for the plots
  for (i in seq_len(m)) {
    g_i <- (function(i) {
      function(xi) {
        tmp <- theta_new
        tmp[i] <- xi
        ker(tmp)
      }
    })(i)
    bnds <- effective.support(g_i, tol = tol/m, scale0 = scale0)
    bounds_list[[i]] <- c(a = bnds$lower, b = bnds$upper)
    f_list[[i]] <- bnds$f         # normalized density over R (or approx support)
  }
  
  ## Sampling
  for (t in seq_len(total_iter)) {
    ## ONE u for the full sweep
    u <- runif(1, 0, ker(theta_new))
    for (i in seq_len(m)) {
      g_i <- (function(i) {
        function(xi) {
          tmp <- theta_new
          tmp[i] <- xi
          ker(tmp)
        }
      })(i)
      
      ## ---- [MODIFIED] check_bracket = TRUE: extreme roots of K = u
      ## (supervisor's check). check_bracket = FALSE: the fixed
      ## effective-support box computed once above (m = 1, so the
      ## conditional kernel never changes). ----
      box_i <- list(lower = unname(bounds_list[[i]]["a"]),
                    upper = unname(bounds_list[[i]]["b"]))
      upd <- asg_coordinate_update(g_i, u, theta_new[i],
                                   check_bracket = check_bracket, box = box_i,
                                   tol = tol / m, scale0 = scale0,
                                   n_grid = n_grid, bisect_e = bisect_e)
      theta_new[i] <- upd$x
      n_retained <- n_retained + upd$retained
    }
    theta_chain[t, ] <- theta_new
  }
  
  sample_idx <- seq(from = burn_in + 1L, to = total_iter, by = thin)
  theta_samples <- theta_chain[sample_idx, , drop = FALSE]
  end_time <- Sys.time()
  time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  ## -----------------------
  ## Plots (per dimension)
  ## -----------------------
  if (make_plots) {
    for (i in seq_len(m)) {
      chain_i <- theta_chain[, i]
      samples_i <- theta_samples[, i]
      a <- bounds_list[[i]]["a"]; b <- bounds_list[[i]]["b"]
      f <- f_list[[i]]
      
      ## 1) Histogram + (initial-conditional) target line + bounds
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
      
      ## 2) Trace
      df_trace <- data.frame(iter = seq_along(chain_i), val = chain_i)
      p_trace <- ggplot(df_trace, aes(iter, val)) +
        geom_line(alpha = 0.8, size = 1.5) +
        geom_vline(xintercept = burn_in, linetype = "dashed", color = "red") +
        labs(title = paste("Dim", i, ": Trace"), x = "Iteration", y = paste0("theta_", i)) +
        theme_minimal(base_size = 12)
      
      ## 3) ACF (lag.max = 40)
      acf_vals <- acf(chain_i, plot = FALSE, lag.max = 50)
      df_acf <- data.frame(lag = as.numeric(acf_vals$lag),
                           acf = as.numeric(acf_vals$acf))
      p_acf <- ggplot(df_acf, aes(lag, acf)) +
        geom_hline(yintercept = 0) +
        geom_segment(aes(xend = lag, yend = 0)) +
        labs(title = paste("Dim", i, ": ACF"), x = "Lag", y = "ACF") +
        theme_minimal(base_size = 12)
      
      ## 4) Running mean
      rmn <- cumsum(chain_i) / seq_along(chain_i)
      df_rm <- data.frame(iter = seq_along(chain_i), mean = rmn)
      p_rm <- ggplot(df_rm, aes(iter, mean)) +
        geom_line(size = 2) +
        geom_vline(xintercept = burn_in, linetype = "dashed", color = "red") +
        labs(title = paste("Dim", i, ": Running mean"), x = "Iteration", y = "Mean") +
        theme_minimal(base_size = 12)
      
      ## Print (you can arrange with gridExtra if you like)
      print(p_density); print(p_trace); print(p_acf); print(p_rm)
    }
  }
  list(samples = theta_samples, chain = theta_chain,
       burn_in = burn_in, thin = thin, time_taken = time_taken,
       retain_rate = n_retained / (total_iter * m))    # [MODIFIED]
}


multivariate_gibbs_sample_ASG <- function(ker, n_samples = 10, burn_in = 2, thin = 1,
                                          theta_init, tol = 0.01, scale0 = 1,
                                          bisect_e = 1e-8,
                                          check_bracket = F,   # [MODIFIED] option
                                          n_grid = 1000L,         # [MODIFIED] grid size for check
                                          make_plots = TRUE) {
  start_time <- Sys.time()
  m <- length(theta_init)
  n_retained <- 0L                                                # [MODIFIED]
  total_iter <- burn_in + n_samples * thin
  theta_chain <- matrix(NA_real_, nrow = total_iter, ncol = m)
  colnames(theta_chain) <- paste0("theta_", seq_len(m))
  theta_new <- as.numeric(theta_init)
  
  for (t in seq_len(total_iter)) {
    
    u <- runif(1, 0, ker(theta_new))   ## fresh u -- every sweep
    
    for (i in seq_len(m)) {
      
      g_i <- (function(i) {           ## fresh g_i -- every sweep, every coordinate
        function(xi) {
          tmp <- theta_new
          tmp[i] <- xi
          ker(tmp)
        }
      })(i)
      
      ## ---- [MODIFIED] check_bracket = TRUE: extreme roots of K = u
      ## (supervisor's check; effective.support() is not called).
      ## check_bracket = FALSE: effective-support box of the CURRENT
      ## conditional, recomputed at every (t, i), with tolerance tol/m. ----
      upd <- asg_coordinate_update(g_i, u, theta_new[i],
                                   check_bracket = check_bracket,
                                   tol = tol / m, scale0 = scale0,
                                   n_grid = n_grid, bisect_e = bisect_e)
      theta_new[i] <- upd$x
      n_retained <- n_retained + upd$retained
    }
    
    theta_chain[t, ] <- theta_new
  }
  
  sample_idx <- seq(from = burn_in + 1L, to = total_iter, by = thin)
  theta_samples <- theta_chain[sample_idx, , drop = FALSE]
  end_time <- Sys.time()
  time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  ## -----------------------------------------
  ## Diagnostics plots (single representative series)
  ## -----------------------------------------
  if (make_plots) {
    ## choose 1D summary to monitor
    diag_series <- theta_chain[, 1]
    series_name <- colnames(theta_chain)[1]
    
    # Trace plot
    df_trace <- data.frame(
      iter = seq_along(diag_series),
      val  = diag_series
    )
    p_trace <- ggplot(df_trace, aes(x = iter, y = val)) +
      geom_line(alpha = 0.8, size = 1.5) +
      geom_vline(xintercept = burn_in, linetype = "dashed", color = "red") +
      labs(
        title = "Trace for theta (first coordinate)",
        x = "Iteration",
        y = series_name
      ) +
      theme_minimal(base_size = 12)   # [MODIFIED] removed hard-coded ylim c(-8, 8)
    
    # ACF plot (lag.max = 50)
    acf_vals <- acf(diag_series, plot = FALSE, lag.max = 50)
    df_acf <- data.frame(
      lag = as.numeric(acf_vals$lag),
      acf = as.numeric(acf_vals$acf)
    )
    p_acf <- ggplot(df_acf, aes(x = lag, y = acf)) +
      geom_hline(yintercept = 0, size = 1) +
      geom_segment(aes(xend = lag, yend = 0)) +
      labs(
        title = "ACF for theta (first coordinate)",
        x = "Lag",
        y = "ACF"
      ) +
      theme_minimal(base_size = 12)
    
    # Running mean plot
    running_mean <- cumsum(diag_series) / seq_along(diag_series)
    df_rm <- data.frame(
      iter = seq_along(diag_series),
      mean = running_mean
    )
    p_rm <- ggplot(df_rm, aes(x = iter, y = mean)) +
      geom_line(size = 2) +
      geom_vline(xintercept = burn_in, linetype = "dashed", color = "red") +
      labs(
        title = "Running mean for theta (first coordinate)",
        x = "Iteration",
        y = "Running mean"
      ) +
      theme_minimal(base_size = 12)
    
    print(p_trace)
    print(p_acf)
    print(p_rm)
    
    ## -----------------------------------------
    ## Log-kernel diagnostics for appendix
    ## log K(theta^(t)) from the full multivariate chain
    ## -----------------------------------------
    tiny <- 1e-300
    logK_series <- apply(theta_chain, 1, function(th) {
      val <- ker(as.numeric(th))
      log(pmax(val, tiny))  # log K(x^{(t)})
    })
    
    # Trace of log K(theta^{(t)})
    df_logk_trace <- data.frame(
      iter = seq_along(logK_series),
      logK = logK_series
    )
    p_logk_trace <- ggplot(df_logk_trace, aes(x = iter, y = logK)) +
      geom_line(alpha = 0.8, size = 2) +
      geom_vline(xintercept = burn_in, linetype = "dashed", color = "red") +
      labs(
        title = "Trace of log-kernel: log K(theta^{(t)})",
        x = "Iteration t",
        y = "log K(theta)"
      ) +
      theme_minimal(base_size = 12)
    
    # ACF of log K(theta^{(t)}) (again, use post-burn-in part)
    acf_logk <- acf(logK_series[(burn_in + 1L):length(logK_series)],
                    plot = FALSE, lag.max = 50)
    df_logk_acf <- data.frame(
      lag = as.numeric(acf_logk$lag),
      acf = as.numeric(acf_logk$acf)
    )
    p_logk_acf <- ggplot(df_logk_acf, aes(x = lag, y = acf)) +
      geom_hline(yintercept = 0, size = 2) +
      geom_segment(aes(xend = lag, yend = 0)) +
      labs(
        title = "ACF of log-kernel: log K(theta^{(t)})",
        x = "Lag",
        y = "ACF"
      ) +
      theme_minimal(base_size = 12)
    
    print(p_logk_trace)
    print(p_logk_acf)
  }
  
  ## -----------------------------------------
  ## Return
  ## -----------------------------------------
  list(
    samples     = theta_samples,
    chain       = theta_chain,
    burn_in     = burn_in,
    thin        = thin,
    time_taken  = time_taken,
    retain_rate = n_retained / (total_iter * m)   # [MODIFIED]
  )
}
