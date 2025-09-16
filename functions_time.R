source("effective_support_uni_s.R") 
# Multivariate Gibbs sampling function with density plots, modified to run for a fixed duration
multivariate_gibbs_sample_time <- function(ker, duration_secs , burn_in = 0, thin = 1, tol = 0.1, scale0 = 1, theta_init) {
  start_time <- Sys.time()
  
  m <- length(theta_init)
  
  # Initialize storage for chain
  theta_chain <- list()
  theta_new <- theta_init
  t <- 1
  
  # Lists to store bounds and normalized density for each dimension
  bounds_list <- vector("list", m)
  f_list <- vector("list", m)
  
  # Gibbs sampling loop until duration is reached
  while (as.numeric(difftime(Sys.time(), start_time, units = "secs")) < duration_secs) {
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
    theta_chain[[t]] <- theta_new
    t <- t + 1
  }
  
  # Convert list to matrix for compatibility with downstream code
  total_iter <- length(theta_chain)
  theta_chain <- do.call(rbind, theta_chain)
  
  # Extract samples after burn-in with thinning
  sample_indices <- seq(burn_in + 1, total_iter, by = thin)
  theta_samples <- theta_chain[sample_indices, , drop = FALSE]
  
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
    # Add the time taken on the left side
    p_density <- p_density + annotate("text", x = -Inf, y = Inf, label = sprintf("Time Taken: %.2f sec", time_taken),
                                      hjust = 0.005, vjust = 1.1, size = 5, color = "black")
    
    # 2. Trace Plot
    df_chain <- data.frame(Iteration = 1:total_iter, X = chain_i)
    p1 <- ggplot(df_chain, aes(x = Iteration, y = X)) +
      geom_line(color = "steelblue", alpha = 0.7) +
      labs(title = paste("Trace Plot of MCMC Chain for Dimension", i),
           x = "Iteration", y = paste("theta_", i, sep = "")) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            panel.grid.minor = element_blank()) 
    
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
      labs(title = paste("Running Mean for Dimension", i),
           x = "Iteration", y = "Running Mean") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            panel.grid.minor = element_blank()) 
    
    # Print the plots
    print(p_density)
    print(p1)
    print(p2)
    print(p3)
  }
  
  return(list(samples = theta_samples, chain = theta_chain, burn_in = burn_in, thin = thin, time_taken = time_taken))
}


################### MH algorithm ######################

library(mcmc)
library(coda)


mh_time <- function(initial, log_target, scale = 1, run_time ) {
  start_time <- Sys.time()
  batch_size = 1000
  current <- initial
  samples <- list()
  i <- 1
  
  while (as.numeric(difftime(Sys.time(), start_time, units = "secs")) < run_time) {
    # Do a batch at once
    out <- metrop(obj = log_target, initial = current,
                  nbatch = batch_size, blen = 1, scale = scale)
    
    samples[[i]] <- out$batch
    current <- out$final
    i <- i + 1
  }
  
  end_time <- Sys.time()
  total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  list(
    samples = do.call(rbind, samples),   # all samples
    elapsed = total_time                 # actual runtime in seconds
  )
}

