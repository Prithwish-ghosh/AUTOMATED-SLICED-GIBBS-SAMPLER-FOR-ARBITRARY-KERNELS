library(mcmcse)
# Proposed Method ##################
#source("time fixed function.R") 
############### Results #################################

# Run the multivariate Gibbs sampler
## Univariate
k1 <- function(theta) {
  0.3 * dbeta((theta[1] + 5) / 2, 0.5, 1) / 2 + 
    0.4 * dbeta((theta[1] + 1) / 2, 2, 2) / 2 + 
    0.3 * dbeta((theta[1] - 5) / 2, 1, 0.5) / 2
}



set.seed(2025)

library(coda)

result1_10_p <- multivariate_gibbs_sample_ASG_time(ker = k1,theta_init = 0, duration_secs = 10)
dim(result1_10_p$samples)

result1_10_p_ess = effectiveSize(result1_10_p$samples)/length(result1_10_p$samples)
result1_10_p_ess

result1_20_p <- multivariate_gibbs_sample_ASG_time(ker = k1,theta_init = 0, duration_secs = 20)
result1_20_p_ess = effectiveSize(result1_20_p$samples)/length(result1_20_p$samples)
result1_20_p_ess

result1_30_p <- multivariate_gibbs_sample_ASG_time(ker = k1,theta_init = 0, duration_secs = 30)
result1_30_p_ess = effectiveSize(result1_30_p$samples)/length(result1_30_p$samples)
result1_30_p_ess


log_k1 <- function(theta) {
  dens <- k1(theta[1])
  if (dens <= 0) return(-Inf)
  log(dens)
}

result1_10_mh <- mh_time(initial = 0, log_target = log_k1,run_time = 10)
dim(result1_10_mh$samples)
result1_10_mh_ess = effectiveSize(result1_10_mh$samples)/length(result1_10_mh$samples)
result1_10_mh_ess

result1_20_mh <-  mh_time(initial = 0, log_target = log_k1,run_time = 20)
dim(result1_20_mh$samples)
result1_20_mh_ess = effectiveSize(result1_20_mh$samples)/length(result1_20_mh$samples)
result1_20_mh_ess

result1_30_mh <-  mh_time(initial = 0, log_target = log_k1,run_time = 30)
dim(result1_30_mh$samples)
result1_30_mh_ess = effectiveSize(result1_30_mh$samples)/length(result1_30_mh$samples)
result1_30_mh_ess


################## Plots ready for univariate ##############################

df_proposed_ess_by_samples = c(result1_10_p_ess, result1_20_p_ess, result1_30_p_ess)
df_proposed_samples = c(length(result1_10_p$samples), length(result1_20_p$samples), length(result1_30_p$samples))
df_mh = c(result1_10_mh_ess, result1_20_mh_ess, result1_30_mh_ess)
df_mh_iter = c(length(result1_10_mh$samples),length(result1_20_mh$samples), length(result1_30_mh$samples))



proposed_data = data.frame("Multi ESS/N" = df_proposed_ess_by_samples, "Log(N)" = log10(df_proposed_samples), "time" = c(10,20,30))
proposed_data

MH_data = data.frame("Multi ESS/N" = df_mh, "Log(N)" = log10(df_mh_iter), "time" = c(10,20,30))
MH_data

# Add an identifier column to each dataset
proposed_data$Method <- "ASG"
MH_data$Method <- "MH"

# Combine the datasets
all_data <- rbind(proposed_data, MH_data)
all_data
all_data$kernel = "univariate"





#### Banana Shaped kernel ####

banana_kernel <- function(theta) {
  x1 <- theta[1]
  x2 <- theta[2]
  log_pi <- -x1^2 / 10 - x2^2 / 10 - 2 * (x2 - x1^2)^2
  exp(log_pi)
}


set.seed(2026)

theta_init = c(0,0)


##  10 sec

result1_10_p_b <- multivariate_gibbs_sample_ASG_time(ker = banana_kernel,theta_init , duration_secs = 10)
dim(result1_10_p_b$samples)
ess_avg_b_10 = multiESS(result1_10_p_b$samples)

result1_10_p_ess_b = ess_avg_b_10/length(result1_10_p_b$samples[,1])
result1_10_p_ess_b

#### 20 sec

result1_20_p_b <- multivariate_gibbs_sample_ASG_time(ker = banana_kernel,theta_init , duration_secs = 20)
ess_avg_b_20 = multiESS(result1_10_p_b$samples)

result1_20_p_ess_b = ess_avg_b_20/length(result1_20_p_b$samples[,1])
result1_20_p_ess_b

### 30 sec

result1_30_p_b <- multivariate_gibbs_sample_ASG_time(ker = banana_kernel,theta_init , duration_secs = 30)
ess_avg_b_30 = multiESS(result1_30_p_b$samples)

result1_30_p_ess_b = ess_avg_b_30/length(result1_30_p_b$samples[,1])
result1_30_p_ess_b


# Log-density for MH sampling
log_banana <- function(theta) {
  dens <- banana_kernel(theta)
  if (dens <= 0) return(-Inf)
  log(dens)
}

### 10 Sec

result1_10_mh_b <- mh_time(initial = theta_init, log_target = log_banana,run_time = 10)
dim(result1_10_mh_b$samples)
ess_avg_b_mh_10 = multiESS(result1_10_mh_b$samples)

result1_10_p_ess_b_mh = ess_avg_b_mh_10/length(result1_10_mh_b$samples[,1])
result1_10_p_ess_b_mh

### 20 Sec

result1_20_mh_b <-  mh_time(initial = theta_init, log_target = log_banana ,run_time = 20)
dim(result1_20_mh_b$samples)
ess_avg_b_mh_20 = multiESS(result1_20_mh_b$samples)

result1_20_p_ess_b_mh = ess_avg_b_mh_20/length(result1_20_mh_b$samples[,1])
result1_20_p_ess_b_mh

### 30 sec

result1_30_mh_b <-  mh_time(initial = theta_init, log_target = log_banana ,run_time = 30)
dim(result1_30_mh_b$samples)
ess_avg_b_mh_30 = multiESS(result1_30_mh_b$samples)

result1_30_p_ess_b_mh = ess_avg_b_mh_30/length(result1_30_mh_b$samples[,1])
result1_30_p_ess_b_mh


################## Plots ready for banana shaped kernel ##############################

df_proposed_ess_by_samples_b = c(result1_10_p_ess_b, result1_20_p_ess_b, result1_30_p_ess_b)
df_proposed_samples_b = c(length(result1_10_p_b$samples[,1]), length(result1_20_p_b$samples[,1]), length(result1_30_p_b$samples[,1]))
df_mh_b = c(result1_10_p_ess_b_mh, result1_20_p_ess_b_mh, result1_30_p_ess_b_mh)
df_mh_iter_b = c(length(result1_10_mh_b$samples),length(result1_20_mh_b$samples), length(result1_30_mh_b$samples))



proposed_data_b = data.frame("Multi ESS/N" = df_proposed_ess_by_samples_b, "Log(N)" = log10(df_proposed_samples_b), "time" = c(10,20,30))
proposed_data_b

MH_data_b = data.frame("Multi ESS/N" = df_mh_b, "Log(N)" = log10(df_mh_iter_b), "time" = c(10,20,30))
MH_data_b

# Add an identifier column to each dataset
proposed_data_b$Method <- "ASG"
MH_data_b$Method <- "MH"

# Combine the datasets
all_data_b <- rbind(proposed_data_b, MH_data_b)
all_data_b
all_data_b$kernel = "Rosenbrock Kernel"


all_dataset = rbind(all_data, all_data_b)
all_dataset










#### Ackley function kernel ####

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

# Kernel (unnormalized density): π(θ) ∝ exp( - f(θ) / temp )
# For vector theta returns scalar; for matrix returns vector (rowwise).
# Prefer using ackley_log_kernel in samplers to avoid overflow/underflow.

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

ackley_kernel <- function(theta, temp = 1, box = c(-32.768, 32.768)) {
  lk <- ackley_log_kernel(theta, temp = temp, box = box)
  exp(lk)
}


set.seed(2026)

theta_init = c(0,0)


##  10 sec

result1_10_p_a <- multivariate_gibbs_sample_ASG_time(ker = ackley_kernel,theta_init , duration_secs = 10)
dim(result1_10_p_a$samples)
ess_avg_a_10 = multiESS(result1_10_p_a$samples)

result1_10_p_ess_a = ess_avg_a_10/length(result1_10_p_a$samples[,1])
result1_10_p_ess_a

#### 20 sec

result1_20_p_a <- multivariate_gibbs_sample_ASG_time(ker = ackley_kernel,theta_init , duration_secs = 20)
ess_avg_a_20 = multiESS(result1_20_p_a$samples)

result1_20_p_ess_a = ess_avg_a_20/length(result1_20_p_a$samples[,1])
result1_20_p_ess_a

### 30 sec

result1_30_p_a <- multivariate_gibbs_sample_ASG_time(ker = ackley_kernel,theta_init , duration_secs = 30)
ess_avg_a_30 = multiESS(result1_30_p_a$samples)

result1_30_p_ess_a = ess_avg_a_30/length(result1_30_p_a$samples[,1])
result1_30_p_ess_a


### 10 Sec

result1_10_mh_a <- mh_time(initial = theta_init, log_target = ackley_log_kernel,run_time = 10)
dim(result1_10_mh_a$samples)
ess_avg_a_mh_10 = multiESS(result1_10_mh_a$samples)

result1_10_p_ess_a_mh = ess_avg_a_mh_10/length(result1_10_mh_a$samples[,1])
result1_10_p_ess_a_mh

### 20 Sec

result1_20_mh_a <-  mh_time(initial = theta_init, log_target = ackley_log_kernel ,run_time = 20)
dim(result1_20_mh_a$samples)
ess_avg_a_mh_20 = multiESS(result1_20_mh_a$samples)

result1_20_p_ess_a_mh = ess_avg_a_mh_20/length(result1_20_mh_a$samples[,1])
result1_20_p_ess_a_mh

### 30 sec

result1_30_mh_a <-  mh_time(initial = theta_init, log_target = ackley_log_kernel ,run_time = 30)
dim(result1_30_mh_a$samples)
ess_avg_a_mh_30 = multiESS(result1_30_mh_a$samples)

result1_30_p_ess_a_mh = ess_avg_a_mh_30/length(result1_30_mh_a$samples[,1])
result1_30_p_ess_a_mh


################## Plots ready for banana shaped kernel ##############################

df_proposed_ess_by_samples_a = c(result1_10_p_ess_a, result1_20_p_ess_a, result1_30_p_ess_a)
df_proposed_samples_a = c(length(result1_10_p_a$samples[,1]), length(result1_20_p_a$samples[,1]), length(result1_30_p_a$samples[,1]))
df_mh_a = c(result1_10_p_ess_a_mh, result1_20_p_ess_a_mh, result1_30_p_ess_a_mh)
df_mh_iter_a = c(length(result1_10_mh_a$samples),length(result1_20_mh_a$samples), length(result1_30_mh_a$samples))



proposed_data_a = data.frame("Multi ESS/N" = df_proposed_ess_by_samples_a, "Log(N)" = log10(df_proposed_samples_a), "time" = c(10,20,30))
proposed_data_a

MH_data_a = data.frame("Multi ESS/N" = df_mh_a, "Log(N)" = log10(df_mh_iter_a), "time" = c(10,20,30))
MH_data_a

# Add an identifier column to each dataset
proposed_data_a$Method <- "ASG"
MH_data_a$Method <- "MH"

# Combine the datasets
all_data_a <- rbind(proposed_data_a, MH_data_a)
all_data_a
all_data_a$kernel = "Ackley Kernel"


all_dataset_final = rbind(all_dataset, all_data_a)
all_dataset_final

#### Plot


ggplot(all_dataset_final,
       aes(x = Log.N.,
           y = Multi.ESS.N,
           color = Method,
           shape = factor(time),
           linetype = kernel)) +
  geom_point(size = 5) +
  geom_line(aes(group = interaction(Method, kernel)), linewidth = 1) +
  labs(
    x = "Log10(N)",
    y = "Multi ESS",
    shape = "Time (s)",
    color = "Method",
    linetype = "Kernel"
  ) +
  theme_minimal(base_size = 14)



## Recover N, ESS and ESS per second from the time-fixed results
all_dataset_final$N           <- 10^all_dataset_final$Log.N.
all_dataset_final$ESS         <- all_dataset_final$Multi.ESS.N * all_dataset_final$N
all_dataset_final$ESS_per_sec <- all_dataset_final$ESS / all_dataset_final$time

## ESS per second (time fixed)
ggplot(all_dataset_final,
       aes(x = factor(time), y = ESS_per_sec, fill = Method)) +
  geom_col(position = "dodge") +
  facet_wrap(~ kernel, scales = "free_y") +
  labs(x = "Run time (seconds)", y = "ESS per second") +
  theme_minimal(base_size = 14)

## ESS/N (time fixed), same bar layout
ggplot(all_dataset_final,
       aes(x = factor(time), y = Multi.ESS.N, fill = Method)) +
  geom_col(position = "dodge") +
  facet_wrap(~ kernel, scales = "free_y") +
  labs(x = "Run time (seconds)", y = "Multi ESS / N") +
  theme_minimal(base_size = 14)
