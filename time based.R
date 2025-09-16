# Proposed Method ##################
rm(list = ls())
source("functions_time.R") 
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

result1_30_p <- multivariate_gibbs_sample_time(ker = k1,theta_init = 0, duration_secs = 30)
dim(result1_30_p$samples)

result1_30_p_ess = effectiveSize(result1_30_p$samples)/length(result1_30_p$samples)
result1_30_p_ess

result1_60_p <- multivariate_gibbs_sample_time(ker = k1,theta_init = 0, duration_secs = 60)
result1_60_p_ess = effectiveSize(result1_60_p$samples)/length(result1_60_p$samples)
result1_60_p_ess

result1_90_p <- multivariate_gibbs_sample_time(ker = k1,theta_init = 0, duration_secs = 90)
result1_90_p_ess = effectiveSize(result1_90_p$samples)/length(result1_90_p$samples)
result1_90_p_ess

result1_120_p <- multivariate_gibbs_sample_time(ker = k1,theta_init = 0, duration_secs = 120)
result1_120_p_ess = effectiveSize(result1_120_p$samples)/length(result1_120_p$samples)
result1_120_p_ess


result1_150_p <- multivariate_gibbs_sample_time(ker = k1,theta_init = 0, duration_secs = 150)
result1_150_p_ess = effectiveSize(result1_150_p$samples)/length(result1_150_p$samples)
result1_150_p_ess


log_k1 <- function(theta) {
  dens <- k1(theta[1])
  if (dens <= 0) return(-Inf)
  log(dens)
}

result1_30_mh <- mh_time(initial = 0, log_target = log_k1,run_time = 30)
dim(result1_30_mh$samples)
result1_30_mh_ess = effectiveSize(result1_30_mh$samples)/length(result1_30_mh$samples)
result1_30_mh_ess

result1_60_mh <-  mh_time(initial = 0, log_target = log_k1,run_time = 60)
dim(result1_60_mh$samples)
result1_60_mh_ess = effectiveSize(result1_60_mh$samples)/length(result1_60_mh$samples)
result1_60_mh_ess

result1_90_mh <-  mh_time(initial = 0, log_target = log_k1,run_time = 90)
dim(result1_90_mh$samples)
result1_90_mh_ess = effectiveSize(result1_90_mh$samples)/length(result1_90_mh$samples)
result1_90_mh_ess


result1_120_mh <-  mh_time(initial = 0, log_target = log_k1,run_time = 120)
dim(result1_120_mh$samples)
result1_120_mh_ess = effectiveSize(result1_120_mh$samples)/length(result1_120_mh$samples)
result1_120_mh_ess


result1_150_mh <-  mh_time(initial = 0, log_target = log_k1,run_time = 150)
dim(result1_150_mh$samples)
result1_150_mh_ess = effectiveSize(result1_150_mh$samples)/length(result1_150_mh$samples)
result1_150_mh_ess


################## Plots ##############################

df_proposed_ess_by_samples = c(result1_30_p_ess, result1_60_p_ess, result1_90_p_ess, result1_120_p_ess, result1_150_p_ess)
df_proposed_samples = c(length(result1_30_p$samples), length(result1_60_p$samples), length(result1_90_p$samples), length(result1_120_p$samples), length(result1_150_p$samples))
df_mh = c(result1_30_mh_ess, result1_60_mh_ess, result1_90_mh_ess, result1_120_mh_ess, result1_150_mh_ess)
df_mh_iter = c(length(result1_30_mh$samples),length(result1_60_mh$samples), length(result1_90_mh$samples), length(result1_120_mh$samples), length(result1_150_mh$samples))



proposed_data = data.frame("ESS/N" = df_proposed_ess_by_samples, "Log(N)" = log10(df_proposed_samples), "time" = c(30,60,90,120,150))
proposed_data

MH_data = data.frame("ESS/N" = df_mh, "Log(N)" = log10(df_mh_iter), "time" = c(30,60,90,120,150))
MH_data

# Add an identifier column to each dataset
proposed_data$Method <- "Proposed"
MH_data$Method <- "MH"

# Combine the datasets
all_data <- rbind(proposed_data, MH_data)
all_data
# Plot
ggplot(all_data, aes(x = Log.N., y = ESS.N, color = Method, shape = factor(time))) +
  geom_point(size = 3) +
  labs(
    x = "Log10(N)",
    y = "ESS/N",
    shape = "Time (s)",
    color = "Method"
  ) +
  theme_minimal(base_size = 14)

