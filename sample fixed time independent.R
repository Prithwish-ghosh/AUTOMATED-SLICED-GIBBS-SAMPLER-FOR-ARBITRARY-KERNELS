rm(list = ls())
library(coda)
library(ggplot2)
setwd("~/Downloads/ASG Sampler/code")
source("effective_support_uni_s.R")
source("asg_sampler.R")
source("time based analysis examples and plots.R")

##### 1000 Samples

n_samples = 1000
burn_in = 250
nbatch =  1000

################### ASG case #####################

# Univariate Case
result1_k1_1000 <- multivariate_gibbs_sample_ASG(ker = k1,theta_init = 0,n_samples,burn_in)
k1_1000 = effectiveSize(result1_k1_1000$samples)

# Rosenbrock Kernel
result1_banana_1000 <- multivariate_gibbs_sample_ASG(ker = banana_kernel,theta_init = c(0,0),n_samples,burn_in)
banana_1000 = mean(effectiveSize(result1_banana_1000$samples))
banana_1000

# Ackley Kernel
result1_ackley_1000 <- multivariate_gibbs_sample_ASG(ker = ackley_kernel,theta_init = c(0,0),n_samples,burn_in)
ackley_1000 = mean(effectiveSize(result1_ackley_1000$samples))
ackley_1000

################### MH Case #########################

# Univariate kernle
time_k1_1000_mh <- system.time({
  out_k1_1000_mh <- metrop(obj = log_k1, initial = c(0), nbatch = nbatch + burn_in, blen = 1, scale = 1)
})
k1_1000_mh = effectiveSize(out_k1_1000_mh$batch)

# Banana Kernel 

time_banana_1000_mh <- system.time({
  out_banana_1000_mh <- metrop(obj = log_banana, initial = c(0,0), nbatch = nbatch + burn_in, blen = 1, scale = 1)
})
banana_1000_mh = mean(effectiveSize(out_banana_1000_mh$batch))
banana_1000_mh

# Ackley Kernel

time_ackley_1000_mh <- system.time({
  out_ackley_1000_mh <- metrop(obj = ackley_log_kernel, initial = c(0,0), nbatch = nbatch + burn_in, blen = 1, scale = 1)
})
ackley_1000_mh = mean(effectiveSize(out_ackley_1000_mh$batch))
ackley_1000_mh



############## Plot data creation ###############

# ASG PART
ess_samples_asg = c(k1_1000, banana_1000, ackley_1000)/1000
ess_samples_asg
time_asg = c(result1_k1_1000$time_taken,result1_banana_1000$time_taken, result1_ackley_1000$time_taken)

dataset_1000_asg = data.frame("ESS/N" = ess_samples_asg, "time" = time_asg, "samples" = 1000, "method" = "ASG")
dataset_1000_asg

# MH PART
ess_samples_mh = c(k1_1000_mh, banana_1000_mh, ackley_1000_mh)/1000
ess_samples_mh
time_mh = c(time_k1_1000_mh[3],time_banana_1000_mh[3], time_ackley_1000_mh[3])

dataset_1000_mh = data.frame("ESS/N" = ess_samples_mh, "time" = time_mh, "samples" = 1000, "method" = "MH")
dataset_1000_mh

final_data_1000 = rbind(dataset_1000_asg, dataset_1000_mh)
final_data_1000$kernel = c("univariate", "Rosenbrock", "ackley", "univariate", "Rosenbrock", "ackley")
final_data_1000
###### 5000 samples ##############

n_samples = 5000
burn_in = 500
nbatch =  5000

################### ASG case #####################

# Univariate Case
result1_k1_5000 <- multivariate_gibbs_sample_ASG(ker = k1,theta_init = 0,n_samples,burn_in)
k1_5000 = effectiveSize(result1_k1_5000$samples)

# Rosenbrock Kernel
result1_banana_5000 <- multivariate_gibbs_sample_ASG(ker = banana_kernel,theta_init = c(0,0),n_samples,burn_in)
banana_5000 = mean(effectiveSize(result1_banana_5000$samples))
banana_5000

# Ackley Kernel
result1_ackley_5000 <- multivariate_gibbs_sample_ASG(ker = ackley_kernel,theta_init = c(0,0),n_samples,burn_in)
ackley_5000 = mean(effectiveSize(result1_ackley_5000$samples))
ackley_5000

################### MH Case #########################

# Univariate kernle
time_k1_5000_mh <- system.time({
  out_k1_5000_mh <- metrop(obj = log_k1, initial = c(0), nbatch = nbatch + burn_in, blen = 1, scale = 1)
})
k1_5000_mh = effectiveSize(out_k1_5000_mh$batch)

# Banana Kernel 

time_banana_5000_mh <- system.time({
  out_banana_5000_mh <- metrop(obj = log_banana, initial = c(0,0), nbatch = nbatch + burn_in, blen = 1, scale = 1)
})
banana_5000_mh = mean(effectiveSize(out_banana_5000_mh$batch))
banana_5000_mh

# Ackley Kernel

time_ackley_5000_mh <- system.time({
  out_ackley_5000_mh <- metrop(obj = ackley_log_kernel, initial = c(0,0), nbatch = nbatch + burn_in, blen = 1, scale = 1)
})
ackley_5000_mh = mean(effectiveSize(out_ackley_5000_mh$batch))
ackley_5000_mh


########### Plot data generation ########

# ASG PART
ess_samples_asg = c(k1_5000, banana_5000, ackley_5000)/5000
ess_samples_asg
time_asg = c(result1_k1_5000$time_taken,result1_banana_5000$time_taken, result1_ackley_5000$time_taken)

dataset_5000_asg = data.frame("ESS/N" = ess_samples_asg, "time" = time_asg, "samples" = 5000, "method" = "ASG")
dataset_5000_asg

# MH PART
ess_samples_mh = c(k1_5000_mh, banana_5000_mh, ackley_5000_mh)/5000
ess_samples_mh
time_mh = c(time_k1_5000_mh[3],time_banana_5000_mh[3], time_ackley_5000_mh[3])

dataset_5000_mh = data.frame("ESS/N" = ess_samples_mh, "time" = time_mh, "samples" = 5000, "method" = "MH")
dataset_5000_mh

final_data_5000 = rbind(dataset_5000_asg, dataset_5000_mh)
final_data_5000
final_data_5000$kernel = c("univariate", "Rosenbrock", "ackley", "univariate", "Rosenbrock", "ackley")


### 10000 sampling #########


n_samples = 10000
burn_in = 1000
nbatch =  10000

################### ASG case #####################

# Univariate Case
result1_k1_10000 <- multivariate_gibbs_sample_ASG(ker = k1,theta_init = 0,n_samples,burn_in)
k1_10000 = effectiveSize(result1_k1_10000$samples)

# Rosenbrock Kernel
result1_banana_10000 <- multivariate_gibbs_sample_ASG(ker = banana_kernel,theta_init = c(0,0),n_samples,burn_in)
banana_10000 = mean(effectiveSize(result1_banana_10000$samples))
banana_10000

# Ackley Kernel
result1_ackley_10000 <- multivariate_gibbs_sample_ASG(ker = ackley_kernel,theta_init = c(0,0),n_samples,burn_in)
ackley_10000 = mean(effectiveSize(result1_ackley_10000$samples))
ackley_10000

################### MH Case #########################

# Univariate kernle
time_k1_10000_mh <- system.time({
  out_k1_10000_mh <- metrop(obj = log_k1, initial = c(0), nbatch = nbatch + burn_in, blen = 1, scale = 1)
})
k1_10000_mh = effectiveSize(out_k1_10000_mh$batch)

# Banana Kernel 

time_banana_10000_mh <- system.time({
  out_banana_10000_mh <- metrop(obj = log_banana, initial = c(0,0), nbatch = nbatch + burn_in, blen = 1, scale = 1)
})
banana_10000_mh = mean(effectiveSize(out_banana_10000_mh$batch))
banana_10000_mh

# Ackley Kernel

time_ackley_10000_mh <- system.time({
  out_ackley_10000_mh <- metrop(obj = ackley_log_kernel, initial = c(0,0), nbatch = nbatch + burn_in, blen = 1, scale = 1)
})
ackley_10000_mh = mean(effectiveSize(out_ackley_10000_mh$batch))
ackley_10000_mh

######## Plot data generation ######

# ASG PART
ess_samples_asg = c(k1_10000, banana_10000, ackley_10000)/10000
ess_samples_asg
time_asg = c(result1_k1_10000$time_taken,result1_banana_10000$time_taken, result1_ackley_10000$time_taken)

dataset_10000_asg = data.frame("ESS/N" = ess_samples_asg, "time" = time_asg, "samples" = 10000, "method" = "ASG")
dataset_10000_asg

# MH PART
ess_samples_mh = c(k1_10000_mh, banana_10000_mh, ackley_10000_mh)/10000
ess_samples_mh
time_mh = c(time_k1_10000_mh[3],time_banana_10000_mh[3], time_ackley_10000_mh[3])

dataset_10000_mh = data.frame("ESS/N" = ess_samples_mh, "time" = time_mh, "samples" = 10000, "method" = "MH")
dataset_10000_mh

final_data_10000 = rbind(dataset_10000_asg, dataset_10000_mh)
final_data_10000
final_data_10000$kernel = c("univariate", "Rosenbrock", "ackley", "univariate", "Rosenbrock", "ackley")


#### Final plot dataset #######
plot_final_dataset = rbind(final_data_1000, final_data_5000, final_data_10000)
plot_final_dataset$time = (plot_final_dataset$time)*10

ggplot(plot_final_dataset,
       aes(x = time,
           y = ESS.N,
           color = method,
           shape = factor(samples),
           linetype = kernel)) +
  geom_point(size = 5) +
  geom_line(aes(group = interaction(method, kernel)), linewidth = 1) +
  labs(
    x = "Time X 10",
    y = "ESS/N",
    shape = "Time (s)",
    color = "Method",
    linetype = "Kernel"
  ) +
  theme_minimal(base_size = 14)
