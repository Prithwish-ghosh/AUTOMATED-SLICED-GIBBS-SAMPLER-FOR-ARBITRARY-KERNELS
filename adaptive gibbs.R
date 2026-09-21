library(rhmc)
library(qslice)
library(mcmc)
library(spBayes)
library(coda)
library(mcmcse)
## ----- Kernel: 3-component Beta mixture on R -----
k1 <- function(theta) {
  0.3 * dbeta((theta + 5) / 2, 0.5, 1) / 2 +
    0.4 * dbeta((theta + 1) / 2, 2, 2) / 2 +
    0.3 * dbeta((theta - 5) / 2, 1, 0.5) / 2
}


## ----- Minus log-density for HMC -----
## hmc() wants: f(x) = -log pi(x)  (up to a constant)
f_k1 <- function(x) {
  dens <- k1(x)
  ## add tiny epsilon so log(0) doesn't explode to -Inf
  -log(dens + 1e-300)
}

set.seed(123)

## HMC parameters
numit <- 1250   # total iterations = burn-in (250) + kept (1000)
L     <- 10     # leapfrog steps
eps   <- 0.05   # step size
mass  <- 1      # mass (scalar for 1D)

## initial point (somewhere in a region with non-tiny density)
init  <- 0

time_hmc <- system.time(
  hmc_out <- hmc(
    f     = f_k1,
    init  = init,
    numit = numit,
    L     = L,
    eps   = eps,
    mass  = mass
    # num_grad = 3     # only if your implementation supports it
  )
)

time_hmc
samples_univariate_case = (hmc_out$chain[1,])

acf(samples_univariate_case, lag.max = 50)
effectiveSize(samples_univariate_case)

plot(hmc_out$chain[1,], type = "l", main = "Trace plot for k1", xlab = "Iteration", ylab = "Theta")

mean(samples_univariate_case)

## Ackley
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

# Negative log-kernel for the Ackley target
L_ackley <- function(theta, a = 20, b = 0.2, c = 2*pi, temp = 1) {
  f <- ackley_f(theta, a = a, b = b, c = c)  # your function above
  f / temp
}

time_hmc_ackley <- system.time(
hmc_out <- hmc(
  f     = L_ackley,
  init  = c(0,0),
  numit = numit,
  L     = 10,
  eps   = eps,
  mass  = mass
  # if your implementation has it:
  # num_grad = 3
)
)

time_hmc_ackley

mode_x <- names(which.max(table(hmc_out$chain[1,])))
mode_x
mode_y <- names(which.max(table(hmc_out$chain[2,])))
mode_y
library(coda)

library(coda)

# Your HMC output: 1 x 1000
hmc_chain <- as.matrix(hmc_out$chain)

# Convert to iterations x parameters
hmc_chain <- t(hmc_chain)

# Multivariate ESS
multiESS(hmc_chain)

effectiveSize(hmc_out$chain[1,])
effectiveSize(hmc_out$chain[2,])



acf(hmc_out$chain[1,], lag.max = 50)
acf(hmc_out$chain[2,], lag.max = 50)

# Define the banana-shaped kernel (unnormalized density)
banana_kernel <- function(theta) {
  x1 <- theta[1]
  x2 <- theta[2]
  log_pi <- -x1^2 / 10 - x2^2 / 10 - 2 * (x2 - x1^2)^2
  exp(log_pi)
}


L_banana <- function(theta) {
  x1 <- theta[1]
  x2 <- theta[2]
  # negative log-kernel (Banana / Rosenbrock type)
  (x1^2) / 10 + (x2^2) / 10 + 2 * (x2 - x1^2)^2
}

time_hmc_banana <- system.time(
hmc_out <- hmc(
  f     = L_banana,
  init  = c(0,0),
  numit = numit,
  L     = 10,
  eps   = eps,
  mass  = mass
  # if your implementation has it:
  # num_grad = 3
)
)

mode_x <- names(which.max(table(hmc_out$chain[1,])))
mode_x
mode_y <- names(which.max(table(hmc_out$chain[2,])))
mode_y

time_hmc_banana


library(coda)

# Your HMC output: 1 x 1000
hmc_chain <- as.matrix(hmc_out$chain)

# Convert to iterations x parameters
hmc_chain <- t(hmc_chain)

# Multivariate ESS
multiESS(hmc_chain)

effectiveSize(hmc_out$chain[1,])
effectiveSize(hmc_out$chain[2,])
acf(hmc_out$chain[1,], lag.max = 50)


#### Adabtiive Gibbs ###

library(spBayes)
#?spBayes::adaptMetropGibbs()

#library(spBayes)

## ----- Kernel: 3-component Beta mixture -----
k1 <- function(theta) {
  x <- theta[1]  # adaptMetropGibbs passes a vector; we use the first coord
  0.3 * dbeta((x + 5) / 2, 0.5, 1) / 2 +
    0.4 * dbeta((x + 1) / 2, 2, 2) / 2 +
    0.3 * dbeta((x - 5) / 2, 1, 0.5) / 2
}

## log target density for adaptMetropGibbs
ltd_k1 <- function(theta) {
  val <- k1(theta)
  log(val + 1e-300)   # avoid log(0)
}



set.seed(123)

n_total       <- 1250      # 250 burn-in + 1000 kept
batch         <- 25
batch.length  <- 50        # 25 * 50 = 1250

time_adap_gibbs <- system.time(
  amg_out <- adaptMetropGibbs(
  ltd         = ltd_k1,
  starting    = c(0),      # initial theta
  tuning      = 0.5,       # initial proposal sd (will be adapted)
  accept.rate = 0.5,      # target acceptance
  batch       = batch,
  batch.length = batch.length,
  report      = 5,
  verbose     = TRUE
))

time_adap_gibbs


acf(amg_out$p.theta.samples, lag.max = 50)
effectiveSize(amg_out$p.theta.samples)


mean(amg_out$p.theta.samples)


## Ackley function
## Log-kernel of Ackley (unnormalized log-density)
ackley_log_kernel <- function(theta, a = 20, b = 0.2, c = 2*pi) {
  if (is.null(dim(theta))) {
    theta <- matrix(theta, nrow = 1)
  }
  d <- ncol(theta)
  
  ss <- rowSums(theta^2)
  cos_mean <- rowMeans(cos(c * theta))
  
  # Ackley function value
  f <- -a * exp(-b * sqrt(ss / d)) -
    exp(cos_mean) +
    a + exp(1)
  
  # log kernel (up to additive constant)
  log_pi <- -f
  as.numeric(log_pi)
}

n_total       <- 1250      # 250 burn-in + 1000 kept
batch         <- 25
batch.length  <- 50        # 25 * 50 = 1250

time_adap_gibbs_ackley <- system.time(
  amg_out_ackley <- adaptMetropGibbs(
    ltd         = ackley_log_kernel,
    starting    = c(0,0),      # initial theta
    tuning      = 0.5,       # initial proposal sd (will be adapted)
    accept.rate = 0.5,      # target acceptance
    batch       = batch,
    batch.length = batch.length,
    report      = 5,
    verbose     = TRUE
  ))

acf(amg_out_ackley$p.theta.samples, lag.max = 50)
multiESS(amg_out$p.theta.samples)

time_adap_gibbs_ackley
mean(amg_out$p.theta.samples)


### Banana Shaped Kernel

############################################################
## Banana-shaped Rosenbrock kernel
## Adaptive Metropolis-within-Gibbs (adaptMetropGibbs)
## Single self-contained R script
############################################################

## If needed (first time):
## install.packages("adaptMCMC")
## install.packages("coda")

library(adaptMCMC)
library(coda)

## ---------------------------------------------------------
## 1. Banana (Rosenbrock-type) kernel + log-kernel
## ---------------------------------------------------------
banana_kernel <- function(theta) {
  x1 <- theta[1]
  x2 <- theta[2]
  log_pi <- -x1^2 / 10 - x2^2 / 10 - 2 * (x2 - x1^2)^2
  exp(log_pi)
}

banana_log_kernel <- function(theta) {
  x1 <- theta[1]
  x2 <- theta[2]
  -x1^2 / 10 - x2^2 / 10 - 2 * (x2 - x1^2)^2
}

## log target density for adaptMetropGibbs
## (must return a scalar log density)
ltd_banana <- function(theta) {
  banana_log_kernel(theta)
}

## ---------------------------------------------------------
## 2. Run Adaptive Metropolis-within-Gibbs
## ---------------------------------------------------------
set.seed(123)

## total iterations = batch * batch.length
n_keep        <- 1000
n_burn        <- 200
n_total       <- n_burn + n_keep

batch         <- 50
batch.length  <- n_total / batch
if (batch.length != as.integer(batch.length)) {
  stop("Choose batch and batch.length so that batch * batch.length = n_total exactly.")
}

starting <- c(0, 0)

## initial proposal SDs per coordinate (will be adapted)
tuning_init <- c(0.8, 1.2)

time_adap_gibbs <- system.time(
  amg_out <- adaptMetropGibbs(
    ltd          = ltd_banana,
    starting     = starting,
    tuning       = tuning_init,
    accept.rate  = 0.44,     # typical target for 1D RW updates; 0.4-0.5 is common
    batch        = batch,
    batch.length = batch.length,
    report       = 5,
    verbose      = TRUE
  )
)

print(time_adap_gibbs)

## ---------------------------------------------------------
## 3. Extract samples
## ---------------------------------------------------------
## adaptMetropGibbs stores samples in $p.theta.samples
samps_all <- amg_out$p.theta.samples

if (is.vector(samps_all)) {
  ## just in case it returns a vector (shouldn't for 2D)
  samps_all <- matrix(samps_all, ncol = 2, byrow = TRUE)
}

colnames(samps_all) <- c("x1", "x2")

## burn-in removal
samps <- samps_all[(n_burn + 1):n_total, , drop = FALSE]

cat("\nAcceptance rates (per coordinate, if available):\n")
if (!is.null(amg_out$acceptance.rate)) {
  print(amg_out$acceptance.rate)
} else {
  cat("Not provided by this version/object.\n")
}

cat("\nPosterior means (after burn-in):\n")
print(colMeans(samps))

cat("\nEffective Sample Sizes (after burn-in):\n")
print(multiESS(as.mcmc(samps)))

## ---------------------------------------------------------
## 4. Diagnostics
## ---------------------------------------------------------
par(mfrow = c(2, 2))

## Trace
plot(samps[,1], type="l", main=expression("Trace: " * x[1]), ylab="")
plot(samps[,2], type="l", main=expression("Trace: " * x[2]), ylab="")

## ACF
acf(samps[,1], lag.max = 50, main=expression("ACF: " * x[1]))
acf(samps[,2], lag.max = 50, main=expression("ACF: " * x[2]))

## ---------------------------------------------------------
## 5. Banana scatter
## ---------------------------------------------------------
par(mfrow = c(1, 1))
plot(
  samps[,1], samps[,2],
  pch = 16, cex = 0.4,
  xlab = expression(x[1]),
  ylab = expression(x[2]),
  main = "Banana-shaped Rosenbrock (Adaptive Metropolis-within-Gibbs samples)"
)



### qslice

library(qslice)

## ----- Kernel: 3-component Beta mixture -----
k1 <- function(theta) {
  x <- theta[1]
  0.3 * dbeta((x + 5) / 2, 0.5, 1) / 2 +
    0.4 * dbeta((x + 1) / 2, 2, 2) / 2 +
    0.3 * dbeta((x - 5) / 2, 1, 0.5) / 2
}

## log-target for slice_quantile (scalar x)
lf_k1 <- function(x) {
  val <- k1(x)
  log(val + 1e-300)   # avoid log(0)
}

## ----- Pseudo-target: Uniform(-8, 8) -----
## (covers all three modes: [-5,-3], [-1,1], [5,7])
pseudo <- list(
  ld = function(x) dunif(x, min = -8, max = 8, log = TRUE),
  q  = function(u) qunif(u, min = -8, max = 8)
)

## ----- Quantile slice sampling -----
n_draws <- 1000
draws <- numeric(n_draws)
draws[1] <- 0      # start in the middle mode
nEvaluations <- 0L

time_qslice_gibbs <- system.time(
for (i in 2:n_draws) {
  out <- slice_quantile(
    x          = draws[i - 1],
    log_target = lf_k1,
    pseudo     = pseudo
  )
  draws[i] <- out$x
  nEvaluations <- nEvaluations + out$nEvaluations
}
)

## draws now contains 1000 samples from the k1 mixture via quantile slice
## nEvaluations is total number of target evaluations us

time_qslice_gibbs

draws
mean(draws)
effectiveSize(draws)
acf(draws)
plot(draws, type = "l")


## Ackley Function

## =========================
## Ackley log-kernel (dD) + Gibbs qslice (coordinate updates)
## =========================

## Log-kernel for vector theta
lf_ackley_vec <- function(theta, a = 20, b = 0.2, c = 2*pi) {
  -ackley_f(theta, a = a, b = b, c = c)
}

## Build a 1D conditional log-target for coordinate j
make_cond_logtarget <- function(theta_current, j, a = 20, b = 0.2, c = 2*pi) {
  force(theta_current); force(j)
  function(xj) {
    th <- theta_current
    th[j] <- xj
    lf_ackley_vec(th, a = a, b = b, c = c)
  }
}

## Pseudo for each coordinate: Uniform(-R, R)
R <- 5
pseudo1d <- list(
  ld = function(x) dunif(x, min = -R, max = R, log = TRUE),
  q  = function(u) qunif(u, min = -R, max = R)
)

## ---- Gibbs qslice ----
d <- 2
n_draws <- 1000

theta <- matrix(NA_real_, nrow = n_draws, ncol = d)
theta[1, ] <- rep(0, d)

nEvaluations <- 0L
time_qslice_gibbs <- system.time(
  for (i in 2:n_draws) {
    th <- theta[i - 1, ]
    
    for (j in 1:d) {
      logt_j <- make_cond_logtarget(th, j)
      
      out <- slice_quantile(
        x          = th[j],
        log_target = logt_j,
        pseudo     = pseudo1d
      )
      
      th[j] <- out$x
      nEvaluations <- nEvaluations + out$nEvaluations
    }
    
    theta[i, ] <- th
  }
)

time_qslice_gibbs
nEvaluations

## Quick diagnostics for coordinate 1
plot(theta[,2], type="l", main="Ackley dD Gibbs-qSlice: coord 1 trace", ylab=expression(theta[1]))
acf(theta[,1], main="ACF coord 1")

multiESS(theta)


### worst case scenario
set.seed(123)

## Log-kernel for vector theta
lf_ackley_vec <- function(theta, a = 20, b = 0.2, c = 2*pi) {
  -ackley_f(theta, a = a, b = b, c = c)
}

## Build a 1D conditional log-target for coordinate j
make_cond_logtarget <- function(theta_current, j, a = 20, b = 0.2, c = 2*pi) {
  force(theta_current); force(j)
  function(xj) {
    th <- theta_current
    th[j] <- xj
    lf_ackley_vec(th, a = a, b = b, c = c)
  }
}

## ---- Gibbs qslice ----
d <- 2
n_draws <- 1000

theta <- matrix(NA_real_, nrow = n_draws, ncol = d)
theta[1, ] <- rep(0, d)

nEvaluations <- 0L

# Make exploration terrible
R_global <- 0.05     # very tight global box (kills big moves)
delta    <- 1e-6     # extremely tiny local proposal width (kills local moves)

time_qslice_gibbs <- system.time(
  for (i in 2:n_draws) {
    th <- theta[i - 1, ]
    
    for (j in 1:d) {
      
      # conditional log-target at current th
      logt_j <- make_cond_logtarget(th, j)
      
      # WORST pseudo: only propose inside [x - delta, x + delta],
      # and also clip to a tiny global region [-R_global, R_global]
      lo <- max(-R_global, th[j] - delta)
      hi <- min( R_global, th[j] + delta)
      
      pseudo1d_bad <- list(
        ld = function(x) dunif(x, min = lo, max = hi, log = TRUE),
        q  = function(u) qunif(u, min = lo, max = hi)
      )
      
      out <- slice_quantile(
        x          = th[j],
        log_target = logt_j,
        pseudo     = pseudo1d_bad
      )
      
      th[j] <- out$x
      nEvaluations <- nEvaluations + out$nEvaluations
    }
    
    theta[i, ] <- th
  }
)

time_qslice_gibbs
nEvaluations

# Diagnostics: should show extremely slow decay / near-flat trace
par(mfrow=c(2,2))
plot(theta[,1], type="l", main="Trace: dim 1 (should be sticky)", xlab="iter", ylab="x1")
plot(theta[,2], type="l", main="Trace: dim 2 (should be sticky)", xlab="iter", ylab="x2")
acf(theta[,1], main="ACF dim 1 (slow decay expected)")
acf(theta[,2], main="ACF dim 2 (slow decay expected)")
par(mfrow=c(1,1))

effectiveSize(theta)

##### Banana Kernel ########

############################################################
## Banana-shaped Rosenbrock kernel
## Gibbs Quantile Slice Sampling (qslice)
## Single self-contained R script
############################################################

library(qslice)
library(coda)

## ---------------------------------------------------------
## 1. Banana (Rosenbrock-type) log-kernel
## ---------------------------------------------------------
banana_log_kernel <- function(theta) {
  x1 <- theta[1]
  x2 <- theta[2]
  -x1^2 / 10 - x2^2 / 10 - 2 * (x2 - x1^2)^2
}

## Conditional log-kernels
lf_x1_given_x2 <- function(x1, x2_fixed) {
  banana_log_kernel(c(x1, x2_fixed))
}

lf_x2_given_x1 <- function(x2, x1_fixed) {
  banana_log_kernel(c(x1_fixed, x2))
}

## ---------------------------------------------------------
## 2. Pseudo-targets (Uniform boxes)
## ---------------------------------------------------------
pseudo_x1 <- list(
  ld = function(x) dunif(x, min = -5, max = 5, log = TRUE),
  q  = function(u) qunif(u, min = -5, max = 5)
)

pseudo_x2 <- list(
  ld = function(x) dunif(x, min = -8, max = 8, log = TRUE),
  q  = function(u) qunif(u, min = -8, max = 8)
)

## --------------------------------------------------------- qslice sampler
## ---------------------------------------------------------
set.seed(123)

n_draws <- 1000
theta <- matrix(NA_real_, nrow = n_draws, ncol = 2)
colnames(theta) <- c("x1", "x2")

theta[1, ] <- c(0, 0)   # start near the mode

nEvaluations <- 0L

time_qslice <- system.time(
  for (i in 2:n_draws) {
    
    x1_prev <- theta[i - 1, 1]
    x2_prev <- theta[i - 1, 2]
    
    ## Update x1 | x2
    out1 <- slice_quantile(
      x          = x1_prev,
      log_target = function(x) lf_x1_given_x2(x, x2_prev),
      pseudo     = pseudo_x1
    )
    x1_new <- out1$x
    nEvaluations <- nEvaluations + out1$nEvaluations
    
    ## Update x2 | x1
    out2 <- slice_quantile(
      x          = x2_prev,
      log_target = function(x) lf_x2_given_x1(x, x1_new),
      pseudo     = pseudo_x2
    )
    x2_new <- out2$x
    nEvaluations <- nEvaluations + out2$nEvaluations
    
    theta[i, ] <- c(x1_new, x2_new)
  }
)

## ---------------------------------------------------------
## 4. Output summary
## ---------------------------------------------------------
print(time_qslice)
cat("Total target evaluations:", nEvaluations, "\n")

cat("Posterior means:\n")
print(colMeans(theta))

cat("Effective Sample Sizes:\n")
print(multiESS(theta))

## ---------------------------------------------------------
## 5. Diagnostics
## ---------------------------------------------------------
par(mfrow = c(2, 2))

## Trace plots
plot(theta[,1], type = "l", main = expression("Trace: " * x[1]), ylab = "")
plot(theta[,2], type = "l", main = expression("Trace: " * x[2]), ylab = "")

## ACF plots
acf(theta[,1], main = expression("ACF: " * x[1]))
acf(theta[,2], main = expression("ACF: " * x[2]))

## ---------------------------------------------------------
## 6. Banana-shaped scatter plot
## ---------------------------------------------------------
plot(
  theta[,1], theta[,2],
  pch = 16, cex = 0.4,
  xlab = expression(x[1]),
  ylab = expression(x[2]),
  main = "Banana-shaped Rosenbrock density (qslice samples)"
)

####### Elliptical Slice ###

## Supporting normal N(mu, Sig) for elliptical slice
mu_k1  <- 0
Sig_k1 <- matrix(25, nrow = 1, ncol = 1)  # variance = 25 (sd = 5)

n_iter <- 1000
draws_k1_ess <- numeric(n_iter)
draws_k1_ess[1] <- 0
nEvaluations <- 0L

time_ess_k1 <- system.time(
  for (i in 2:n_iter) {
    out <- slice_elliptical_mv(
      x          = draws_k1_ess[i - 1],  # scalar current state
      log_target = lf_k1,
      mu         = mu_k1,
      Sig        = Sig_k1,
      is_chol    = FALSE
    )
    draws_k1_ess[i] <- out$x
    nEvaluations <- nEvaluations + out$nEvaluations
  }
)

time_ess_k1
effectiveSize(draws_k1_ess)
mean(draws_k1_ess)


### Banana Kerel

set.seed(123)

n_iter <- 1000
draws_banana_ess <- matrix(0, nrow = n_iter, ncol = 2)
draws_banana_ess[1, ] <- c(0, 0)
nEvaluations <- 0L

## Supporting normal N(mu, Sig)
mu_banana  <- c(0, 0)
Sig_banana <- diag(2)   # you can change this covariance if you want

time_ess_banana <- system.time(
  for (i in 2:n_iter) {
    out <- slice_elliptical_mv(
      x          = draws_banana_ess[i - 1, ],
      log_target = banana_log_kernel,  # <-- log-kernel directly
      mu         = mu_banana,
      Sig        = Sig_banana,
      is_chol    = FALSE
    )
    draws_banana_ess[i, ] <- out$x
    nEvaluations <- nEvaluations + out$nEvaluations
  }
)

multiESS(draws_banana_ess)
time_ess_banana


# Ackley Kernel

ackley_f <- function(theta, a = 20, b = 0.2, c = 2*pi) {
  if (is.null(dim(theta))) {
    theta <- matrix(theta, nrow = 1)
  }
  d <- ncol(theta)
  ss <- rowSums(theta^2)
  cos_mean <- rowMeans(cos(c * theta))
  f <- -a * exp(-b * sqrt(ss / d)) - exp(cos_mean) + a + exp(1)
  as.numeric(f)
}

## ---------------------------------------------------------
## Define a LOG-KERNEL from Ackley energy:
## target kernel p(x) ∝ exp( - ackley_f(x) )
ackley_log_kernel <- function(theta, a = 20, b = 0.2, c = 2*pi) {
  theta <- as.numeric(theta)
  -ackley_f(theta, a = a, b = b, c = c)
}


set.seed(123)

D <- 2                 # choose dimension
n_iter <- 1000

draws_ackley_ess <- matrix(0, nrow = n_iter, ncol = D)
draws_ackley_ess[1, ] <- rep(0, D)
nEvaluations <- 0L

## Supporting normal N(mu, Sig)
mu_ackley  <- rep(0, D)
Sig_ackley <- diag(D)   # you can change scale if you want

time_ess_ackley <- system.time(
  for (i in 2:n_iter) {
    out <- slice_elliptical_mv(
      x          = draws_ackley_ess[i - 1, ],
      log_target = function(x) ackley_log_kernel(x, a = 20, b = 0.2, c = 2*pi),
      mu         = mu_ackley,
      Sig        = Sig_ackley,
      is_chol    = FALSE
    )
    draws_ackley_ess[i, ] <- out$x
    nEvaluations <- nEvaluations + out$nEvaluations
  }
)

D <- 2                 # choose dimension
n_iter <- 1000

draws_ackley_ess <- matrix(0, nrow = n_iter, ncol = D)
draws_ackley_ess[1, ] <- rep(0, D)
nEvaluations <- 0L

## Supporting normal N(mu, Sig)
mu_ackley  <- rep(0, D)
Sig_ackley <- diag(D)   # you can change scale if you want

time_ess_ackley <- system.time(
  for (i in 2:n_iter) {
    out <- slice_elliptical_mv(
      x          = draws_ackley_ess[i - 1, ],
      log_target = function(x) ackley_log_kernel(x, a = 20, b = 0.2, c = 2*pi),
      mu         = mu_ackley,
      Sig        = Sig_ackley,
      is_chol    = FALSE
    )
    draws_ackley_ess[i, ] <- out$x
    nEvaluations <- nEvaluations + out$nEvaluations
  }
)

multiESS(draws_ackley_ess)
time_ess_ackley



### Flaws

D <- 2
mu_ackley  <- c(50, -50)  

rho <- 0.999999
s1 <- 1
s2 <- 1
Sig_ackley <- matrix(c(s1^2,       rho*s1*s2,
                       rho*s1*s2,  s2^2), nrow = 2, byrow = TRUE)


time_ess_ackley_flaw <- system.time(
  for (i in 2:n_iter) {
    out <- slice_elliptical_mv(
      x          = draws_ackley_ess[i - 1, ],
      log_target = function(x) ackley_log_kernel(x, a = 20, b = 0.2, c = 2*pi),
      mu         = mu_ackley,
      Sig        = Sig_ackley,
      is_chol    = FALSE
    )
    draws_ackley_ess[i, ] <- out$x
    nEvaluations <- nEvaluations + out$nEvaluations
  }
)

effectiveSize(draws_ackley_ess)
