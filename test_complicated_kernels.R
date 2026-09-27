## ==================================================================
## test_complicated_kernels.R
##
## Stress tests for ASG on harder kernels, comparing the three bracket
## options:  "box" (effective support, paper's main method),
##           "grid" (extreme-roots grid scan), "adaptive" (level-set
##           end points by adaptive stepping, level_set_endpoints.R).
## Every kernel has a known truth, computed independently of ASG, so
## accuracy is checked directly and not just by ESS.
##
## Usage:  Rscript test_complicated_kernels.R
## Output: asg_test_results.csv, asg_test_plots.pdf, printed summary
## ==================================================================

source("effective_support_uni_s.R")
source("level_set_endpoints.R")
source("asg_sampler.R")

set.seed(2026)
N1   <- as.integer(Sys.getenv("N1", 3000))   # retained samples, 1-D kernels
N2   <- as.integer(Sys.getenv("N2", 1500))   # retained samples, 2-D kernels
BURN <- as.integer(Sys.getenv("BURN", 300))
BRACKETS <- c("box", "grid", "adaptive")

ess1 <- function(x) as.numeric(coda::effectiveSize(x))

## ------------------------------------------------------------------
## 1-D kernels: each has an exact CDF and known mode weights.
## ------------------------------------------------------------------
tests1 <- list(
  beta_mixture = list(
    desc = "Beta mixture (paper, eq. 11): unbounded spikes, gaps in support",
    K = function(x) 0.3/2*dbeta((x+5)/2, .5, 1) + 0.4/2*dbeta((x+1)/2, 2, 2) +
                    0.3/2*dbeta((x-5)/2, 1, .5),
    F = function(x) 0.3*pbeta((x+5)/2, .5, 1) + 0.4*pbeta((x+1)/2, 2, 2) +
                    0.3*pbeta((x-5)/2, 1, .5),
    init = 0, cuts = c(-2, 2), w = c(0.3, 0.4, 0.3)),

  five_narrow_modes = list(
    desc = "Five narrow, far-apart Gaussian modes with unequal weights",
    K = function(x) 0.10*dnorm(x, -15, 0.3) + 0.30*dnorm(x, -4, 0.5) +
                    0.20*dnorm(x, 0, 0.1) + 0.25*dnorm(x, 6, 1) +
                    0.15*dnorm(x, 25, 0.2),
    F = function(x) 0.10*pnorm(x, -15, 0.3) + 0.30*pnorm(x, -4, 0.5) +
                    0.20*pnorm(x, 0, 0.1) + 0.25*pnorm(x, 6, 1) +
                    0.15*pnorm(x, 25, 0.2),
    init = 0, cuts = c(-10, -2, 2, 15), w = c(0.10, 0.30, 0.20, 0.25, 0.15)),

  cauchy_plus_spike = list(
    desc = "Heavy-tailed Cauchy bulk plus a very narrow distant spike",
    K = function(x) 0.7*dcauchy(x) + 0.3*dnorm(x, 10, 0.05),
    F = function(x) 0.7*pcauchy(x) + 0.3*pnorm(x, 10, 0.05),
    init = 0, cuts = 9.5,
    w = c(0.7*pcauchy(9.5), 0.7*(1 - pcauchy(9.5)) + 0.3)),

  oscillating = list(
    desc = "Gaussian envelope times sin^2: many modes separated by exact zeros",
    K = function(x) exp(-x^2/50) * sin(2*x)^2,
    F = NULL,        # computed numerically below
    init = 0.8, cuts = 0, w = c(0.5, 0.5))
)
local({
  xs <- seq(-40, 40, length.out = 400001)
  k  <- tests1$oscillating$K(xs)
  cdf <- cumsum(k); cdf <- cdf / cdf[length(cdf)]
  tests1$oscillating$F <<- approxfun(xs, cdf, yleft = 0, yright = 1)
})

## ------------------------------------------------------------------
## 2-D kernels: truth by numerical integration on a fine grid
## ------------------------------------------------------------------
tests2 <- list(
  rosenbrock = list(
    desc = "Rosenbrock banana (paper, eq. 12)",
    K = function(x) exp(-x[1]^2/10 - x[2]^2/10 - 2*(x[2] - x[1]^2)^2),
    init = c(0, 0), lim = list(c(-8, 8), c(-6, 20)),
    region = function(x1, x2) x1 > 0),

  bimodal_axis = list(
    desc = "Two well-separated modes along an axis, unequal weights",
    K = function(x) 0.3*exp(-((x[1]+4)^2 + x[2]^2)/(2*0.5^2)) +
                    0.7*exp(-((x[1]-4)^2 + x[2]^2)/(2*0.5^2)),
    init = c(-4, 0), lim = list(c(-8, 8), c(-4, 4)),
    region = function(x1, x2) x1 > 0),

  bimodal_diagonal = list(
    desc = "Two modes on a diagonal: hard for any coordinate-wise sampler",
    K = function(x) 0.5*exp(-((x[1]+3)^2 + (x[2]+3)^2)/(2*0.7^2)) +
                    0.5*exp(-((x[1]-3)^2 + (x[2]-3)^2)/(2*0.7^2)),
    init = c(-3, -3), lim = list(c(-7, 7), c(-7, 7)),
    region = function(x1, x2) x1 + x2 > 0),

  ring = list(
    desc = "Thin ring (radius 3, width 0.2): conditionals have two pieces",
    K = function(x) exp(-(sqrt(x[1]^2 + x[2]^2) - 3)^2 / (2*0.2^2)),
    init = c(3, 0), lim = list(c(-4.5, 4.5), c(-4.5, 4.5)),
    region = function(x1, x2) x1 > 0)
)

truth2 <- function(tt, n = 1201) {
  x1 <- seq(tt$lim[[1]][1], tt$lim[[1]][2], length.out = n)
  x2 <- seq(tt$lim[[2]][1], tt$lim[[2]][2], length.out = n)
  G  <- expand.grid(x1 = x1, x2 = x2)
  k  <- apply(G, 1, tt$K); w <- k / sum(k)
  c(m1 = sum(w*G$x1), m2 = sum(w*G$x2),
    s1 = sqrt(sum(w*G$x1^2) - sum(w*G$x1)^2),
    s2 = sqrt(sum(w*G$x2^2) - sum(w*G$x2)^2),
    p  = sum(w[tt$region(G$x1, G$x2)]))
}

## ------------------------------------------------------------------
## Run
## ------------------------------------------------------------------
res  <- list()
keep <- list()

for (nm in names(tests1)) {
  tt <- tests1[[nm]]
  for (br in BRACKETS) {
    fit <- tryCatch(
      dim1_gibbs_sample_ASG(tt$K, n_samples = N1, burn_in = BURN,
                            theta_init = tt$init, bracket = br,
                            check_bracket = FALSE, make_plots = FALSE),
      error = function(e) { message(nm, "/", br, ": ", conditionMessage(e)); NULL })
    if (is.null(fit)) next
    s  <- as.numeric(fit$samples)
    ks <- max(abs(ecdf(s)(sort(s)) - tt$F(sort(s))))
    wh <- as.numeric(table(factor(findInterval(s, tt$cuts),
                                  levels = seq_along(tt$w) - 1))) / length(s)
    res[[length(res) + 1]] <- data.frame(
      kernel = nm, dim = 1, bracket = br,
      KS = round(ks, 3),
      max_weight_err = round(max(abs(wh - tt$w)), 3),
      ESS = round(ess1(s)), time_s = round(fit$time_taken, 1),
      ESS_per_s = round(ess1(s) / fit$time_taken, 1),
      retain = round(fit$retain_rate, 4))
    keep[[paste(nm, br)]] <- s
    cat(sprintf("%-18s %-9s KS=%.3f  weight err=%.3f  ESS=%5.0f  %5.1fs\n",
                nm, br, ks, max(abs(wh - tt$w)), ess1(s), fit$time_taken))
  }
}

for (nm in names(tests2)) {
  tt <- tests2[[nm]]; tr <- truth2(tt)
  for (br in BRACKETS) {
    fit <- tryCatch(
      multivariate_gibbs_sample_ASG(tt$K, n_samples = N2, burn_in = BURN,
                                    theta_init = tt$init, bracket = br,
                                    make_plots = FALSE),
      error = function(e) { message(nm, "/", br, ": ", conditionMessage(e)); NULL })
    if (is.null(fit)) next
    S  <- fit$samples
    est <- c(mean(S[, 1]), mean(S[, 2]), sd(S[, 1]), sd(S[, 2]),
             mean(tt$region(S[, 1], S[, 2])))
    ess <- min(ess1(S[, 1]), ess1(S[, 2]))
    res[[length(res) + 1]] <- data.frame(
      kernel = nm, dim = 2, bracket = br,
      KS = NA,
      max_weight_err = round(abs(est[5] - tr["p"]), 3),
      ESS = round(ess), time_s = round(fit$time_taken, 1),
      ESS_per_s = round(ess / fit$time_taken, 1),
      retain = round(fit$retain_rate, 4))
    keep[[paste(nm, br)]] <- S
    cat(sprintf("%-18s %-9s mean err=%.3f  sd err=%.3f  region err=%.3f  minESS=%5.0f  %5.1fs\n",
                nm, br, max(abs(est[1:2] - tr[1:2])), max(abs(est[3:4] - tr[3:4])),
                abs(est[5] - tr["p"]), ess, fit$time_taken))
  }
}

res <- do.call(rbind, res)
write.csv(res, "asg_test_results.csv", row.names = FALSE)

## ------------------------------------------------------------------
## Plots
## ------------------------------------------------------------------
pdf("asg_test_plots.pdf", width = 11, height = 3.8)
for (nm in names(tests1)) {
  tt <- tests1[[nm]]
  par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
  for (br in BRACKETS) {
    s <- keep[[paste(nm, br)]]; if (is.null(s)) { plot.new(); next }
    lo <- quantile(s, 0.001); hi <- quantile(s, 0.999)
    xs <- seq(lo, hi, length.out = 3000)
    hist(s[s >= lo & s <= hi], breaks = 120, freq = FALSE, col = "skyblue",
         border = NA, main = paste0(nm, " - ", br), xlab = "x", cex.main = 0.95)
    dens <- diff(tt$F(xs)) / diff(xs) / (tt$F(hi) - tt$F(lo))
    lines(xs[-1], dens, col = "red", lwd = 1.5)
  }
}
for (nm in names(tests2)) {
  tt <- tests2[[nm]]
  x1 <- seq(tt$lim[[1]][1], tt$lim[[1]][2], length.out = 200)
  x2 <- seq(tt$lim[[2]][1], tt$lim[[2]][2], length.out = 200)
  Z  <- outer(x1, x2, Vectorize(function(a, b) tt$K(c(a, b))))
  par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
  for (br in BRACKETS) {
    S <- keep[[paste(nm, br)]]; if (is.null(S)) { plot.new(); next }
    contour(x1, x2, Z, drawlabels = FALSE, col = "grey60",
            main = paste0(nm, " - ", br), xlab = "x1", ylab = "x2", cex.main = 0.95)
    points(S, pch = 16, cex = 0.25, col = adjustcolor("red", 0.4))
  }
}

cat("\nSaved asg_test_results.csv and asg_test_plots.pdf\n")
print(res, row.names = FALSE)
