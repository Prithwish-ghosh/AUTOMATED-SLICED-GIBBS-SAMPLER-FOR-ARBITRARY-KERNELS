library(LaplacesDemon)
source("effective_support_uni_s-2.R")
source("level_set_endpoints.R")
source("asg_sample bracket.R")

set.seed(2026)

## 3x3 grid of narrow, equal-weight Gaussian bumps at all combinations
## of {-6, 0, 6} x {-6, 0, 6}. Population covariance of this target is
## EXACTLY diagonal (X1, X2 independent uniform over {-6,0,6}), so PCA
## of the full grid gives the coordinate axes -- but any method whose
## exploration is confined to fewer than both axes (e.g. only the two
## diagonals) structurally cannot reach the 4 "edge" points.
centers <- expand.grid(c1 = c(-6,0,6), c2 = c(-6,0,6))
sigma <- 0.4
K <- function(x) sum(exp(-((x[1]-centers$c1)^2 + (x[2]-centers$c2)^2)/(2*sigma^2)))

nearest_cell <- function(S) {
  S <- as.matrix(S); colnames(S) <- c("c1","c2")
  C <- as.matrix(centers)
  d2 <- outer(1:nrow(S), 1:nrow(C), Vectorize(function(i,j) sum((S[i,]-C[j,])^2)))
  apply(d2, 1, which.min)
}
n_visited <- function(S) length(unique(nearest_cell(S)))

N <- 1500; BURN <- 300; ITER <- N + BURN

cat("=== 3x3 grid of narrow Gaussian bumps (9 modes, spacing 6, sigma 0.4) ===\n\n")

## ASG
for (br in c("box","grid","adaptive")) {
  t0 <- Sys.time()
  fit <- multivariate_gibbs_sample_ASG(K, n_samples=N, burn_in=BURN,
                                       theta_init=c(0,0), bracket=br, make_plots=FALSE)
  el <- as.numeric(Sys.time()-t0, units="secs")
  cat(sprintf("ASG-%-8s  modes visited: %d/9   time=%.1fs\n", br, n_visited(fit$samples), el))
}

## LaplacesDemon AFSS and ESS
Model <- function(parm, Data) {
  LP <- log(K(parm)); if (!is.finite(LP)) LP <- -1e10
  list(LP=LP, Dev=-2*LP, Monitor=LP, yhat=parm, parm=parm)
}
MyData <- list(mon.names="LP", parm.names=c("x1","x2"))

t0 <- Sys.time()
fit_afss <- LaplacesDemon(Model, Data=MyData, Initial.Values=c(6,6),
                          Covar=NULL, Iterations=ITER, Status=ITER+1, Thinning=1,
                          Algorithm="AFSS", Specs=list(A=BURN, B=NULL, m=100, n=0, w=1))
el <- as.numeric(Sys.time()-t0, units="secs")
S <- fit_afss$Posterior1[(BURN+1):ITER,,drop=FALSE]
cat(sprintf("LD-AFSS      modes visited: %d/9   time=%.1fs\n", n_visited(S), el))

t0 <- Sys.time()
fit_ess <- LaplacesDemon(Model, Data=MyData, Initial.Values=c(6,6),
                         Covar=diag(2)*36, Iterations=ITER, Status=ITER+1, Thinning=1,
                         Algorithm="ESS", Specs=list(B=NULL))
el <- as.numeric(Sys.time()-t0, units="secs")
S <- fit_ess$Posterior1[(BURN+1):ITER,,drop=FALSE]
cat(sprintf("LD-ESS       modes visited: %d/9   time=%.1fs\n", n_visited(S), el))

png("grid_kernel_comparison.png", width=1400, height=500, res=120)
par(mfrow=c(1,3), mar=c(4,4,3,1))
xs <- seq(-9,9,length.out=250)
Z <- outer(xs,xs,Vectorize(function(a,b) K(c(a,b))))
fitb <- multivariate_gibbs_sample_ASG(K, n_samples=N, burn_in=BURN, theta_init=c(0,0), bracket="grid", make_plots=FALSE)
plots <- list("ASG-grid"=fitb$samples,
              "LD-AFSS"=fit_afss$Posterior1[(BURN+1):ITER,],
              "LD-ESS"=fit_ess$Posterior1[(BURN+1):ITER,])
for (nm in names(plots)) {
  contour(xs,xs,Z, drawlabels=FALSE, col="grey70", main=paste0(nm," (",n_visited(plots[[nm]]),"/9 modes)"), xlab="x1", ylab="x2")
  points(plots[[nm]], pch=16, cex=0.3, col=adjustcolor("red",0.4))
}
library(mcmcse)
multiESS(fitb$samples)
multiESS(fit_afss$Posterior1[(BURN+1):ITER,])
multiESS(fit_ess$Posterior1[(BURN+1):ITER,])





dev.off()
