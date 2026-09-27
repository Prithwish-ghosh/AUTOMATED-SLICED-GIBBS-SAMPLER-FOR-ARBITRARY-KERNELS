## ==================================================================
## dim_sweep_and_plot.R
##
## Runs ASG, LaplacesDemon-AFSS and LaplacesDemon-ESS on the same
## lattice-of-modes kernel at m = 2..5, records the % of modes each
## one actually visited, and plots those recorded values directly --
## no numbers are typed in separately from the run.
## ==================================================================
library(LaplacesDemon)
library(mcmcse)
source("effective_support_uni_s-2.R")
source("level_set_endpoints.R")
source("asg_sample bracket.R")

set.seed(2026)
vals  <- c(-6, 0, 6)
sigma <- 0.4
N <- 800; BURN <- 200; ITER <- N + BURN
dims <- c(2, 3, 4, 5)

## results accumulate here as the loop runs -- this IS the data the
## plot is built from, nothing is retyped afterward
results <- data.frame(m = integer(), n_modes = integer(),
                      pct_ASG = numeric(), pct_AFSS = numeric(), pct_ESS = numeric())

for (m in dims) {
  centers <- as.matrix(do.call(expand.grid, rep(list(vals), m)))
  n_modes <- nrow(centers)
  K <- function(x) sum(exp(-rowSums((matrix(x, nrow(centers), m, byrow=TRUE) - centers)^2)/(2*sigma^2)))
  n_visited <- function(S) {
    d2 <- apply(centers, 1, function(c) rowSums((S - matrix(c, nrow(S), m, byrow=TRUE))^2))
    length(unique(apply(d2, 1, which.min)))
  }
  init <- rep(0.01, m)
  
  ## --- ASG ---
  fit <- multivariate_gibbs_sample_ASG(K, n_samples=N, burn_in=BURN, theta_init=init,
                                       bracket="adaptive", make_plots=FALSE)
  v_asg <- n_visited(fit$samples)
  sample_Q_asg <- multiESS(fit$samples)
  
  ## --- LaplacesDemon AFSS / ESS ---
  Model <- function(parm, Data) { LP<-log(K(parm)); if(!is.finite(LP)) LP<--1e10
  list(LP=LP,Dev=-2*LP,Monitor=LP,yhat=parm,parm=parm) }
  MyData <- list(mon.names="LP", parm.names=paste0("x",1:m))
  
  fit_a <- LaplacesDemon(Model, Data=MyData, Initial.Values=init, Covar=NULL,
                         Iterations=ITER, Status=ITER+1, Thinning=1, Algorithm="AFSS",
                         Specs=list(A=BURN,B=NULL,m=100,n=0,w=1))
  v_afss <- n_visited(fit_a$Posterior1[(BURN+1):ITER,,drop=FALSE])
  sample_Q_afss <- multiESS(fit_a$Posterior1[(BURN+1):ITER,,drop=FALSE])
  
  fit_e <- LaplacesDemon(Model, Data=MyData, Initial.Values=init, Covar=diag(m)*36,
                         Iterations=ITER, Status=ITER+1, Thinning=1, Algorithm="ESS",
                         Specs=list(B=NULL))
  v_ess <- n_visited(fit_e$Posterior1[(BURN+1):ITER,,drop=FALSE])
  sample_Q_ess <- multiESS(fit_e$Posterior1[(BURN+1):ITER,,drop=FALSE])
  
  results <- rbind(results, data.frame(
    m = m, n_modes = n_modes,
    pct_ASG  = 100 * v_asg  / n_modes,
    pct_AFSS = 100 * v_afss / n_modes,
    pct_ESS  = 100 * v_ess  / n_modes))
  
  cat(sprintf("m=%d  modes=%4d | ASG %5.1f%% | AFSS %5.1f%% | ESS %5.1f%%\n",
              m, n_modes, results$pct_ASG[nrow(results)],
              results$pct_AFSS[nrow(results)], results$pct_ESS[nrow(results)]))
}
sample_Q_afss
sample_Q_asg
sample_Q_ess
write.csv(results, "dim_sweep_results.csv", row.names = FALSE)

## --- plot built directly from the `results` data frame above ---
png("dimension_curse_comparison.png", width=900, height=600, res=120)
plot(results$m, results$pct_ASG, type="b", pch=16, col="darkgreen", lwd=2, ylim=c(0,105),
     xlab="dimension m (3^m modes on a coordinate lattice)",
     ylab="% of modes visited", xaxt="n",
     main="ASG vs. LaplacesDemon AFSS/ESS on a lattice of narrow modes")
axis(1, at=results$m, labels=paste0("m=",results$m,"\n(3^",results$m,"=",results$n_modes,")"))
lines(results$m, results$pct_ESS,  type="b", pch=16, col="steelblue", lwd=2)
lines(results$m, results$pct_AFSS, type="b", pch=16, col="firebrick", lwd=2)
legend("bottomleft", c("ASG-adaptive","LD-ESS","LD-AFSS"),
       col=c("darkgreen","steelblue","firebrick"), lwd=2, pch=16, bty="n")
dev.off()

cat("\nSaved dim_sweep_results.csv and dimension_curse_comparison.png\n")
print(results, row.names = FALSE)
