rm(list = ls())

# Function to compute running mean
running_mean <- function(x) {
  cumsum(x) / seq_along(x)
}

# Function to assess mixing quality
assess_mixing <- function(accept_rate, ess, acf_vals, trace, is_multimodal = FALSE) {
  cat("Mixing Assessment:\n")
  if (accept_rate >= 0.2 && accept_rate <= 0.5) {
    cat("- Acceptance rate:", accept_rate, "(Good: within 0.2–0.5)\n")
  } else if (accept_rate < 0.2) {
    cat("- Acceptance rate:", accept_rate, "(Low: consider reducing scale)\n")
  } else {
    cat("- Acceptance rate:", accept_rate, "(High: consider increasing scale)\n")
  }
  if (is.vector(ess)) {
    for (i in seq_along(ess)) {
      if (ess[i] > 500) {
        cat("- ESS for dimension", i, ":", ess[i], "(Good: >500)\n")
      } else if (ess[i] > 100) {
        cat("- ESS for dimension", i, ":", ess[i], "(Moderate: 100–500)\n")
      } else {
        cat("- ESS for dimension", i, ":", ess[i], "(Poor: <100, consider more iterations or adjusting scale)\n")
      }
    }
  } else {
    if (ess > 500) {
      cat("- ESS:", ess, "(Good: >500)\n")
    } else if (ess > 100) {
      cat("- ESS:", ess, "(Moderate: 100–500)\n")
    } else {
      cat("- ESS:", ess, "(Poor: <100, consider more iterations or adjusting scale)\n")
    }
  }
  acf_lag10 <- acf_vals$acf[11]
  if (is.vector(acf_lag10)) {
    for (i in seq_along(acf_lag10)) {
      if (abs(acf_lag10[i]) < 0.1) {
        cat("- ACF at lag 10 for dimension", i, ":", acf_lag10[i], "(Good: near 0)\n")
      } else {
        cat("- ACF at lag 10 for dimension", i, ":", acf_lag10[i], "(High: slow decay, consider more )\n")
      }
    }
  } else {
    if (abs(acf_lag10) < 0.1) {
      cat("- ACF at lag 10:", acf_lag10, "(Good: near 0)\n")
    } else {
      cat("- ACF at lag 10:", acf_lag10, "(High: slow decay, consider more)\n")
    }
  }
  if (is_multimodal) {
    cat("- Trace: Check plot for exploration of all modes (multimodal distribution)\n")
  } else {
    cat("- Trace: Check plot for no obvious trends or sticking\n")
  }
  cat("\n")
}

