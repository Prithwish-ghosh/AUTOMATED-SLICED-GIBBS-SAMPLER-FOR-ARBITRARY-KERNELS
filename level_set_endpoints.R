## ==================================================================
## level_set_endpoints.R
##
## Adaptive search for the OUTERMOST end points of a level set
##   S_u = { x : g(x) > u }
## following S. K. Ghosh, "On Finding Level Set End Points":
## directional stepping (secant-type adaptive step) + bisection.
##
## Illustration / GitHub only -- not part of the paper's main method.
## In ASG it can be used as an optional alternative to the grid-scan
## check extreme_roots_bracket() in asg_sampler.R (bracket = "adaptive").
##
## Changes relative to the note's code (see comments marked [CHANGE]):
##   * slice is g > u (strict), as in the paper
##   * bisection returns the OUTSIDE end of each bracket, so [a, b]
##     never cuts off slice mass (up to tol)
##   * non-finite g values (singularities) are treated safely
##   * reports 'hit_edge' when the slice reaches the end of (0, 1),
##     i.e. the level set is unbounded on the original scale
##   * max_iter guard on the stepping loops
## ==================================================================


## ------------------------------------------------------------------
## 1. Core routine on (0, 1)
## g  : function on (0, 1);  u : slice height;  x0 : point with g(x0) > u
## ------------------------------------------------------------------
find_outer_boundaries <- function(g, u, x0, min_step = 1e-4, max_step = 0.05,
                                  tol = 1e-7, eps = 1e-9, max_iter = 1e5L) {
  safety_factor <- 0.5                      # take 50% of estimated distance
  inside <- function(x) {                   # [CHANGE] strict, NaN-safe
    gx <- g(x)
    isTRUE(gx > u)                          # Inf counts as inside, NaN as outside
  }
  gval <- function(x) { gx <- g(x); if (is.finite(gx)) gx else if (isTRUE(gx > 0)) Inf else 0 }
  if (!inside(x0)) stop("x0 must satisfy g(x0) > u")

  ## one directional search: dir = +1 (from the left end) or -1 (from the right end)
  search_side <- function(dir) {
    start <- if (dir > 0) eps else 1 - eps
    if (inside(start)) return(list(end = start, hit_edge = TRUE))

    curr <- start; g_curr <- gval(curr)
    step <- max_step / 2
    nxt <- curr
    for (it in seq_len(max_iter)) {
      nxt <- if (dir > 0) min(curr + step, x0) else max(curr - step, x0)
      if (inside(nxt)) break                # bracketed: curr outside, nxt inside
      g_nxt <- gval(nxt)
      delta_g <- g_nxt - g_curr
      if (is.finite(delta_g) && delta_g > 0) {
        next_step <- step * ((u - g_nxt) / delta_g) * safety_factor   # climbing
      } else {
        next_step <- step * 1.5                                        # falling / flat
      }
      step <- min(max(next_step, min_step), max_step)
      curr <- nxt; g_curr <- g_nxt
    }
    ## bisection on the bracket, keeping out = outside, inn = inside
    out <- curr; inn <- nxt
    while (abs(inn - out) > tol) {
      mid <- (out + inn) / 2
      if (inside(mid)) inn <- mid else out <- mid
    }
    list(end = out, hit_edge = FALSE)       # [CHANGE] outside end
  }

  L <- search_side(+1); R <- search_side(-1)
  list(a = L$end, b = R$end, hit_edge = L$hit_edge || R$hit_edge)
}


## ------------------------------------------------------------------
## 2. ASG wrapper: search on the Cauchy scale v = F_C(x; scale0)
## g : conditional kernel on the natural scale (as in asg_sampler.R)
## Returns the same fields as extreme_roots_bracket(): lower, upper, ok
## ------------------------------------------------------------------
## NOTE on max_step: any slice component narrower than max_step (on the
## Cauchy scale) can be stepped over, and then [a, b] misses it without
## any warning -- the retain rule does NOT protect against this. Slice
## components shrink as u grows, so max_step must be small. On the Beta
## mixture of the paper, max_step = 0.05 or 0.01 visibly biases the
## chain (a mode is under-visited); 0.002 matches the grid check.
adaptive_outer_bracket <- function(g, u, x0, scale0 = 1, min_step = 1e-4,
                                   max_step = 0.002, tol = 1e-8) {
  G  <- function(v) g(qcauchy(v, scale = scale0))
  v0 <- pcauchy(x0, scale = scale0)
  if (!isTRUE(G(v0) > u))
    return(list(lower = NA_real_, upper = NA_real_, ok = FALSE))
  br <- find_outer_boundaries(G, u, v0, min_step = min_step,
                              max_step = max_step, tol = tol)
  a <- qcauchy(br$a, scale = scale0); b <- qcauchy(br$b, scale = scale0)
  list(lower = a, upper = b,
       ok = is.finite(a) && is.finite(b) && a < x0 && x0 < b && !br$hit_edge)
}


## ------------------------------------------------------------------
## 3. Drop-in replacement for asg_coordinate_update() in asg_sampler.R
##    bracket = "box"      -> effective-support box (paper's main method)
##              "grid"     -> extreme_roots_bracket() (grid scan + bisection)
##              "adaptive" -> adaptive_outer_bracket() (this file)
##    For backward compatibility, check_bracket = TRUE/FALSE still works
##    and maps to "grid"/"box".
##    Requires effective.support(), extreme_roots_bracket() and
##    slice1d_fixed_u_x0() from asg_sampler.R.
## ------------------------------------------------------------------
asg_coordinate_update <- function(g, u, x0, check_bracket = FALSE, box = NULL,
                                  tol = 0.01, scale0 = 1,
                                  n_grid = 1000L, bisect_e = 1e-8,
                                  bracket = if (check_bracket) "grid" else "box",
                                  min_step = 1e-4, max_step = 0.002) {
  bracket <- match.arg(bracket, c("box", "grid", "adaptive"))
  if (bracket == "box") {
    if (is.null(box)) box <- effective.support(g, tol = tol, scale0 = scale0)
    a <- box$lower; b <- box$upper
  } else {
    br <- if (bracket == "grid") {
      extreme_roots_bracket(g, u, x0, scale0 = scale0, n_grid = n_grid, e = bisect_e)
    } else {
      adaptive_outer_bracket(g, u, x0, scale0 = scale0, min_step = min_step,
                             max_step = max_step, tol = bisect_e)
    }
    if (!br$ok) return(list(x = x0, retained = TRUE))
    a <- br$lower; b <- br$upper
  }
  slice1d_fixed_u_x0(g, u, a, b, x0)
}
## To use it inside dim1_gibbs_sample_ASG / multivariate_gibbs_sample_ASG,
## pass `bracket = bracket` in their asg_coordinate_update() calls and add
## `bracket = "box"` to their argument lists.


## ==================================================================
## 4. Demo: the three test functions from the note
## [CHANGE] u is set BELOW g(x0) (u = g(x0) - 0.02), as the method
## requires g(x0) > u -- in ASG this always holds since u ~ Unif(0, K(x)).
## [CHANGE] the dotted line marks x0 (the note plotted g(x0) on the x-axis).
## ==================================================================
if (sys.nframe() == 0L) {
  tests <- list(
    "Multimodal"                   = list(g = function(x) x * sin(5 * pi * x)^2 + 0.5,
                                          x0 = c(0.05, 0.35, 0.65, 0.87)),
    "Multimodal, singularity at 0" = list(g = function(x) 1 / x^0.1 + 10 * x * sin(10 * pi * x)^2,
                                          x0 = c(0.05, 0.35, 0.65, 0.85)),
    "Multimodal, singularity at 1" = list(g = function(x) 1 / (1 - x)^0.1 +
                                            10 * (1 - x) * sin(10 * pi * x)^2,
                                          x0 = c(0.03, 0.45, 0.65, 0.85))
  )

  ## brute-force check: true outermost end points on a fine grid
  ## (log-spaced near 0 and 1 so that thin pieces at a singularity are seen)
  true_ends <- function(g, u) {
    edge <- 10^seq(-9, -3, length.out = 2000)
    xs <- sort(c(edge, seq(1e-3, 1 - 1e-3, length.out = 2e5), 1 - edge))
    ins <- which(g(xs) > u)
    c(xs[min(ins)], xs[max(ins)])
  }

  pdf("level_set_endpoints_demo.pdf", width = 9, height = 7)
  for (nm in names(tests)) {
    g <- tests[[nm]]$g
    par(mfrow = c(2, 2), lwd = 2, mar = c(4, 4, 2.5, 1))
    cat("\n==", nm, "==\n")
    for (x0 in tests[[nm]]$x0) {
      u  <- g(x0) - 0.02
      br <- find_outer_boundaries(g, u, x0)
      te <- true_ends(g, u)
      cat(sprintf("x0 = %.2f  u = %.3f  a = %.6f  b = %.6f  (grid: %.6f, %.6f)%s\n",
                  x0, u, br$a, br$b, te[1], te[2],
                  if (br$hit_edge) "  [reaches edge]" else ""))
      curve(g, 1e-4, 1 - 1e-4, n = 5000, ylab = "g(x)",
            main = sprintf("%s, x0 = %.2f", nm, x0), cex.main = 0.85)
      abline(h = u, col = "blue")
      abline(v = x0, col = "grey50", lty = 3)
      abline(v = c(br$a, br$b), col = "red", lty = 2)
    }
  }
  dev.off()
  cat("\nPlots written to level_set_endpoints_demo.pdf\n")
}
