# ── Internal utilities (not exported) ─────────────────────────────────────────

#' Build the pre-factorized Cholesky of the HP system matrix
#'
#' Constructs the second-difference matrix `D` and factorizes the HP system
#' `Q = I + lambda * D'D` once, so the (expensive) factorization can be reused
#' across every bootstrap replicate and every bHP inner iteration instead of
#' being recomputed on each [Matrix::solve()] call.
#'
#' @param n Integer length of the series (must be >= 3).
#' @param lambda Numeric smoothing parameter.
#' @return A sparse Cholesky factor (from [Matrix::Cholesky()]).
#' @keywords internal
#' @noRd
.build_hp_cholesky <- function(n, lambda) {
  D <- Matrix::bandSparse(
    n - 2L, n,
    k         = c(0L, 1L, 2L),
    diagonals = list(rep(1, n - 2L), rep(-2, n - 2L), rep(1, n - 2L))
  )
  Q <- Matrix::Diagonal(n) + lambda * Matrix::crossprod(D)
  Matrix::Cholesky(Matrix::forceSymmetric(Q))
}

#' Fast HP cycle for internal calibration and bootstrap
#'
#' Returns only the cycle (`y - trend`) of the HP filter. No input validation,
#' no class wrapping, no timing. Pass a pre-built Cholesky factor `CH` to skip
#' the per-call factorization (the hot path inside `.boot_engine`); otherwise
#' the factor is built from `lambda`, which itself falls back to the
#' Ravn-Uhlig rule via `stats::frequency(x)`.
#'
#' @param x Numeric vector or `ts` object.
#' @param lambda Optional smoothing parameter (used only when `CH` is `NULL`).
#' @param CH Optional pre-built Cholesky factor from `.build_hp_cholesky()`.
#' @return Numeric vector of cyclical residuals, same length as `x`.
#' @keywords internal
#' @noRd
.hp_fast <- function(x, lambda = NULL, CH = NULL) {
  n <- length(x)
  if (n < 3L) return(rep(0, n))          # bandSparse needs n-2 >= 1 rows
  y <- as.double(x)

  if (is.null(CH)) {
    if (is.null(lambda)) {
      freq   <- stats::frequency(x)      # 1 for numeric/annual, 4/12 for ts
      lambda <- 6.25 * freq^4            # Ravn-Uhlig rule
    }
    CH <- .build_hp_cholesky(n, lambda)
  }

  y - as.numeric(Matrix::solve(CH, y))
}


#' Fast Hamilton cycle (regression residuals)
#'
#' Bare-metal version of [hamilton_filter()] for use inside `.boot_engine`.
#' Builds the lag matrix with `embed()` and solves the OLS via `lm.fit()`;
#' returns only the cycle, full length, with the first `h + p - 1` lead-in
#' positions left as `NA` (package convention).
#'
#' @param y Numeric vector.
#' @param h Integer horizon.
#' @param p Integer number of lags.
#' @return Numeric cycle vector, same length as `y` (NA lead-in).
#' @keywords internal
#' @noRd
.hamilton_fast <- function(y, h, p) {
  n       <- length(y)
  n_valid <- n - h - p + 1L
  E       <- stats::embed(y, p)
  X       <- cbind(1, E[seq_len(n_valid), , drop = FALSE])
  target  <- y[(p + h):n]
  fit     <- stats::lm.fit(X, target)

  cycle <- rep(NA_real_, n)
  cycle[(p + h):n] <- fit$residuals
  cycle
}


#' Fast boosted-HP cycle (fixed iteration count)
#'
#' Bare-metal version of [bhp_filter()] for use inside `.boot_engine`. Runs a
#' tight `for` loop that re-applies the HP smoother to the running residual a
#' fixed number of times (no BIC/ADF stopping) and returns the accumulated
#' cycle. Conditioning on a fixed `iter` is what makes the bootstrap refit use
#' the *same estimator* as the base fit.
#'
#' @param y Numeric vector (bootstrap draw).
#' @param CH Pre-built Cholesky factor from `.build_hp_cholesky()`.
#' @param iter Integer number of boosting passes (the base fit's final count).
#' @return Numeric cycle vector, same length as `y`.
#' @keywords internal
#' @noRd
.bhp_fast <- function(y, CH, iter) {
  # Boosted HP re-smooths the RESIDUAL each pass: the cycle of the cycle.
  # Iterating `.hp_fast` (which returns the cycle) reproduces the base loop's
  # `u <- u - HP_trend(u)` exactly, while reusing the cached factor CH.
  res <- as.double(y)
  for (i in seq_len(iter)) {
    res <- .hp_fast(res, CH = CH)
  }
  res                                          # final cycle after `iter` passes
}


#' Resolve the bootstrap block size
#'
#' Shared "auto" logic so every public filter resolves `block_size`
#' identically: two full cycles (`2 * frequency`) for `"auto"`, otherwise the
#' user value, both capped at `floor(n_eff / 3)` to keep at least three blocks.
#'
#' @param x Original input (for `stats::frequency`).
#' @param block_size `"auto"` or a positive integer.
#' @param n_eff Effective series length the engine will resample over.
#' @return Positive integer block size.
#' @keywords internal
#' @noRd
.resolve_block_size <- function(x, block_size, n_eff) {
  cap <- max(1L, floor(n_eff / 3L))
  if (identical(block_size, "auto")) {
    max(1L, min(2L * as.integer(stats::frequency(x)), cap))
  } else {
    max(1L, min(as.integer(block_size), cap))
  }
}


#' Stripped MBH trend for bootstrap iterations
#'
#' Runs mboost with a reduced iteration budget so that applying it hundreds
#' of times in `.boot_engine` is tractable. Uses the same spline geometry as
#' the base fit (same `knots`, `df`, `boundary.knots`) so bootstrap trends
#' live in the same function space.
#'
#' @param y      Numeric vector (synthetic bootstrap draw).
#' @param d_val  Calibrated Huber threshold from the base fit.
#' @param knots  Clamped knot count from the base fit.
#' @param mstop_boot Reduced iteration budget (typically `max(50, mstop %/% 5)`).
#' @param nu     Learning rate from the base fit.
#' @param df     P-spline degrees of freedom from the base fit.
#' @param boundary.knots Boundary knots from the base fit (or `NULL`).
#' @return Numeric trend vector, same length as `y`.
#' @keywords internal
#' @noRd
.mbh_fast_trend <- function(y, d_val, knots, mstop_boot, nu, df,
                             boundary.knots) {
  n        <- length(y)
  time_idx <- seq_len(n)

  bl_linear <- mboost::bols(time_idx, intercept = TRUE)
  bl_smooth <- mboost::bbs(time_idx, knots = knots, degree = 3L,
                            differences = 2L, df = df,
                            boundary.knots = boundary.knots)

  # Native silencing via trace = FALSE, avoids the per-iteration text
  # connection that capture.output() opens and tears down on every call.
  ctrl <- mboost::boost_control(mstop = mstop_boot, nu = nu, trace = FALSE)

  suppressWarnings({
    mod <- mboost::mboost(y ~ bl_linear + bl_smooth,
                          data = data.frame(y = y, time_idx = time_idx),
                          family = mboost::Huber(d = d_val),
                          control = ctrl)
  })
  as.numeric(stats::fitted(mod))
}


#' Universal block-bootstrap engine (Circular Block Bootstrap)
#'
#' Builds a long `data.table` of synthetic series (one per iteration) by
#' circular-block resampling of the pseudo-residuals, applies `filter_func`
#' to each block via `data.table`'s `by` grouping, and returns a 95%
#' normal-approximation band (`trend_base +/- 1.96 * bootstrap sd`) by time
#' point. The sd is used instead of empirical percentiles because it varies
#' smoothly and converges with far fewer iterations.
#'
#' Uses the Circular Block Bootstrap (Politis & Romano, 1992): block starts
#' are drawn uniformly from `1:n` and blocks wrap around the series end with
#' modular arithmetic. Every observation therefore appears in exactly
#' `block_size` candidate blocks, eliminating the end-point under-representation
#' of the plain moving-block bootstrap — critical because the trend's
#' end-of-sample uncertainty is precisely what we want to quantify.
#'
#' Index generation uses pure vector algebra:
#'   1. Sample `blk_n` block starts per iteration (all at once, from `1:n`).
#'   2. Expand each start with `rep()` + offset arithmetic.
#'   3. Wrap overshoot with `(idx - 1) %% n + 1` (circular, preserves variance).
#'   4. Flatten into the long `data.table` via positional arithmetic.
#'
#' @param filter_func Function `f(y_numeric) -> numeric trend`, same length.
#' @param trend_base  Numeric trend from the base fit (length `n`).
#' @param cycle_base  Numeric cycle (pseudo-residuals) from the base fit.
#' @param boot_iter   Positive integer number of iterations.
#' @param block_size  Block length (already resolved and validated by caller).
#' @return Named list with `lower` and `upper` numeric vectors (length `n`).
#' @keywords internal
#' @noRd
.boot_engine <- function(filter_func, trend_base, cycle_base,
                          boot_iter, block_size) {
  n       <- length(trend_base)
  blk_n   <- ceiling(n / block_size)
  raw_len <- blk_n * block_size            # >= n

  # Step 1-3: All indices for all iterations
  all_starts <- sample.int(n,
                            size    = boot_iter * blk_n,
                            replace = TRUE)

  offsets <- 0L:(block_size - 1L)
  raw_idx <- rep(all_starts, each = block_size) +
             rep(offsets, times  = boot_iter * blk_n)
  raw_idx <- (raw_idx - 1L) %% n + 1L      # circular wrap (no variance collapse)

  # Step 4: Map (iter, time) pairs → positions in raw_idx
  iter_id  <- rep(seq_len(boot_iter), each = n)
  time_pos <- rep(seq_len(n), times = boot_iter)
  flat_pos <- (iter_id - 1L) * raw_len + time_pos

  # Step 5: Build long data.table
  DT <- data.table::data.table(
    iter     = iter_id,
    time_idx = time_pos,
    y_boot   = trend_base[time_pos] + cycle_base[raw_idx[flat_pos]]
  )

  # Step 6: Apply filter; per time point, get the bootstrap standard deviation
  DT[, boot_trend := filter_func(y_boot), by = iter]

  s_DT <- DT[, list(s = stats::sd(boot_trend)), by = time_idx]
  data.table::setorder(s_DT, time_idx)

  # Step 7: Normal-approximation band centred on the point estimate.
  # Empirical 2.5/97.5 percentiles need hundreds of iterations to be smooth:
  # with a practical boot_iter the extreme order statistics are noisy and
  # leave visible kinks, worst at the high-variance end points. The bootstrap
  # standard deviation converges far faster and varies smoothly, so
  # trend_base +/- z * sd gives a clean band that is centred on the trend,
  # still widens honestly at the boundaries, and is free of jitter. Assumes
  # an approximately symmetric bootstrap law per time point (CLT for a smooth
  # functional), which holds well for a trend estimate.
  z <- stats::qnorm(0.975)
  list(
    lower = trend_base - z * s_DT$s,
    upper = trend_base + z * s_DT$s
  )
}
