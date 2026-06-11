# ── MacroBoost Hybrid Filter (MBH) ────────────────────────────────────────────
# The flagship algorithm. Combines P-Splines with Component-wise Gradient
# Boosting and Huber loss to achieve robustness against structural breaks.

#' MacroBoost Hybrid (MBH) Filter
#'
#' Decomposes a time series into trend and cycle using a robust boosting
#' algorithm. Unlike the HP filter, MBH uses the Huber loss function to
#' automatically downweight outliers (like the COVID-19 shock), preventing
#' them from distorting the trend.
#'
#' @param x Numeric vector, `ts`, `xts`, or `zoo` object.
#' @param knots Integer.
#'   Number of interior knots for the P-Spline.
#'   If `NULL` (default), it is calculated as `min(max(20, floor(n / 2)), 250)`.
#'   High knot density keeps the trend flexible, while the cap of 250 keeps the
#'   B-spline basis bounded for long / high-frequency series: in a P-spline the
#'   smoothness is governed by the difference penalty (via `df`, `mstop`, `nu`),
#'   not by the knot count, so beyond a few hundred knots the extra basis
#'   columns only inflate memory and runtime without adding useful flexibility.
#' @param mstop Integer.
#'   Maximum number of boosting iterations (default 500).
#'   If `select_mstop = TRUE` this is the upper bound; the actual stopping
#'   point is chosen by AICc.
#'
#'   **Under-smoothing warning (`mstop` vs `d`):** when `d` is small relative
#'   to the trend's range -- the typical case for long log-level series, where
#'   the cycle (hence the auto-calibrated `d`) is tiny but the trend spans a
#'   large range -- the Huber loss caps the gradient from the first iteration,
#'   so each boosting step advances the trend only slightly. Reducing `mstop`
#'   then leaves the trend unable to climb its full range: it collapses to a
#'   nearly flat curve while the cycle absorbs the long-run variation. Keep the
#'   default `mstop = 500` (or higher) for such series; lower it only for short
#'   or high-cycle-variance inputs.
#' @param d Numeric or `"auto"`.
#'   The delta parameter for Huber loss. If `"auto"` (default), it is
#'   calibrated as `stats::mad(.hp_fast(x))`, i.e. the MAD of the HP cyclical
#'   residual. This anchors the threshold to the output-gap scale rather than
#'   the growth-rate scale, avoiding the under-truncation failure mode of the
#'   legacy `mad(diff(y))` heuristic. A `message()` is emitted reporting the
#'   exact value chosen. Supply an explicit positive numeric to override.
#' @param boot_iter Non-negative integer.
#'   Number of block-bootstrap iterations for uncertainty quantification
#'   (default `0`, bootstrap disabled). When `> 0`, the function adds
#'   `$trend_lower` and `$trend_upper`: a 95% normal-approximation band,
#'   `trend +/- 1.96 * sd(bootstrap trends)`, centred on the estimated trend.
#'   The bootstrap sd is used instead of empirical percentiles because it is
#'   smooth and stable at practical `boot_iter`. Each bootstrap refit uses the
#'   same `mstop` as the base fit, so larger `boot_iter` raises cost linearly.
#'   See also `block_size`.
#' @param block_size Positive integer or `"auto"`.
#'   Block length for the moving-block bootstrap (used only when
#'   `boot_iter > 0`). If `"auto"` (default), it is set to
#'   `2 * stats::frequency(x)` (two full cycles), bounded above by
#'   `floor(length(x) / 3)` to keep at least three blocks.
#' @param hp_lambda Numeric or `NULL`.
#'   Smoothing parameter for the internal HP filter used to auto-calibrate
#'   `d` (only relevant when `d = "auto"`). If `NULL` (default), it is derived
#'   from `stats::frequency(x)` via the Ravn-Uhlig rule. **Supply this when
#'   `x` is a plain numeric vector whose true frequency is not annual**, since
#'   `frequency()` returns `1` for unclassed vectors and would otherwise
#'   under-smooth the calibration cycle (e.g. monthly data: `hp_lambda = 129600`).
#' @param nu Numeric.
#'   The learning rate (shrinkage) for boosting (default 0.1).
#' @param df Integer.
#'   Effective degrees of freedom per boosting step for the P-Spline base
#'   learner (default 4). This enforces the *weak-learner* constraint of
#'   Bühlmann & Hothorn (2007): each boosting step contributes only a small,
#'   smooth update so that the trend is built up gradually over many
#'   iterations rather than fitted in one pass.
#'
#'   **End-point instability warning:** Higher `df` values cause the B-spline
#'   basis matrix to shift drastically when the sample size changes by even
#'   one observation (the "rubber-band effect"). The last few data points pull
#'   the estimated trend non-smoothly, producing unreliable end-of-sample
#'   estimates. Keep `df = 4` (the default) unless you have a specific reason
#'   to deviate.
#' @param select_mstop Logical.
#'   If `TRUE`, the optimal number of boosting iterations is selected
#'   automatically via AICc (corrected AIC), following Bühlmann & Hothorn
#'   (2007). The `mstop` argument acts as the search upper bound. Default
#'   `FALSE`.
#'
#'   **AICc underfitting warning:** In the combination of Huber
#'   quasi-likelihood + P-splines, AICc penalises model complexity
#'   hyper-aggressively. In practice the algorithm stops at iteration ~5–15
#'   instead of the intended ~500. The resulting trend is nearly a straight
#'   line; all long-run variance is pushed into the cycle component, defeating
#'   the purpose of the filter. Treat `select_mstop = TRUE` as an
#'   experimental option and validate visually before relying on it.
#' @param boundary.knots A numeric vector of length 2 specifying the global
#'   domain for the B-spline basis (e.g., `c(1, T_max)`). If `NULL` (default),
#'   the range of `time_idx` is used. For real-time stability, fix this to the
#'   full-sample domain so the basis does not shift as the sample grows.
#'
#' @details
#' The model estimated is an additive model:
#' \deqn{y_t = \text{Linear}(t) + \text{Smooth}(t) + \epsilon_t}
#'
#' It is fitted using [mboost::mboost()] with:
#' \itemize{
#'   \item **Base Learners:** A linear time trend ([mboost::bols()]) to capture
#'     the global path, plus a B-spline ([mboost::bbs()]) to capture local
#'     curvature.
#'   \item **Loss Function:** Huber loss ([mboost::Huber()]) with parameter
#'     `d`.  This is the key to robustness.
#' }
#'
#' The default parameters (`knots = min(n/2, 250)`, `mstop = 500`) are
#' calibrated to mimic the flexibility of a standard HP filter while retaining
#' the robustness of the Huber loss.
#'
#' @section Calibration Guidance:
#' Three failure modes were discovered through empirical stress-testing.
#' The defaults guard against all three:
#'
#' \describe{
#'   \item{1. Huber delta scale mismatch (`d`)}{
#'     The automatic fallback `mad(diff(y))` operates on the scale of
#'     growth rates, not the output gap. For log-level input this sets `d`
#'     one to two orders of magnitude too small, causing ordinary
#'     business-cycle swings to be treated as outliers. If the estimated
#'     cycle looks implausibly large or the trend is nearly linear, override
#'     with `d = mad(hp_filter(x)$cycle)` as a starting point.
#'   }
#'   \item{2. AICc underfitting (`select_mstop`)}{
#'     AICc + Huber quasi-likelihood + P-splines stops boosting at
#'     iteration ~5--15. The trend degenerates to a near-straight line and
#'     the cycle absorbs all long-run variance. Leave `select_mstop = FALSE`
#'     (the default) and set `mstop` explicitly instead.
#'   }
#'   \item{3. End-point instability (`df`)}{
#'     Values above 4 shift the B-spline basis matrix non-smoothly as
#'     the sample grows, producing a "rubber-band" distortion in the final
#'     observations. Keep `df = 4` (the default) for real-time applications.
#'   }
#' }
#'
#' @return A list of class `c("macrofilter", "list")` with:
#'   \describe{
#'     \item{`$trend`}{Numeric trend vector.}
#'     \item{`$cycle`}{Numeric cycle vector.}
#'     \item{`$data`}{Original input as numeric.}
#'     \item{`$meta`}{Named list: `method`, `knots`, `d`, `mstop`, `nu`,
#'       `df`, `select_mstop`, `compute_time`.}
#'     \item{`$trend_lower`, `$trend_upper`}{95% normal-approximation bootstrap
#'       band (`trend +/- 1.96 * sd`). Present only when `boot_iter > 0`.}
#'   }
#'
#' @export
#' @importFrom mboost mboost bbs bols Huber boost_control mstop
#' @importFrom stats fitted AIC frequency mad sd qnorm
#' @importFrom data.table data.table setorder
#'
#' @examples
#' # Fast example with reduced series and iterations
#' set.seed(42)
#' y <- ts(cumsum(rnorm(80)), start = c(2000, 1), frequency = 4)
#' result <- mbh_filter(y, mstop = 100L)
#' print(result)
#'
#' \donttest{
#' # Full example with default parameters
#' y2 <- ts(cumsum(rnorm(200)), start = c(2000, 1), frequency = 4)
#' result2 <- mbh_filter(y2)
#' print(result2)
#' }
mbh_filter <- function(x, d = "auto", boot_iter = 0, block_size = "auto",
                       knots = NULL, mstop = 500L, nu = 0.1, df = 4L,
                       select_mstop = FALSE, boundary.knots = NULL,
                       hp_lambda = NULL) {

  # 1. Ingest ----------------------------------------------------------------
  inputs <- ensure_computable(x)
  y <- inputs$y
  n <- length(y)

  # 2. Length validation ------------------------------------------------------
  if (n < 6L) {
    stop("Series must have at least 6 observations for MBH filter.", call. = FALSE)
  }

  # 3. Auto-calibrate d (Huber delta) ----------------------------------------
  # MAD of the HP cycle anchors the threshold to the output-gap scale,
  # avoiding the scale-mismatch failure of the legacy mad(diff(y)) heuristic.
  # hp_lambda guards calibration for frequency-less numeric input (see @param).
  if (identical(d, "auto")) {
    d_val <- stats::mad(.hp_fast(x, lambda = hp_lambda))
    if (d_val < 1e-6) d_val <- 0.01
    message(sprintf(
      "Info: Huber threshold automatically calibrated to d = %.6f via HP cyclical MAD.",
      d_val
    ))
  } else {
    d_val <- as.double(d)
    if (d_val < 1e-6) d_val <- 0.01
  }

  # 4. Knots heuristic (aggressive) ------------------------------------------
  # High flexibility (~1 knot every 2 obs) because the Huber loss will be the
  # smoothness constraint, not the spline basis.
  if (is.null(knots)) {
    knots <- min(max(20L, floor(n / 2L)), 250L)
  }
  # Clamp knots to avoid singular spline bases with very short series
  max_knots <- max(1L, n - 4L)
  knots <- min(as.integer(knots), max_knots)
  mstop <- as.integer(mstop)
  df    <- as.integer(df)

  # 5. Timer -----------------------------------------------------------------
  t0 <- proc.time()

  # 6. Base learners ---------------------------------------------------------
  time_idx  <- seq_len(n)
  df_boost  <- data.frame(y = y, time_idx = time_idx)
  bl_linear <- mboost::bols(time_idx, intercept = TRUE)
  bl_smooth <- mboost::bbs(time_idx, knots = knots, degree = 3,
                            differences = 2, df = df,
                            boundary.knots = boundary.knots)

  # 7. Fit -------------------------------------------------------------------
  fam  <- mboost::Huber(d = d_val)
  # Native silencing via trace = FALSE instead of capturing stdout.
  ctrl <- mboost::boost_control(mstop = mstop, nu = nu, trace = FALSE)

  suppressWarnings(
    mod <- mboost::mboost(
      y ~ bl_linear + bl_smooth,
      data = df_boost, family = fam, control = ctrl
    )
  )

  # 7.5. Optional AICc-based mstop selection ---------------------------------
  mstop_final <- mstop
  if (select_mstop) {
    aic_obj <- tryCatch(
      stats::AIC(mod, method = "corrected"),
      error = function(e) NULL
    )
    if (!is.null(aic_obj)) {
      best_stop   <- mboost::mstop(aic_obj)
      mod         <- mod[best_stop]
      mstop_final <- as.integer(best_stop)
    }
  }

  # 8. Extract ---------------------------------------------------------------
  trend_num <- as.numeric(fitted(mod))
  cycle_num <- y - trend_num

  elapsed <- (proc.time() - t0)[["elapsed"]]

  # 9. Optional block bootstrap -----------------------------------------------
  ci_bands <- NULL
  if (boot_iter > 0L) {
    bs <- .resolve_block_size(x, block_size, n)
    # Bootstrap refits MUST use the same estimator (identical mstop) as the
    # base fit, or the band width is biased. boot_iter is the speed dial.
    mstop_boot <- mstop_final
    ff <- function(y_b) {
      .mbh_fast_trend(y_b, d_val, knots, mstop_boot, nu, df, boundary.knots)
    }
    ci_bands <- .boot_engine(
      filter_func = ff,
      trend_base  = trend_num,
      cycle_base  = cycle_num,
      boot_iter   = as.integer(boot_iter),
      block_size  = bs
    )
  }

  # 10. Build S3 result -------------------------------------------------------
  result <- list(
    trend = trend_num,
    cycle = cycle_num,
    data  = y,
    meta  = list(
      method         = "MBH",
      knots          = knots,
      d              = d_val,
      mstop          = mstop_final,
      mstop_initial  = mstop,
      nu             = nu,
      df             = df,
      select_mstop   = select_mstop,
      boundary.knots = boundary.knots,
      compute_time   = elapsed,
      # Temporal identity: stored once so trend, cycle AND bands can all be
      # mapped back to real dates downstream (e.g. autoplot x-axis).
      ts_class       = inputs$class,
      tsp            = inputs$tsp,
      idx            = inputs$idx
    )
  )
  if (!is.null(ci_bands)) {
    result$trend_lower <- ci_bands$lower
    result$trend_upper <- ci_bands$upper
  }
  class(result) <- c("macrofilter", "list")

  # 11. Validate --------------------------------------------------------------
  validate_macrofilter(result)

  result
}
