# ── Hamilton Filter (Vectorized OLS) ──────────────────────────────────────────
# James Hamilton (2018) regression-based alternative to the HP filter.
# Regresses y_{t+h} on (1, y_t, y_{t-1}, ..., y_{t-p+1}) via OLS and defines
# the cycle as the residual.

#' Hamilton Filter
#'
#' Decomposes a time series into trend and cycle components using the
#' regression-based filter proposed by Hamilton (2018). The trend is the
#' fitted value from an OLS regression of \eqn{y_{t+h}} on
#' \eqn{(1, y_t, y_{t-1}, \ldots, y_{t-p+1})}, and the cycle is the residual.
#'
#' @param x Numeric vector, `ts`, `xts`, or `zoo` object.
#' @param h Integer horizon (number of periods ahead). If `NULL` (default),
#'   auto-detected from the series frequency using Hamilton's rule:
#'   annual = 2, quarterly = 8, monthly = 24.
#' @param p Integer number of lags in the regression (default 4).
#'
#' @details
#' Hamilton (2018) proposes replacing the HP filter with a simple regression:
#' \deqn{y_{t+h} = \beta_0 + \beta_1 y_t + \beta_2 y_{t-1} + \cdots +
#'   \beta_p y_{t-p+1} + v_{t+h}}
#' The fitted values \eqn{\hat{y}_{t+h}} define the trend and the residuals
#' \eqn{\hat{v}_{t+h}} define the cycle.
#'
#' The first \eqn{h + p - 1} observations have no computable trend or cycle
#' and are filled with `NA`.
#'
#' The lag matrix is constructed vectorized via `embed()` and the
#' regression is solved with [stats::lm.fit()] for speed.
#'
#' @return A `macrofilter` object with `trend`, `cycle`, `data`, and `meta`
#'   components. `meta` includes `h`, `p`, `coefficients`, and `compute_time`.
#'
#' @references
#' Hamilton, J.D. (2018). Why You Should Never Use the Hodrick-Prescott
#'   Filter. *Review of Economics and Statistics*, 100(5), 831--843.
#'
#' @export
#' @importFrom stats embed lm.fit
#'
#' @examples
#' # Quarterly GDP-like series
#' y <- ts(cumsum(rnorm(200)), start = c(2000, 1), frequency = 4)
#' result <- hamilton_filter(y)
#' print(result)
hamilton_filter <- function(x, h = NULL, p = 4L) {

  # 1. Ingest ----------------------------------------------------------------
  inputs <- ensure_computable(x)
  y <- inputs$y
  n <- length(y)

  # 2. Auto-h (horizon) ------------------------------------------------------
  if (is.null(h)) {
    if (!is.null(inputs$tsp)) {
      freq_detected <- inputs$tsp[3L]
    } else {
      freq_detected <- 4
      warning(
        "Cannot determine series frequency; assuming quarterly (freq = 4). ",
        "Pass `h` explicitly to silence this warning.",
        call. = FALSE
      )
    }
    h <- switch(as.character(freq_detected),
      "1"  = 2L,
      "4"  = 8L,
      "12" = 24L,
      as.integer(2 * freq_detected)
    )
  }
  h <- as.integer(h)
  p <- as.integer(p)

  # 3. Length validation ------------------------------------------------------
  min_obs <- h + p
  if (n <= min_obs) {
    stop(
      sprintf(
        "Series too short for selected horizon/lags (n = %d, need > %d = h + p).",
        n, min_obs
      ),
      call. = FALSE
    )
  }

  # 4. Construct regression matrix (vectorized) -------------------------------
  t0 <- proc.time()

  # embed(y, p) gives an (n - p + 1) x p matrix.
  # Row i corresponds to time t = p + i - 1 (1-based).
  # Column 1 = y_t, column 2 = y_{t-1}, ..., column p = y_{t-p+1}.
  E <- embed(y, p)

  # Keep only rows where the target y[t+h] exists: t + h <= n
  # Row i maps to t = p + i - 1, so we need i <= n - h - p + 1.
  n_valid <- n - h - p + 1L
  X <- cbind(1, E[seq_len(n_valid), , drop = FALSE])

  # Target: y[t+h] for t = p, p+1, ..., n-h
  target <- y[(p + h):n]

  # 5. OLS via lm.fit --------------------------------------------------------
  fit <- lm.fit(X, target)

  # 6. Reconstruct with NA padding --------------------------------------------
  trend_num <- rep(NA_real_, n)
  cycle_num <- rep(NA_real_, n)

  # Fitted values / residuals correspond to times (p + h) through n
  fitted_idx <- (p + h):n
  trend_num[fitted_idx] <- fit$fitted.values
  cycle_num[fitted_idx] <- fit$residuals

  elapsed <- (proc.time() - t0)[["elapsed"]]

  # 7. Package into macrofilter -----------------------------------------------
  result <- new_macrofilter(
    cycle = cycle_num,
    trend = trend_num,
    data  = y,
    meta  = list(
      method       = "Hamilton",
      h            = h,
      p            = p,
      coefficients = fit$coefficients,
      compute_time = elapsed
    )
  )
  validate_macrofilter(result)

  result
}
