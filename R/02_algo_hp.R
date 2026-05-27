# ── HP Filter (Sparse Matrix Optimized) ──────────────────────────────────────
# Hodrick-Prescott filter using sparse Cholesky decomposition via the Matrix
# package.  Solves (I + lambda * D'D) t = y where D is the (n-2) x n second-
# difference matrix.

#' Hodrick-Prescott Filter (Sparse Matrix Implementation)
#'
#' Decomposes a time series into trend and cycle components by solving the
#' HP penalized least-squares problem using a sparse Cholesky factorization.
#' This avoids the dense O(n^3) inversion used by other implementations and
#' scales linearly in the number of observations.
#'
#' @param x Numeric vector, `ts`, `xts`, or `zoo` object.
#' @param lambda Smoothing parameter.
#'   If `NULL` (default), it is auto-detected using the Ravn-Uhlig rule
#'   (`6.25 * freq^4`).
#' @param freq Numeric frequency override (1 = annual, 4 = quarterly,
#'   12 = monthly). Used only when `lambda` is `NULL` and the frequency cannot
#'   be inferred from `x`.
#'
#' @details
#' The HP filter minimises
#' \deqn{\sum (y_t - \tau_t)^2 + \lambda \sum (\Delta^2 \tau_t)^2}
#' which admits the closed-form solution
#' \deqn{(I + \lambda D'D)\,\tau = y}
#' where \eqn{D} is the second-difference operator.
#'
#' The implementation builds \eqn{D} as a banded sparse matrix
#' ([Matrix::bandSparse()]) and solves the symmetric positive-definite system
#' with a sparse Cholesky decomposition ([Matrix::solve()]).
#'
#' When `lambda` is not supplied the Ravn-Uhlig (2002) rule is applied:
#' `lambda = 6.25 * freq^4`, yielding 6.25 (annual), 1600 (quarterly), and
#' 129 600 (monthly).
#'
#' @return A `macrofilter` object with `trend`, `cycle`, `data`, and `meta`
#'   components.
#'
#' @references
#' Hodrick, R.J. and Prescott, E.C. (1997). Postwar U.S. Business Cycles: An
#'   Empirical Investigation. *Journal of Money, Credit and Banking*, 29(1),
#'   1--16.
#'
#' Ravn, M.O. and Uhlig, H. (2002). On Adjusting the Hodrick-Prescott Filter
#'   for the Frequency of Observations. *Review of Economics and Statistics*,
#'   84(2), 371--376.
#'
#' @export
#' @importFrom Matrix bandSparse Diagonal crossprod solve
#'
#' @examples
#' # Quarterly GDP-like series
#' y <- ts(cumsum(rnorm(200)), start = c(2000, 1), frequency = 4)
#' result <- hp_filter(y)
#' print(result)
hp_filter <- function(x, lambda = NULL, freq = NULL) {

  # 1. Ingest ----------------------------------------------------------------
  inputs <- ensure_computable(x)
  y <- inputs$y
  n <- length(y)

  # 2. Length validation ------------------------------------------------------
  if (n < 3L) {
    stop("Series must have at least 3 observations for HP filter.", call. = FALSE)
  }

  # 3. Auto-lambda (Ravn-Uhlig heuristic) ------------------------------------
  if (is.null(lambda)) {
    if (!is.null(inputs$tsp)) {
      freq_detected <- inputs$tsp[3L]
    } else if (!is.null(freq)) {
      freq_detected <- freq
    } else {
      freq_detected <- 4
      warning(
        "Cannot determine series frequency; assuming quarterly (freq = 4). ",
        "Pass `lambda` or `freq` explicitly to silence this warning.",
        call. = FALSE
      )
    }
    lambda <- 6.25 * freq_detected^4
  }

  # 4-6. Sparse solve (timed) ------------------------------------------------
  t0 <- proc.time()

  # 4. Second-difference matrix D: (n-2) x n, bands at k = 0, 1, 2
  D <- Matrix::bandSparse(
    n    = n - 2L,
    m    = n,
    k    = c(0L, 1L, 2L),
    diagonals = list(
      rep(1,  n - 2L),
      rep(-2, n - 2L),
      rep(1,  n - 2L)
    )
  )

  # 5. Pentadiagonal system matrix Q = I + lambda * D'D
  Q <- Matrix::Diagonal(n) + lambda * Matrix::crossprod(D)

  # 6. Solve for trend
  trend_num <- as.numeric(Matrix::solve(Q, y))
  cycle_num <- y - trend_num

  elapsed <- (proc.time() - t0)[["elapsed"]]

  # 7. Package into macrofilter -----------------------------------------------
  result <- new_macrofilter(
    cycle = cycle_num,
    trend = trend_num,
    data  = y,
    meta  = list(
      method       = "HP",
      lambda       = lambda,
      compute_time = elapsed
    )
  )
  validate_macrofilter(result)

  result
}
