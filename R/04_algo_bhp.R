# ── Boosted HP Filter (Phillips & Shi, 2021) ──────────────────────────────────
# Iteratively applies the HP filter on residuals to capture stochastic trends
# that a single HP pass misses.  Three stopping rules: BIC (default), ADF,
# and fixed iteration count.

#' Boosted HP Filter
#'
#' Iteratively applies the Hodrick-Prescott filter on residuals to better
#' capture stochastic trends.  At each iteration the HP smoother is applied to
#' the current residual and the resulting trend increment is added to the
#' cumulative trend estimate.  Iteration stops according to one of three rules:
#' BIC minimisation (default), ADF stationarity test on residuals, or a fixed
#' number of iterations.
#'
#' @param x Numeric vector, `ts`, `xts`, or `zoo` object.
#' @param lambda Smoothing parameter.
#'   If `NULL` (default), it is auto-detected using the Ravn-Uhlig rule
#'   (`6.25 * freq^4`).
#' @param iter_max Integer.
#'   Maximum number of boosting iterations (default 100).
#' @param stopping Character.
#'   Stopping rule: `"bic"` (default), `"adf"`, or `"fixed"`.
#' @param sig_level Numeric.
#'   Significance level for the ADF test when `stopping = "adf"`
#'   (default 0.05).
#' @param freq Numeric frequency override (1 = annual, 4 = quarterly,
#'   12 = monthly). Used only when `lambda` is `NULL` and the frequency cannot
#'   be inferred from `x`.
#'
#' @details
#' The boosted HP filter starts from the standard HP solution and then
#' re-applies the same HP smoother to the residual (cycle) component.  The
#' trend increment from each pass is accumulated, and the procedure stops
#' when one of the following criteria is met:
#'
#' \describe{
#'   \item{`"bic"`}{Schwarz information criterion computed as
#'     \eqn{n \log(\hat\sigma^2) + \log(n)\,\mathrm{tr}(S^m)},
#'     where \eqn{S^m} is the iterated smoother.  Iteration stops when the BIC
#'     increases relative to the previous best.}
#'   \item{`"adf"`}{Augmented Dickey-Fuller test on the residual.  Iteration
#'     stops when the residual is stationary at level `sig_level`.}
#'   \item{`"fixed"`}{Runs exactly `iter_max` iterations.}
#' }
#'
#' @return A `macrofilter` object with `trend`, `cycle`, `data`, and `meta`
#'   components.  The `meta` list contains `method = "bHP"`, `lambda`,
#'   `iterations`, `stopping_rule`, and `compute_time`.
#'
#' @references
#' Phillips, P.C.B. and Shi, Z. (2021). Boosting: Why You Can Use the HP
#'   Filter. *International Economic Review*, 62(2), 521--570.
#'
#' @export
#' @importFrom Matrix bandSparse Diagonal crossprod solve
#' @importFrom tseries adf.test
#'
#' @examples
#' # Quarterly GDP-like series
#' y <- ts(cumsum(rnorm(200)), start = c(2000, 1), frequency = 4)
#' result <- bhp_filter(y)
#' print(result)
bhp_filter <- function(x, lambda = NULL, iter_max = 100L,
                       stopping = c("bic", "adf", "fixed"),
                       sig_level = 0.05, freq = NULL) {

  # 1. Ingest ----------------------------------------------------------------
  inputs <- ensure_computable(x)
  y <- inputs$y
  n <- length(y)

  # 2. Length validation ------------------------------------------------------
  if (n < 3L) {
    stop("Series must have at least 3 observations for bHP filter.", call. = FALSE)
  }

  # 3. Match stopping rule ---------------------------------------------------
  stopping <- match.arg(stopping)
  iter_max <- as.integer(iter_max)

  # 4. Auto-lambda (Ravn-Uhlig heuristic) ------------------------------------
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

  # 5. Build sparse HP system once -------------------------------------------
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
  Q <- Matrix::Diagonal(n) + lambda * Matrix::crossprod(D)

  # 6. BIC precomputation (eigenvalues of Q) ---------------------------------
  if (stopping == "bic") {
    if (n > 5000L) {
      warning(
        "BIC stopping with n > 5000 requires O(n^3) eigendecomposition. ",
        "Consider `stopping = \"adf\"` or `stopping = \"fixed\"` for large series.",
        call. = FALSE
      )
    }
    eig <- eigen(as.matrix(Q), symmetric = TRUE, only.values = TRUE)$values
    # ratios[i] = (eig[i] - 1) / eig[i]  -- these are (1 - s_i)
    ratios <- (eig - 1) / eig
  }

  # 7. Start timer -----------------------------------------------------------
  t0 <- proc.time()

  # 8. Boosting loop ---------------------------------------------------------
  # Iteration 1
  res <- .hp_solve(y, Q)
  f_hat <- res$trend
  u     <- res$cycle

  final_iter  <- 1L
  final_trend <- f_hat
  final_cycle <- u

  if (stopping == "bic") {
    best_bic   <- n * log(sum(u^2) / n) + log(n) * (n - sum(ratios^1))
    best_iter  <- 1L
    best_trend <- f_hat
    best_cycle <- u
  }

  # Iterations 2..iter_max
  if (iter_max > 1L) {
    for (k in 2L:iter_max) {
      res_k <- .hp_solve(u, Q)
      f_hat <- f_hat + res_k$trend
      u     <- res_k$cycle

      if (stopping == "fixed") {
        final_iter  <- k
        final_trend <- f_hat
        final_cycle <- u
        next
      }

      if (stopping == "adf") {
        adf_p <- tseries::adf.test(u)$p.value
        if (adf_p < sig_level) {
          final_iter  <- k
          final_trend <- f_hat
          final_cycle <- u
          break
        }
        final_iter  <- k
        final_trend <- f_hat
        final_cycle <- u
        next
      }

      if (stopping == "bic") {
        tr_Sm <- n - sum(ratios^k)
        bic_k <- n * log(sum(u^2) / n) + log(n) * tr_Sm
        if (bic_k > best_bic) {
          # BIC increased -- stop and use previous best
          final_iter  <- best_iter
          final_trend <- best_trend
          final_cycle <- best_cycle
          break
        }
        best_bic   <- bic_k
        best_iter  <- k
        best_trend <- f_hat
        best_cycle <- u
        final_iter  <- k
        final_trend <- f_hat
        final_cycle <- u
      }
    }
  }

  elapsed <- (proc.time() - t0)[["elapsed"]]

  # 9. Package into macrofilter -----------------------------------------------
  result <- new_macrofilter(
    cycle = final_cycle,
    trend = final_trend,
    data  = y,
    meta  = list(
      method        = "bHP",
      lambda        = lambda,
      iterations    = final_iter,
      stopping_rule = stopping,
      compute_time  = elapsed
    )
  )

  # 10. Validate --------------------------------------------------------------
  validate_macrofilter(result)

  result
}

# ── Internal helper ─────────────────────────────────────────────────────────
#' Apply one HP smoothing pass given a precomputed system matrix
#'
#' @param y Numeric vector to smooth.
#' @param Q Precomputed sparse system matrix (I + lambda * D'D).
#' @return A list with `trend` and `cycle` numeric vectors.
#' @noRd
#' @keywords internal
.hp_solve <- function(y, Q) {
  trend <- as.numeric(Matrix::solve(Q, y))
  list(trend = trend, cycle = y - trend)
}
