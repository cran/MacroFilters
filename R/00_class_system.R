# ── MacroFilters S3 Class System ──────────────────────────────────────────────
# Low-level constructor, validator, and print method for the `macrofilter` class.

# ── Constructor (low-level) ──────────────────────────────────────────────────

#' Create a macrofilter object (low-level)
#'
#' @param cycle Numeric vector. The cyclical component.
#' @param trend Numeric vector. The trend component.
#' @param data  Numeric vector. The original series (immutable).
#' @param meta  A named list with at least `method` (character) and optionally
#'   `lambda`, `lags`, `d`, `iter`, `compute_time`.
#' @return An object of class `macrofilter`.
#' @keywords internal
#' @noRd
new_macrofilter <- function(cycle, trend, data, meta) {
  stopifnot(is.numeric(cycle))
  stopifnot(is.numeric(trend))
  stopifnot(is.numeric(data))
  stopifnot(is.list(meta), !is.null(meta$method))

  structure(
    list(
      cycle = cycle,
      trend = trend,
      data  = data,
      meta  = meta
    ),
    class = "macrofilter"
  )
}

# ── Validator ────────────────────────────────────────────────────────────────

#' Validate a macrofilter object
#'
#' Checks that `trend + cycle` reconstructs the original data within floating
#' point tolerance.
#'
#' @param x A `macrofilter` object.
#' @param tol Numeric tolerance for the identity check (default `1e-9`).
#' @return `x`, invisibly, if valid.
#' @keywords internal
#' @noRd
validate_macrofilter <- function(x, tol = 1e-9) {
  n_cycle <- length(x$cycle)
  n_trend <- length(x$trend)
  n_data  <- length(x$data)

  if (n_cycle != n_data || n_trend != n_data) {
    stop(
      sprintf(
        "Length mismatch: data has %d obs, trend has %d, cycle has %d.",
        n_data, n_trend, n_cycle
      ),
      call. = FALSE
    )
  }

  residual <- abs(x$trend + x$cycle - x$data)
  # Ignore positions where any component is NA (e.g. Hamilton filter pads with NAs)
  check_idx <- !is.na(x$data) & !is.na(x$trend) & !is.na(x$cycle)

  if (!any(check_idx)) {
    stop("No complete cases to validate (all positions contain NA).", call. = FALSE)
  }

  max_err <- max(residual[check_idx])

  if (max_err > tol) {
    stop(
      sprintf(
        "Identity violation: max|trend + cycle - data| = %.2e (tol = %.2e).",
        max_err, tol
      ),
      call. = FALSE
    )
  }

  invisible(x)
}

# ── Print method ─────────────────────────────────────────────────────────────

#' @export
#' @importFrom stats sd
print.macrofilter <- function(x, ...) {
  m <- x$meta
  n <- length(x$data)

  cat(sprintf("-- MacroFilter [%s] --\n", m$method))
  cat(sprintf("   Observations : %d\n", n))


  # Show relevant scalar parameters per method
  params <- setdiff(names(m), c("method", "compute_time"))
  params <- Filter(function(p) length(m[[p]]) == 1L, params)
  if (length(params) > 0L) {
    vals <- vapply(params, function(p) format(m[[p]], digits = 4), character(1))
    cat(sprintf("   Parameters   : %s\n", paste(params, vals, sep = " = ", collapse = ", ")))
  }

  # Cycle summary (skip NAs)
  cy <- x$cycle[!is.na(x$cycle)]
  if (length(cy) > 0L) {
    cat(sprintf(
      "   Cycle range  : [%s, %s]  sd = %s\n",
      format(min(cy), digits = 4),
      format(max(cy), digits = 4),
      format(sd(cy),  digits = 4)
    ))
  }

  if (!is.null(m$compute_time)) {
    cat(sprintf("   Compute time : %.3f s\n", as.numeric(m$compute_time)))
  }

  invisible(x)
}