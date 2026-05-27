# ── Data Ingestion: Adapter Pattern ───────────────────────────────────────────
# ensure_computable() normalizes any time-series input to a clean numeric vector
# plus metadata. restore_class() reverses the operation.

# ── Normalizer ───────────────────────────────────────────────────────────────

#' Normalize a time-series input to a clean numeric vector
#'
#' Detects the class of `x`, extracts values and temporal metadata, and
#' validates that the series contains no `NA` or `Inf` values.
#'
#' @param x A `numeric`, `ts`, `xts`, or `zoo` object.
#' @return A list with components:
#'   \describe{
#'     \item{y}{Clean `double` vector of values.}
#'     \item{idx}{Time index vector or `NULL` (for `numeric`/`ts`).}
#'     \item{tsp}{`tsp` attribute or `NULL` (only for `ts`).}
#'     \item{class}{Character: the original class name.}
#'     \item{xts_atts}{Named list of extra `xts` attributes, or `NULL`.}
#'   }
#' @keywords internal
#' @noRd
ensure_computable <- function(x) {
  cls <- class(x)[1L]

  handler <- switch(cls,
    "numeric" = .handle_numeric(x),
    "integer" = .handle_numeric(x),
    "ts"      = .handle_ts(x),
    "xts"     = .handle_xts(x),
    "zoo"     = .handle_zoo(x),
    # data.frame / data.table / anything else -> error
    stop(
      sprintf(
        "Unsupported input class '%s'. Pass a numeric vector, ts, xts, or zoo object.\n  Hint: for data.frames use df$column.",
        cls
      ),
      call. = FALSE
    )
  )

  .validate_values(handler$y)
  handler
}

# ── Class-specific handlers (internal) ───────────────────────────────────────

.handle_numeric <- function(x) {
  list(
    y        = as.double(x),
    idx      = NULL,
    tsp      = NULL,
    class    = "numeric",
    xts_atts = NULL
  )
}

.handle_ts <- function(x) {
  list(
    y        = as.double(x),
    idx      = NULL,
    tsp      = tsp(x),
    class    = "ts",
    xts_atts = NULL
  )
}

.handle_xts <- function(x) {
  if (!requireNamespace("zoo", quietly = TRUE)) {
    stop("Package 'zoo' is required to handle xts objects.", call. = FALSE)
  }
  # Preserve xtsAttributes (user-set metadata on the xts object)
  xa <- if (requireNamespace("xts", quietly = TRUE)) xts::xtsAttributes(x) else NULL

  list(
    y        = as.double(zoo::coredata(x)),
    idx      = zoo::index(x),
    tsp      = NULL,
    class    = "xts",
    xts_atts = xa
  )
}

.handle_zoo <- function(x) {
  if (!requireNamespace("zoo", quietly = TRUE)) {
    stop("Package 'zoo' is required to handle zoo objects.", call. = FALSE)
  }
  list(
    y        = as.double(zoo::coredata(x)),
    idx      = zoo::index(x),
    tsp      = NULL,
    class    = "zoo",
    xts_atts = NULL
  )
}

# ── Value validation ─────────────────────────────────────────────────────────

.validate_values <- function(y) {
  n_nan <- sum(is.nan(y))
  n_na  <- sum(is.na(y)) - n_nan   # is.na() counts NaN too
  n_inf <- sum(is.infinite(y))

  issues <- character(0L)
  if (n_na > 0L)  issues <- c(issues, sprintf("%d NA", n_na))
  if (n_nan > 0L) issues <- c(issues, sprintf("%d NaN", n_nan))
  if (n_inf > 0L) issues <- c(issues, sprintf("%d Inf/-Inf", n_inf))

  if (length(issues) > 0L) {
    stop(
      sprintf(
        "Input contains non-finite values (%s). Clean the series before filtering.",
        paste(issues, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

# ── Restorer ─────────────────────────────────────────────────────────────────

#' Restore a numeric result to the original time-series class
#'
#' @param y Numeric vector (trend or cycle result from an algorithm).
#' @param input_handler The list returned by `ensure_computable()`.
#' @return An object matching the original class of the input.
#' @importFrom stats ts tsp
#' @keywords internal
#' @noRd
restore_class <- function(y, input_handler) {
  switch(input_handler$class,
    "numeric" = as.double(y),
    "ts"      = .restore_ts(y, input_handler),
    "xts"     = .restore_xts(y, input_handler),
    "zoo"     = .restore_zoo(y, input_handler),
    as.double(y)  # fallback
  )
}

.restore_ts <- function(y, h) {
  ts(y, start = h$tsp[1L], frequency = h$tsp[3L])
}

.restore_xts <- function(y, h) {
  if (!requireNamespace("xts", quietly = TRUE)) {
    stop("Package 'xts' is required to restore xts objects.", call. = FALSE)
  }
  out <- xts::xts(y, order.by = h$idx)
  if (length(h$xts_atts) > 0L) {
    xts::xtsAttributes(out) <- h$xts_atts
  }
  out
}

.restore_zoo <- function(y, h) {
  if (!requireNamespace("zoo", quietly = TRUE)) {
    stop("Package 'zoo' is required to restore zoo objects.", call. = FALSE)
  }
  zoo::zoo(y, order.by = h$idx)
}