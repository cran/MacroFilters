# ── autoplot method for macrofilter objects ───────────────────────────────────

# Re-export the ggplot2 generic so users can call autoplot() on a macrofilter
# without attaching ggplot2 first.
#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot

#' Plot a macrofilter decomposition with ggplot2
#'
#' Visualises any `macrofilter` object: the observed series (attenuated) over
#' the estimated trend (emphasised). When the object carries bootstrap bands
#' (`$trend_lower` / `$trend_upper`), a 95% confidence ribbon is drawn behind
#' the lines automatically.
#'
#' @param object A `macrofilter` object.
#' @param ... Currently ignored; present for S3 generic compatibility.
#'
#' @return A `ggplot` object.
#'
#' @export
#' @method autoplot macrofilter
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon scale_color_manual
#' @importFrom ggplot2 labs theme_minimal theme element_blank
#' @importFrom data.table data.table
#'
#' @examples
#' \donttest{
#' y <- ts(cumsum(rnorm(120)), start = c(2000, 1), frequency = 4)
#' autoplot(mbh_filter(y, mstop = 100L, boot_iter = 50L))
#' }
autoplot.macrofilter <- function(object, ...) {

  # 4.1 Time axis + wide data ------------------------------------------------
  if (!inherits(object, "macrofilter")) {
    stop("`object` must be a 'macrofilter' object.", call. = FALSE)
  }

  if (!is.null(object$meta$tsp)) {
    tsp      <- object$meta$tsp
    time_vec <- seq(from = tsp[1L], to = tsp[2L], by = 1 / tsp[3L])
  } else {
    time_vec <- seq_along(object$trend)
  }

  DT <- data.table::data.table(
    time  = time_vec,
    Data  = as.numeric(object$data),
    Trend = as.numeric(object$trend)
  )

  has_ci <- !is.null(object$trend_lower) && !is.null(object$trend_upper)
  if (has_ci) {
    DT[, lower := as.numeric(object$trend_lower)]
    DT[, upper := as.numeric(object$trend_upper)]
  }

  # 4.2 Layered ggplot (ribbon at the back) ----------------------------------
  p <- ggplot2::ggplot(DT, ggplot2::aes(x = time))

  if (has_ci) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      fill = "#B0C4DE", alpha = 0.4, na.rm = TRUE
    )
  }

  p <- p +
    ggplot2::geom_line(
      ggplot2::aes(y = Data, color = "Observed Data"),
      alpha = 0.5, linewidth = 0.6, na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = Trend, color = "Structural Trend"),
      linewidth = 0.9, na.rm = TRUE
    )

  # 4.3 Aesthetics + metadata ------------------------------------------------
  m_name <- if (!is.null(object$meta$method)) object$meta$method else "Filter"

  p +
    ggplot2::scale_color_manual(
      values = c("Observed Data" = "gray40", "Structural Trend" = "#0047AB")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title    = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title    = paste(m_name, "Trend Extraction"),
      subtitle = if (has_ci) "With 95% Bootstrap Confidence Intervals" else NULL,
      x = "Time",
      y = "Value"
    )
}
