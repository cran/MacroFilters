## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 7,
  fig.height = 4.5,
  out.width  = "100%"
)

## ----setup, message=FALSE-----------------------------------------------------
library(MacroFilters)
library(ggplot2)
data("fr_gdp", package = "MacroFilters")

# France real GDP (log level) as a quarterly ts; includes the 2020 Q2 COVID shock
d0 <- fr_gdp$date[1]
fr <- ts(
  fr_gdp$gdp_log,
  start     = c(as.integer(format(d0, "%Y")),
                (as.integer(format(d0, "%m")) - 1L) %/% 3L + 1L),
  frequency = 4
)

## ----basic-band---------------------------------------------------------------
fit <- mbh_filter(fr, boot_iter = 50L)   # mstop defaults to 500
str(fit[c("trend_lower", "trend_upper")], max.level = 1)

## ----basic-plot---------------------------------------------------------------
autoplot(fit)

## ----block-size---------------------------------------------------------------
# Quarterly data -> auto block size = 2 * 4 = 8
fit_b <- mbh_filter(fr, boot_iter = 50L, block_size = 8L)

## ----all-filters, fig.height=3.6----------------------------------------------
autoplot(hp_filter(fr,       boot_iter = 50L))
autoplot(bhp_filter(fr,      boot_iter = 50L))
autoplot(hamilton_filter(fr, boot_iter = 50L))
autoplot(mbh_filter(fr,      boot_iter = 50L))

