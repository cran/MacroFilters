## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>",
  fig.width = 8, fig.height = 4.5, out.width = "100%"
)
library(MacroFilters)
library(data.table)
library(ggplot2)
data("us_gdp_vintage", package = "MacroFilters")

## ----nu-equiv-----------------------------------------------------------------
y <- us_gdp_vintage$gdp_log

res_default <- mbh_filter(y, mstop = 500,  nu = 0.10)
res_equiv   <- mbh_filter(y, mstop = 1000, nu = 0.05)

max_diff <- max(abs(res_default$trend - res_equiv$trend))
cat(sprintf("Max trend difference (mstop×nu equivalence): %.2e\n", max_diff))

## ----scale-invariance---------------------------------------------------------
y_level <- us_gdp_vintage$gdp_real   # billions USD (~20 000 scale)
y_log   <- us_gdp_vintage$gdp_log    # natural log (~10 scale)

d_level <- stats::mad(diff(y_level))
d_log   <- stats::mad(diff(y_log))

cat(sprintf("d (level series) : %.4f\n", d_level))
cat(sprintf("d (log series)   : %.6f\n", d_log))
cat(sprintf("Ratio d_level / mean(level): %.6f\n", d_level / mean(y_level)))
cat(sprintf("Ratio d_log   / mean(log)  : %.6f\n", d_log   / mean(y_log)))

## ----d-sensitivity------------------------------------------------------------
y_growth <- diff(us_gdp_vintage$gdp_log)   # quarterly log-differences

res_auto   <- mbh_filter(y_growth)
res_strict <- mbh_filter(y_growth, d = 0.005)
res_lenient <- mbh_filter(y_growth, d = 0.02)

cat(sprintf("Auto d = %.6f\n", res_auto$meta$d))

## ----d-sensitivity-plot-------------------------------------------------------
dt_growth <- data.table(
  t        = us_gdp_vintage$date[-1],
  observed = y_growth,
  auto     = res_auto$trend,
  strict   = res_strict$trend,
  lenient  = res_lenient$trend
)

dt_long <- melt(dt_growth,
                id.vars      = "t",
                measure.vars = c("auto", "strict", "lenient"),
                variable.name = "delta",
                value.name    = "trend")

# Human-readable labels
auto_label <- sprintf("Auto (d=%.4f)", res_auto$meta$d)
# data.table::melt() returns variable.name as factor; fcase() returns character.
# Assigning character to a factor column via := raises a type mismatch error,
# so coerce to character first.
dt_long[, delta := as.character(delta)]
dt_long[, delta := fcase(
  delta == "auto",    auto_label,
  delta == "strict",  "Strict (d=0.005)",
  delta == "lenient", "Lenient (d=0.020)"
)]

colour_vals <- c("#0072B2", "#009E73", "#E69F00")
names(colour_vals) <- c("Strict (d=0.005)", auto_label, "Lenient (d=0.020)")

p_d <- ggplot() +
  geom_line(
    data = dt_growth,
    aes(x = t, y = observed),
    colour = "grey70", linewidth = 0.5
  ) +
  geom_line(
    data = dt_long,
    aes(x = t, y = trend, colour = delta),
    linewidth = 0.9
  ) +
  annotate("rect",
           xmin = as.Date("2020-01-01"), xmax = as.Date("2020-10-01"),
           ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "firebrick") +
  annotate("text", x = as.Date("2020-04-01"), y = Inf,
           label = "COVID Q2\n-9% q-o-q", vjust = 1.4,
           size = 3.2, colour = "firebrick") +
  scale_colour_manual(values = colour_vals) +
  labs(
    title    = "MBH Trend Sensitivity to Huber Delta d",
    subtitle = "Data: US quarterly GDP growth rates (log-diff)",
    x        = NULL, y = "Log-difference", colour = "d setting"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

print(p_d)

## ----benchmark, cache=TRUE----------------------------------------------------
y          <- us_gdp_vintage$gdp_log
mstop_grid <- seq(100L, 1000L, by = 100L)   # 10 evenly-spaced points

bench_dt <- rbindlist(lapply(mstop_grid, function(m) {
  t0       <- proc.time()
  res      <- mbh_filter(y, mstop = m)
  elapsed  <- (proc.time() - t0)[["elapsed"]]
  cycle_sd <- sd(res$cycle)
  data.table(
    mstop       = m,
    elapsed_sec = round(elapsed, 3),
    cycle_sd    = round(cycle_sd, 6)
  )
}))

knitr::kable(
  bench_dt,
  col.names = c("mstop", "Wall time (s)", "Cycle SD"),
  caption   = "MBH computational benchmark — US log GDP (316 obs)"
)

## ----benchmark-plot-----------------------------------------------------------
# Dual-axis layout: wall time (left) + cycle_sd convergence (right)
# Use a secondary-axis trick by normalising cycle_sd to the time scale
time_range  <- range(bench_dt$elapsed_sec)
sd_range    <- range(bench_dt$cycle_sd)
# Guard against division by zero if cycle_sd converges to a flat line
if (diff(sd_range)   < 1e-10) sd_range   <- sd_range   + c(-1e-5,  1e-5)
if (diff(time_range) < 1e-10) time_range <- time_range + c(-1e-5,  1e-5)
sd_to_time  <- function(x) (x - sd_range[1]) / diff(sd_range) * diff(time_range) + time_range[1]
time_to_sd  <- function(x) (x - time_range[1]) / diff(time_range) * diff(sd_range) + sd_range[1]

p_bench <- ggplot(bench_dt, aes(x = mstop)) +
  geom_line(aes(y = elapsed_sec), colour = "#0072B2", linewidth = 1) +
  geom_point(aes(y = elapsed_sec), colour = "#CC0000", size = 3) +
  geom_line(aes(y = sd_to_time(cycle_sd)),
            colour = "#E69F00", linewidth = 0.9, linetype = "dashed") +
  geom_point(aes(y = sd_to_time(cycle_sd)),
             colour = "#E69F00", size = 2.5) +
  scale_x_continuous(breaks = mstop_grid) +
  scale_y_continuous(
    name     = "Wall time (s)  [blue / red points]",
    sec.axis = sec_axis(~ time_to_sd(.), name = "Cycle SD  [orange dashed]",
                        labels = scales::label_number(accuracy = 0.0001))
  ) +
  labs(
    title    = "Wall Time vs Boosting Iterations",
    subtitle = "US Real GDP log level (316 obs). Cycle SD plateaus well before mstop = 500.",
    x        = "mstop"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.y.left  = element_text(colour = "#0072B2"),
    axis.title.y.right = element_text(colour = "#E69F00")
  )

print(p_bench)

## ----summary-table, echo=FALSE------------------------------------------------
summary_tbl <- data.table(
  Parameter = c("`mstop`", "`nu`", "`knots`", "`d`"),
  Default   = c("500", "0.1", "`max(20, n/2)`", "auto via MAD"),
  `When to increase` = c(
    "Publication accuracy required",
    "Very long series; computational budget tight",
    "Highly nonlinear trend",
    "Series has frequent large spikes"
  ),
  `When to decrease` = c(
    "Exploratory / fast iteration",
    "Stability preferred over speed",
    "Short series or near-linear trend",
    "Series is log-level (use `mad(hp$cycle)` instead)"
  )
)
knitr::kable(summary_tbl, caption = "MBH hyperparameter quick-reference")

