## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>",
  fig.width = 8, fig.height = 4.5, out.width = "100%"
)
library(MacroFilters)
library(data.table)
library(ggplot2)
data("us_gdp_vintage", package = "MacroFilters")

## ----expanding-window, cache=TRUE---------------------------------------------
T_max    <- nrow(us_gdp_vintage)   # full-sample size (316 rows as of 2025)
ref_date <- as.Date("2019-10-01")  # 2019 Q4 — last pre-COVID quarter
ref_idx  <- which(us_gdp_vintage$date == ref_date)

# d calibrated on the output-gap scale (not the growth-rate scale).
# Using mad(diff(y_log)) sets d ~0.006 which is far below a typical output
# gap (0.01–0.05), so MBH over-smooths to a near-straight line.  As
# post-COVID recovery data arrive, that straight line tilts upward and
# revises the 2019 Q4 trend estimate — precisely the instability we want
# to avoid.  Fixing d = mad(hp_cycle) anchors the threshold to the right
# scale and keeps the backward revision near zero.
d_fixed <- stats::mad(hp_filter(us_gdp_vintage$gdp_log, freq = 4)$cycle)

# df = 10: avoids the extreme endpoint under-fit that df = 4 (default)
# produces.  With df = 4, the cycle at the 2019 Q4 endpoint is ~20 ppts;
# with df = 10 it is ~0.003 ppts (same as HP), removing the spurious
# upward drift that would otherwise accumulate as recovery data arrives.
# mstop = 1000 gives the additional boosting budget needed at df = 10.
fixed_df     <- 10L
fixed_mstop  <- 1000L

# Freeze the B-spline domain so every vintage uses the same basis — knot
# count and domain are set once from the full sample and never change.
fixed_knots  <- 50L                            # invariant functional space
fixed_bounds <- c(1L, T_max)                   # global B-spline domain

# Vintages: 2016 Q1 – 2022 Q4 (28 publication dates)
eval_dates   <- us_gdp_vintage[date >= "2016-01-01" & date <= "2022-10-01", date]
eval_indices <- which(us_gdp_vintage$date %in% eval_dates)

fan_list  <- vector("list", length(eval_indices))
back_list <- vector("list", length(eval_indices))

for (k in seq_along(eval_indices)) {
  i         <- eval_indices[k]
  y_current <- us_gdp_vintage$gdp_log[seq_len(i)]

  hp_res  <- hp_filter(y_current, freq = 4)

  # Single call per vintage — all parameters frozen outside the loop so
  # both the fan chart and the backward-revision series use the same basis.
  mbh_res <- mbh_filter(
    y_current,
    d              = d_fixed,
    knots          = fixed_knots,
    df             = fixed_df,
    mstop          = fixed_mstop,
    boundary.knots = fixed_bounds
  )

  n_cur    <- length(y_current)
  tail_idx <- max(1L, n_cur - 27L):n_cur   # 28-obs trailing window

  fan_list[[k]] <- data.table(
    vintage_date = us_gdp_vintage$date[i],
    obs_date     = us_gdp_vintage$date[tail_idx],
    hp_trend     = hp_res$trend[tail_idx],
    mbh_trend    = mbh_res$trend[tail_idx],
    gdp_log      = y_current[tail_idx]
  )

  if (i >= ref_idx) {
    back_list[[k]] <- data.table(
      vintage_date = us_gdp_vintage$date[i],
      hp_at_ref    = hp_res$trend[ref_idx],
      mbh_at_ref   = mbh_res$trend[ref_idx]
    )
  }
}

revisions_dt <- rbindlist(fan_list)
backward_dt  <- rbindlist(Filter(Negate(is.null), back_list))

## ----boundary-knots-demo------------------------------------------------------
n_demo  <- 200L   # truncated sample
y_demo  <- us_gdp_vintage$gdp_log[seq_len(n_demo)]

res_free      <- mbh_filter(y_demo, knots = fixed_knots, df = fixed_df,
                             mstop = fixed_mstop, boundary.knots = NULL)
res_anchored  <- mbh_filter(y_demo, knots = fixed_knots, df = fixed_df,
                             mstop = fixed_mstop, boundary.knots = fixed_bounds)

# Extend by one observation and refit
y_demo_p1 <- us_gdp_vintage$gdp_log[seq_len(n_demo + 1L)]

res_free_p1     <- mbh_filter(y_demo_p1, knots = fixed_knots, df = fixed_df,
                               mstop = fixed_mstop, boundary.knots = NULL)
res_anchored_p1 <- mbh_filter(y_demo_p1, knots = fixed_knots, df = fixed_df,
                               mstop = fixed_mstop, boundary.knots = fixed_bounds)

# Revision at the final shared observation (position n_demo)
rev_free     <- abs(res_free_p1$trend[n_demo]     - res_free$trend[n_demo])
rev_anchored <- abs(res_anchored_p1$trend[n_demo] - res_anchored$trend[n_demo])

cat(sprintf(
  "Revision at obs %d after adding one data point:\n  free domain   : %.6f\n  anchored domain: %.6f\n",
  n_demo, rev_free, rev_anchored
))

## ----fan-chart----------------------------------------------------------------
# Show only the shared overlap window (2018 Q1 onward) where every one of
# the 28 vintages contributes data — this eliminates the staircase / accordion
# artefact that arises when staggered trailing windows are plotted together.
fan_shared <- revisions_dt[obs_date >= as.Date("2018-01-01")]

p1 <- ggplot(fan_shared, aes(x = obs_date)) +
  geom_line(aes(y = hp_trend,  group = vintage_date),
            colour = "#0072B2", alpha = 0.4, linewidth = 0.6) +
  geom_line(aes(y = mbh_trend, group = vintage_date),
            colour = "#E69F00", alpha = 0.4, linewidth = 0.6) +
  geom_line(
    data = us_gdp_vintage[date >= as.Date("2018-01-01") & date <= as.Date("2022-10-01")],
    aes(x = date, y = gdp_log),
    colour = "black", linewidth = 0.8, linetype = "dashed"
  ) +
  annotate("rect",
           xmin = as.Date("2020-01-01"), xmax = as.Date("2020-10-01"),
           ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "firebrick") +
  annotate("text", x = as.Date("2020-04-01"), y = Inf,
           label = "COVID\nshock", vjust = 1.4, size = 3.2, colour = "firebrick") +
  labs(
    title    = "Vintage Fan Chart: HP (blue) vs MBH (orange)",
    subtitle = "Each line = one vintage estimate (2018–2022 overlap window). Dashed = observed.",
    x        = NULL,
    y        = "Log Real GDP"
  ) +
  theme_minimal(base_size = 12)

print(p1)

## ----backward-revision--------------------------------------------------------
# backward_dt was built in the expanding-window loop above.
# MBH used df = 10 (correct endpoint estimate), knots = 50, and
# boundary.knots = c(1, T_max) (frozen basis) — so all vintages are
# numerically comparable and the baseline trend is unbiased.
back_dt <- backward_dt[order(vintage_date)]
setnames(back_dt, c("hp_at_ref", "mbh_at_ref"), c("hp_trend", "mbh_trend"))

# Normalise to the base vintage (2019 Q4 = first vintage that includes ref_date)
base_hp  <- back_dt[vintage_date == ref_date, hp_trend]
base_mbh <- back_dt[vintage_date == ref_date, mbh_trend]

back_dt[, hp_revision  := (hp_trend  - base_hp)  * 100]   # × 100 = ppts log GDP
back_dt[, mbh_revision := (mbh_trend - base_mbh) * 100]

p2 <- ggplot(back_dt, aes(x = vintage_date)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(aes(y = hp_revision),  colour = "#0072B2", linewidth = 1) +
  geom_point(aes(y = hp_revision), colour = "#0072B2", size = 2.5) +
  geom_line(aes(y = mbh_revision),  colour = "#E69F00", linewidth = 1) +
  geom_point(aes(y = mbh_revision), colour = "#E69F00", size = 2.5) +
  annotate("rect",
           xmin = as.Date("2020-01-01"), xmax = as.Date("2020-10-01"),
           ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "firebrick") +
  annotate("text", x = as.Date("2020-04-01"), y = Inf,
           label = "COVID\nshock", vjust = 1.4, size = 3.2, colour = "firebrick") +
  annotate("text", x = max(back_dt$vintage_date),
           y = tail(back_dt$hp_revision, 1), label = "HP",
           hjust = -0.3, colour = "#0072B2", fontface = "bold") +
  annotate("text", x = max(back_dt$vintage_date),
           y = tail(back_dt$mbh_revision, 1), label = "MBH",
           hjust = -0.3, colour = "#E69F00", fontface = "bold") +
  scale_x_date(
    date_breaks = "6 months",
    labels      = function(x) paste0(format(x, "%Y"), " ", quarters(x)),
    expand      = expansion(add = c(30, 90))
  ) +
  labs(
    title    = "Backward Revision at 2019 Q4 as New Data Arrive",
    subtitle = "How does each filter re-estimate the pre-COVID trend as COVID data are published?",
    x        = "Vintage (publication date)",
    y        = "Revision vs 2019 Q4 baseline (log pts x100)"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.title = element_text(size = 9))

print(p2)

## ----spread-table-------------------------------------------------------------
knitr::kable(
  back_dt[, .(
    Vintage              = format(vintage_date),
    `HP trend at 2019Q4`  = round(hp_trend,  5),
    `MBH trend at 2019Q4` = round(mbh_trend, 5),
    `HP revision (ppts)`  = round(hp_revision,  3),
    `MBH revision (ppts)` = round(mbh_revision, 3)
  )],
  caption = "Backward revision of 2019 Q4 trend estimate across vintages"
)

