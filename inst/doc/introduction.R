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

## ----input-agnosticism--------------------------------------------------------
set.seed(7)
y_raw <- cumsum(rnorm(60)) + (1:60) * 0.3   # a simple integrated series

# As plain numeric
hp_num <- hp_filter(y_raw)
class(hp_num$trend)   # numeric

# As a monthly ts object
y_ts   <- ts(y_raw, start = c(2019, 1), frequency = 12)
hp_ts  <- hp_filter(y_ts)
class(hp_ts$trend)    # ts — output matches input

## ----hp-filter----------------------------------------------------------------
set.seed(42)
n  <- 100
y  <- ts(100 + 0.4 * (1:n) + 5 * sin(2 * pi * (1:n) / 20) + rnorm(n, sd = 2),
         start = c(2000, 1), frequency = 4)

hp <- hp_filter(y)
hp

## ----hp-plot, echo=FALSE, fig.height=3----------------------------------------
df_hp <- data.frame(
  t     = as.numeric(time(y)),
  data  = as.numeric(y),
  trend = as.numeric(hp$trend),
  cycle = as.numeric(hp$cycle)
)

ggplot(df_hp, aes(x = t)) +
  geom_line(aes(y = data,  colour = "Observed"), linewidth = 0.6, linetype = "dashed") +
  geom_line(aes(y = trend, colour = "Trend"),    linewidth = 1.1) +
  scale_colour_manual(values = c("Observed" = "grey60", "Trend" = "#0072B2")) +
  labs(title = "HP Filter", subtitle = paste0("\u03bb = ", hp$meta$lambda),
       x = "Year", y = "Value", colour = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

## ----hp-cycle-plot, echo=FALSE, fig.height=2.8--------------------------------
ggplot(df_hp, aes(x = t, y = cycle)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_line(colour = "#0072B2", linewidth = 0.8) +
  labs(title = "HP Cycle", x = "Year", y = "Cycle") +
  theme_minimal(base_size = 12)

## ----hamilton-filter----------------------------------------------------------
ham <- hamilton_filter(y)      # auto-detects h = 8 for quarterly
ham

## ----bhp-filter---------------------------------------------------------------
bhp <- bhp_filter(y, stopping = "bic")
bhp

## ----mbh-demo-----------------------------------------------------------------
set.seed(42)
n <- 80
t <- 1:n

# 1. Define simulation parameters
sd_noise <- 1.8
trend    <- 100 + 0.3 * t + 0.005 * t^2
cycle    <- 2.5 * sin(2 * pi * t / 28) # 7-year business cycle

# 2. Simulate GDP index
gdp_num <- trend + cycle + rnorm(n, sd = sd_noise)

# 3. Inject a structural shock (e.g., COVID-19 lockdown)
# Expressed as extreme standard deviation events
gdp_num[60] <- gdp_num[60] - (16 * sd_noise) # Massive crash
gdp_num[61] <- gdp_num[61] - (9 * sd_noise)  # Partial recovery
gdp_num[62] <- gdp_num[62] - (4 * sd_noise)  # V Stabilization
gdp_num[62] <- gdp_num[62] - (2 * sd_noise)  # Stabilization

gdp <- ts(gdp_num, start = c(2001, 1), frequency = 4)

# Extract trends
hp_res  <- hp_filter(gdp)

# MBH Filter: Auto-calibrated threshold (d) based on MAD of differences
mbh_res <- mbh_filter(gdp) 

mbh_res

## ----mbh-plot, echo=FALSE-----------------------------------------------------
df_mbh <- data.frame(
  t    = as.numeric(time(gdp)),
  data = as.numeric(gdp),
  hp   = as.numeric(hp_res$trend),
  mbh  = as.numeric(mbh_res$trend)
)

ggplot(df_mbh, aes(x = t)) +
  geom_line(aes(y = data, colour = "Observed"),   linewidth = 0.6, linetype = "dashed") +
  geom_line(aes(y = hp,   colour = "HP trend"),   linewidth = 1.0) +
  geom_line(aes(y = mbh,  colour = "MBH trend"),  linewidth = 1.1) +
  annotate("rect",
           xmin = df_mbh$t[59], xmax = df_mbh$t[63],
           ymin = -Inf, ymax = Inf,
           alpha = 0.12, fill = "firebrick") +
  annotate("text",
           x = df_mbh$t[61], y = max(df_mbh$data, na.rm = TRUE),
           label = "COVID shock", vjust = -0.5,
           size  = 3.5, colour = "firebrick") +
  scale_colour_manual(
    values = c("Observed" = "grey60", "HP trend" = "#0072B2", "MBH trend" = "#E69F00")
  ) +
  labs(
    title    = "HP vs MBH under a Structural Shock",
    subtitle = "MBH trend (orange) stays smooth; HP trend (blue) is pulled down",
    x = "Year", y = "GDP Index", colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

## ----s3-print-----------------------------------------------------------------
mbh_res

## ----s3-access----------------------------------------------------------------
# Trend and cycle as plain vectors
head(mbh_res$trend, 8)
head(mbh_res$cycle, 8)

# Verify the fundamental identity: trend + cycle == data
max(abs((mbh_res$trend + mbh_res$cycle) - mbh_res$data))  # should be < 1e-9

## ----s3-meta------------------------------------------------------------------
str(mbh_res$meta)

## ----cycle-comparison, fig.height=5-------------------------------------------
df_cycle <- data.frame(
  t          = as.numeric(time(gdp)),
  HP_cycle   = as.numeric(hp_res$cycle),
  MBH_cycle  = as.numeric(mbh_res$cycle)
)

ggplot(df_cycle, aes(x = t)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_line(aes(y = HP_cycle,  colour = "HP cycle"),  linewidth = 0.8) +
  geom_line(aes(y = MBH_cycle, colour = "MBH cycle"), linewidth = 0.8) +
  annotate("rect",
           xmin = df_cycle$t[59], xmax = df_cycle$t[63],
           ymin = -Inf, ymax = Inf,
           alpha = 0.12, fill = "firebrick") +
  scale_colour_manual(values = c("HP cycle" = "#0072B2", "MBH cycle" = "#E69F00")) +
  labs(
    title    = "Cyclical Components",
    subtitle = "HP cycle absorbs the shock; MBH cycle faithfully records it",
    x = "Year", y = "Cycle", colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

