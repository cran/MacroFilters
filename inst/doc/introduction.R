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
data("es_gdp", package = "MacroFilters")

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

## ----fred-download, eval=FALSE------------------------------------------------
# # FRED public endpoint — no API key needed.
# # See data-raw/intl_gdp.R for the full reproducible download script.
# read_fred <- function(id) {
#   url <- sprintf("https://fred.stlouisfed.org/graph/fredgraph.csv?id=%s", id)
#   dt  <- read.csv(url, col.names = c("date", "gdp_real"), na.strings = ".")
#   dt$date    <- as.Date(dt$date)
#   dt$gdp_log <- log(as.numeric(dt$gdp_real))
#   dt[!is.na(dt$gdp_real), ]
# }
# fr_raw <- read_fred("CLVMNACSCAB1GQFR")
# es_raw <- read_fred("CLVMNACSCAB1GQES")

## ----mbh-demo-----------------------------------------------------------------
# Apply HP + MBH per country.
# For log-level series, auto d (MAD of diff) is too tight — calibrate d on the
# cycle scale instead (see vignette "Hyperparameter Tuning for the MBH Filter").
make_trend_df <- function(raw, country) {
  dt  <- raw[raw$date >= as.Date("2000-01-01"), ]
  g   <- ts(dt$gdp_log, start = c(2000, 1), frequency = 4)
  hp  <- hp_filter(g)
  mbh <- mbh_filter(g, d = mad(hp$cycle))
  data.frame(country  = country,
             t        = as.numeric(time(g)),
             observed = as.numeric(g),
             hp       = as.numeric(hp$trend),
             mbh      = as.numeric(mbh$trend))
}

df_plot <- rbind(
  make_trend_df(fr_gdp, "France"),
  make_trend_df(es_gdp, "Spain")
)

# Keep Spain filter objects for the S3 class examples in Section 5
dt_es   <- es_gdp[es_gdp$date >= as.Date("2000-01-01"), ]
gdp     <- ts(dt_es$gdp_log, start = c(2000, 1), frequency = 4)
hp_res  <- hp_filter(gdp)
mbh_res <- mbh_filter(gdp, d = mad(hp_res$cycle))

mbh_res

## ----mbh-plot, echo=FALSE-----------------------------------------------------
ggplot(df_plot, aes(x = t)) +
  geom_line(aes(y = observed, colour = "Observed"),  linewidth = 0.6, linetype = "dashed") +
  geom_line(aes(y = hp,       colour = "HP trend"),  linewidth = 1.0) +
  geom_line(aes(y = mbh,      colour = "MBH trend"), linewidth = 1.1) +
  annotate("rect",
           xmin = 2020.00, xmax = 2020.75,
           ymin = -Inf, ymax = Inf,
           alpha = 0.12, fill = "firebrick") +
  annotate("text",
           x = 2020.375, y = Inf, vjust = 1.4,
           label = "COVID-19\n2020 Q2", size = 3.2, colour = "firebrick") +
  scale_colour_manual(
    values = c("Observed" = "grey60", "HP trend" = "#0072B2", "MBH trend" = "#E69F00")
  ) +
  facet_wrap(~country, scales = "free_y") +
  labs(
    title    = "HP vs MBH under a Structural Shock",
    subtitle = "MBH trend (orange) stays smooth;\nHP trend (blue) is pulled down by the COVID crash",
    x = "Year", y = "Log Real GDP", colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

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
           xmin = 2020.00, xmax = 2020.75,
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

