#' US Real GDP — FRED Vintage
#'
#' Quarterly US Real Gross Domestic Product from the Federal Reserve Bank of
#' St. Louis (FRED) public data API (series **GDPC1**), expressed in billions
#' of chained 2017 US dollars, seasonally adjusted annual rate.
#'
#' @format A `data.table` with one row per quarter and three columns:
#' \describe{
#'   \item{`date`}{`Date`. Quarter start date (e.g. `1947-01-01` = 1947 Q1).}
#'   \item{`gdp_real`}{`numeric`. Real GDP level, billions of chained 2017 USD.}
#'   \item{`gdp_log`}{`numeric`. Natural logarithm of `gdp_real`, pre-computed
#'     for convenience.}
#' }
#'
#' @details
#' The dataset covers 1947 Q1 through the latest vintage available at download
#' time (approximately 316 rows as of 2025). Rows with `NA` in `gdp_real` are
#' excluded.
#'
#' The log-level column (`gdp_log`) is particularly useful for trend-cycle
#' decomposition because log-differences approximate quarter-on-quarter
#' percentage growth rates:
#' \deqn{\Delta \log(\text{GDP}_t) \approx g_t}
#'
#' @source
#' Federal Reserve Bank of St. Louis — FRED Economic Data,
#' series GDPC1. Downloaded via the public CSV endpoint
#' `https://fred.stlouisfed.org/graph/fredgraph.csv?id=GDPC1`.
#' See `data-raw/us_gdp_vintage.R` for the reproducible download script.
#'
#' @examples
#' data("us_gdp_vintage", package = "MacroFilters")
#' head(us_gdp_vintage)
#' plot(us_gdp_vintage$date, us_gdp_vintage$gdp_log,
#'      type = "l", xlab = "Date", ylab = "Log Real GDP",
#'      main = "US Real GDP (log level)")
"us_gdp_vintage"
