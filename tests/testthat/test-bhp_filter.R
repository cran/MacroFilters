# ── Tests: bhp_filter() ────────────────────────────────────────────────────────

# ── 1. Identity: trend + cycle == data ───────────────────────────────────────

test_that("Identity holds with stopping = 'fixed'", {
  set.seed(1)
  y <- cumsum(rnorm(100))
  result <- bhp_filter(y, lambda = 1600, iter_max = 10L, stopping = "fixed")
  expect_equal(result$trend + result$cycle, result$data, tolerance = 1e-9)
})

test_that("Identity holds with stopping = 'adf'", {
  set.seed(2)
  y <- cumsum(rnorm(100))
  result <- bhp_filter(y, lambda = 1600, stopping = "adf")
  expect_equal(result$trend + result$cycle, result$data, tolerance = 1e-9)
})

test_that("Identity holds with stopping = 'bic'", {
  set.seed(3)
  y <- cumsum(rnorm(100))
  result <- bhp_filter(y, lambda = 1600, stopping = "bic")
  expect_equal(result$trend + result$cycle, result$data, tolerance = 1e-9)
})

# ── 2. Invariance: same result regardless of input class ─────────────────────

test_that("bHP produces identical results across numeric and ts", {
  set.seed(42)
  y_num <- cumsum(rnorm(100))
  y_ts  <- ts(y_num, start = c(2000, 1), frequency = 4)

  res_num <- bhp_filter(y_num, lambda = 1600, iter_max = 5L, stopping = "fixed")
  res_ts  <- bhp_filter(y_ts, lambda = 1600, iter_max = 5L, stopping = "fixed")

  expect_equal(res_num$trend, res_ts$trend, tolerance = 1e-12)
  expect_equal(res_num$cycle, res_ts$cycle, tolerance = 1e-12)
})

test_that("bHP produces identical results with xts input", {
  skip_if_not_installed("xts")
  set.seed(42)
  y_num <- cumsum(rnorm(100))
  dates <- seq(as.Date("2000-01-01"), by = "quarter", length.out = 100)
  y_xts <- xts::xts(y_num, order.by = dates)

  res_num <- bhp_filter(y_num, lambda = 1600, iter_max = 5L, stopping = "fixed")
  res_xts <- bhp_filter(y_xts, lambda = 1600, iter_max = 5L, stopping = "fixed")

  expect_equal(res_num$trend, res_xts$trend, tolerance = 1e-12)
  expect_equal(res_num$cycle, res_xts$cycle, tolerance = 1e-12)
})

test_that("bHP produces identical results with zoo input", {
  skip_if_not_installed("zoo")
  set.seed(42)
  y_num <- cumsum(rnorm(100))
  dates <- seq(as.Date("2000-01-01"), by = "quarter", length.out = 100)
  y_zoo <- zoo::zoo(y_num, order.by = dates)

  res_num <- bhp_filter(y_num, lambda = 1600, iter_max = 5L, stopping = "fixed")
  res_zoo <- bhp_filter(y_zoo, lambda = 1600, iter_max = 5L, stopping = "fixed")

  expect_equal(res_num$trend, res_zoo$trend, tolerance = 1e-12)
  expect_equal(res_num$cycle, res_zoo$cycle, tolerance = 1e-12)
})

# ── 3. Stress tests ──────────────────────────────────────────────────────────

test_that("bHP errors when n < 3", {
  expect_error(bhp_filter(c(1, 2), lambda = 1600), "at least 3")
  expect_error(bhp_filter(1, lambda = 1600), "at least 3")
})

test_that("bHP works with minimum-length series (n = 3, 4, 5)", {
  for (n in 3:5) {
    y <- seq_len(n) * 1.0
    result <- bhp_filter(y, lambda = 1600, iter_max = 3L, stopping = "fixed")
    expect_equal(result$trend + result$cycle, result$data, tolerance = 1e-9)
  }
})

test_that("bHP handles long series with stopping = 'fixed'", {
  set.seed(99)
  y <- cumsum(rnorm(10000))
  result <- bhp_filter(y, lambda = 1600, iter_max = 5L, stopping = "fixed")
  expect_equal(result$trend + result$cycle, result$data, tolerance = 1e-9)
})

test_that("iter_max = 1 with fixed stopping equals hp_filter", {
  set.seed(10)
  y <- cumsum(rnorm(100))
  res_bhp <- bhp_filter(y, lambda = 1600, iter_max = 1L, stopping = "fixed")
  res_hp  <- hp_filter(y, lambda = 1600)
  expect_equal(res_bhp$trend, res_hp$trend, tolerance = 1e-12)
  expect_equal(res_bhp$cycle, res_hp$cycle, tolerance = 1e-12)
})

test_that("ADF stops before iter_max on stationary series", {
  set.seed(11)
  # White noise is stationary; ADF should reject unit root quickly
  y <- rnorm(200)
  result <- bhp_filter(y, lambda = 1600, iter_max = 100L, stopping = "adf")
  expect_lt(result$meta$iterations, 100L)
})

test_that("BIC selects reasonable iteration count", {
  set.seed(12)
  y <- cumsum(rnorm(200))
  result <- bhp_filter(y, lambda = 1600, iter_max = 100L, stopping = "bic")
  expect_gte(result$meta$iterations, 1L)
  expect_lte(result$meta$iterations, 100L)
})

test_that("Auto-lambda from ts frequency works correctly", {
  y_ts <- ts(1:50 * 1.0, frequency = 4)
  result <- bhp_filter(y_ts, iter_max = 3L, stopping = "fixed")
  expect_equal(result$meta$lambda, 6.25 * 4^4)
})

test_that("Warning when frequency cannot be determined", {
  y <- cumsum(rnorm(50))
  expect_warning(bhp_filter(y, iter_max = 1L, stopping = "fixed"),
                 "Cannot determine series frequency")
})

# ── 4. Output format ─────────────────────────────────────────────────────────

test_that("bHP returns macrofilter with correct meta fields", {
  set.seed(7)
  y <- cumsum(rnorm(100))
  result <- bhp_filter(y, lambda = 1600, iter_max = 5L, stopping = "fixed")

  expect_s3_class(result, "macrofilter")
  expect_equal(result$meta$method, "bHP")
  expect_equal(result$meta$lambda, 1600)
  expect_equal(result$meta$stopping_rule, "fixed")
  expect_true(is.integer(result$meta$iterations))
  expect_true(is.numeric(result$meta$compute_time))
})

test_that("iterations is integer and stopping_rule is character", {
  set.seed(8)
  y <- cumsum(rnorm(100))
  result <- bhp_filter(y, lambda = 1600, stopping = "bic")

  expect_type(result$meta$iterations, "integer")
  expect_type(result$meta$stopping_rule, "character")
})
