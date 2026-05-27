# ── Tests: hamilton_filter() ───────────────────────────────────────────────────

# ── 1. Identity: trend + cycle == data (non-NA positions) ─────────────────────

test_that("Hamilton trend + cycle reconstructs original data", {
  set.seed(1)
  y <- cumsum(rnorm(100))
  result <- hamilton_filter(y, h = 8, p = 4)
  valid <- !is.na(result$trend)
  expect_equal(
    result$trend[valid] + result$cycle[valid],
    result$data[valid],
    tolerance = 1e-9
  )
})

test_that("Identity holds for ts input with auto-h", {
  set.seed(2)
  y_ts <- ts(cumsum(rnorm(80)), start = c(2000, 1), frequency = 4)
  result <- hamilton_filter(y_ts)
  valid <- !is.na(result$trend)
  expect_equal(
    result$trend[valid] + result$cycle[valid],
    result$data[valid],
    tolerance = 1e-9
  )
})

# ── 2. Invariance: same result regardless of input class ─────────────────────

test_that("Hamilton filter produces identical results across input classes", {
  set.seed(42)
  y_num <- cumsum(rnorm(100))
  y_ts  <- ts(y_num, start = c(2000, 1), frequency = 4)

  res_num <- hamilton_filter(y_num, h = 8, p = 4)
  res_ts  <- hamilton_filter(y_ts, h = 8, p = 4)

  expect_equal(res_num$trend, res_ts$trend, tolerance = 1e-12)
  expect_equal(res_num$cycle, res_ts$cycle, tolerance = 1e-12)
})

test_that("Hamilton filter produces identical results with xts input", {
  skip_if_not_installed("xts")
  set.seed(42)
  y_num <- cumsum(rnorm(100))
  dates <- seq(as.Date("2000-01-01"), by = "quarter", length.out = 100)
  y_xts <- xts::xts(y_num, order.by = dates)

  res_num <- hamilton_filter(y_num, h = 8, p = 4)
  res_xts <- hamilton_filter(y_xts, h = 8, p = 4)

  expect_equal(res_num$trend, res_xts$trend, tolerance = 1e-12)
  expect_equal(res_num$cycle, res_xts$cycle, tolerance = 1e-12)
})

test_that("Hamilton filter produces identical results with zoo input", {
  skip_if_not_installed("zoo")
  set.seed(42)
  y_num <- cumsum(rnorm(100))
  dates <- seq(as.Date("2000-01-01"), by = "quarter", length.out = 100)
  y_zoo <- zoo::zoo(y_num, order.by = dates)

  res_num <- hamilton_filter(y_num, h = 8, p = 4)
  res_zoo <- hamilton_filter(y_zoo, h = 8, p = 4)

  expect_equal(res_num$trend, res_zoo$trend, tolerance = 1e-12)
  expect_equal(res_num$cycle, res_zoo$cycle, tolerance = 1e-12)
})

# ── 3. Stress tests ──────────────────────────────────────────────────────────

test_that("Hamilton filter errors when series is too short", {
  y_short <- 1:12 * 1.0  # n = 12, need > h + p = 12
  expect_error(hamilton_filter(y_short, h = 8, p = 4), "too short")
  expect_error(hamilton_filter(1:5 * 1.0, h = 2, p = 4), "too short")
})

test_that("Hamilton filter works at minimum viable length", {
  # n = 13, h = 8, p = 4 -> n > h + p = 12 -> n_valid = n - h - p + 1 = 2
  y <- 1:13 * 1.0
  result <- hamilton_filter(y, h = 8, p = 4)
  valid <- !is.na(result$trend)
  expect_equal(sum(valid), 2L)
  expect_equal(
    result$trend[valid] + result$cycle[valid],
    result$data[valid],
    tolerance = 1e-9
  )
})

test_that("Hamilton filter handles long series (n = 10000)", {
  set.seed(99)
  y <- cumsum(rnorm(10000))
  result <- hamilton_filter(y, h = 8, p = 4)
  valid <- !is.na(result$trend)
  expect_equal(
    result$trend[valid] + result$cycle[valid],
    result$data[valid],
    tolerance = 1e-9
  )
})

test_that("First h + p - 1 positions are NA, rest are computed", {
  set.seed(10)
  y <- cumsum(rnorm(100))
  h <- 8L
  p <- 4L
  result <- hamilton_filter(y, h = h, p = p)

  n_na <- h + p - 1L
  expect_true(all(is.na(result$trend[1:n_na])))
  expect_true(all(is.na(result$cycle[1:n_na])))
  expect_true(all(!is.na(result$trend[(n_na + 1):length(y)])))
  expect_true(all(!is.na(result$cycle[(n_na + 1):length(y)])))
})

test_that("Auto-h for quarterly ts gives h = 8", {
  y_ts <- ts(1:50 * 1.0, frequency = 4)
  result <- hamilton_filter(y_ts)
  expect_equal(result$meta$h, 8L)
})

test_that("Auto-h for monthly ts gives h = 24", {
  y_ts <- ts(1:100 * 1.0, frequency = 12)
  result <- hamilton_filter(y_ts)
  expect_equal(result$meta$h, 24L)
})

test_that("Auto-h for annual ts gives h = 2", {
  y_ts <- ts(1:30 * 1.0, frequency = 1)
  result <- hamilton_filter(y_ts)
  expect_equal(result$meta$h, 2L)
})

test_that("Warning emitted when frequency cannot be determined", {
  y <- cumsum(rnorm(50))
  expect_warning(hamilton_filter(y), "Cannot determine series frequency")
})

# ── 4. Output format ─────────────────────────────────────────────────────────

test_that("Hamilton filter returns macrofilter object with correct meta", {
  set.seed(7)
  y <- cumsum(rnorm(100))
  result <- hamilton_filter(y, h = 8, p = 4)

  expect_s3_class(result, "macrofilter")
  expect_equal(result$meta$method, "Hamilton")
  expect_equal(result$meta$h, 8L)
  expect_equal(result$meta$p, 4L)
  expect_true(!is.null(result$meta$coefficients))
  expect_length(result$meta$coefficients, 5L)  # intercept + 4 lags
  expect_true(!is.null(result$meta$compute_time))
  expect_true(is.numeric(result$meta$compute_time))
})

test_that("Lengths of trend, cycle, and data all match", {
  set.seed(3)
  y <- cumsum(rnorm(100))
  result <- hamilton_filter(y, h = 8, p = 4)
  expect_equal(length(result$trend), length(result$data))
  expect_equal(length(result$cycle), length(result$data))
})
