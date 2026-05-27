# ── Tests: mbh_filter() ────────────────────────────────────────────────────────

# ── 1. Identity: trend + cycle == data ───────────────────────────────────────

test_that("Identity holds with explicit params", {
  set.seed(1)
  y <- cumsum(rnorm(100))
  result <- mbh_filter(y, knots = 20L, mstop = 50L, d = 0.01)
  expect_equal(result$trend + result$cycle, result$data, tolerance = 1e-9)
})

test_that("Identity holds with auto-knots", {
  set.seed(2)
  y <- cumsum(rnorm(200))
  result <- mbh_filter(y)
  expect_equal(result$trend + result$cycle, result$data, tolerance = 1e-9)
})

# ── 2. Invariance: same result regardless of input class ─────────────────────

test_that("MBH produces identical results across numeric and ts", {
  set.seed(42)
  y_num <- cumsum(rnorm(100))
  y_ts  <- ts(y_num, start = c(2000, 1), frequency = 4)

  res_num <- mbh_filter(y_num, knots = 20L, mstop = 50L, d = 0.01)
  res_ts  <- mbh_filter(y_ts,  knots = 20L, mstop = 50L, d = 0.01)

  expect_equal(res_num$trend, res_ts$trend, tolerance = 1e-12)
  expect_equal(res_num$cycle, res_ts$cycle, tolerance = 1e-12)
})

test_that("MBH produces identical results with xts input", {
  skip_if_not_installed("xts")
  set.seed(42)
  y_num <- cumsum(rnorm(100))
  dates <- seq(as.Date("2000-01-01"), by = "quarter", length.out = 100)
  y_xts <- xts::xts(y_num, order.by = dates)

  res_num <- mbh_filter(y_num, knots = 20L, mstop = 50L, d = 0.01)
  res_xts <- mbh_filter(y_xts, knots = 20L, mstop = 50L, d = 0.01)

  expect_equal(res_num$trend, res_xts$trend, tolerance = 1e-12)
  expect_equal(res_num$cycle, res_xts$cycle, tolerance = 1e-12)
})

test_that("MBH produces identical results with zoo input", {
  skip_if_not_installed("zoo")
  set.seed(42)
  y_num <- cumsum(rnorm(100))
  dates <- seq(as.Date("2000-01-01"), by = "quarter", length.out = 100)
  y_zoo <- zoo::zoo(y_num, order.by = dates)

  res_num <- mbh_filter(y_num, knots = 20L, mstop = 50L, d = 0.01)
  res_zoo <- mbh_filter(y_zoo, knots = 20L, mstop = 50L, d = 0.01)

  expect_equal(res_num$trend, res_zoo$trend, tolerance = 1e-12)
  expect_equal(res_num$cycle, res_zoo$cycle, tolerance = 1e-12)
})

# ── 3. Stress tests ──────────────────────────────────────────────────────────

test_that("MBH errors when n < 6", {
  expect_error(mbh_filter(c(1, 2, 3, 4, 5)), "at least 6")
  expect_error(mbh_filter(c(1, 2, 3)),        "at least 6")
  expect_error(mbh_filter(1),                  "at least 6")
})

test_that("MBH works with minimum-length series (n = 6, 10, 20)", {
  for (n in c(6L, 10L, 20L)) {
    y <- seq_len(n) * 1.0
    result <- mbh_filter(y, mstop = 10L)
    expect_equal(result$trend + result$cycle, result$data, tolerance = 1e-9)
  }
})

test_that("MBH handles long series (n = 2000)", {
  set.seed(99)
  y <- cumsum(rnorm(2000))
  result <- mbh_filter(y, mstop = 50L)
  expect_equal(result$trend + result$cycle, result$data, tolerance = 1e-9)
})

test_that("Auto-knots heuristic: n=100 gives knots = max(20, 50) = 50", {
  set.seed(50)
  y <- cumsum(rnorm(100))
  result <- mbh_filter(y, mstop = 10L)
  expect_equal(result$meta$knots, 50L)
})

test_that("Custom d parameter is stored correctly in meta", {
  set.seed(51)
  y <- cumsum(rnorm(100))
  result <- mbh_filter(y, d = 0.5, mstop = 10L)
  expect_equal(result$meta$d, 0.5)
})

test_that("Auto-d computes correctly for normal series", {
  set.seed(10)
  y <- cumsum(rnorm(100))
  result <- mbh_filter(y, mstop = 20L)
  expect_true(is.numeric(result$meta$d))
  expect_gt(result$meta$d, 0)
})

test_that("Custom d overrides Auto-d", {
  set.seed(11)
  y <- cumsum(rnorm(100))
  result <- mbh_filter(y, d = 5.5, mstop = 20L)
  expect_equal(result$meta$d, 5.5)
})

test_that("Auto-d falls back to 0.01 when MAD of diffs is zero", {
  # Perfectly linear series: diff(y) is constant, so mad(diff(y)) == 0
  y <- seq(1, 20) * 1.0
  result <- mbh_filter(y, mstop = 10L)
  expect_equal(result$meta$d, 0.01)
})

# ── 4. Output format ─────────────────────────────────────────────────────────

test_that("MBH returns macrofilter with correct meta$method = 'MBH'", {
  set.seed(7)
  y <- cumsum(rnorm(100))
  result <- mbh_filter(y, mstop = 20L)

  expect_s3_class(result, "macrofilter")
  expect_equal(result$meta$method, "MBH")
})

test_that("Meta contains all params: knots, d, mstop, nu, compute_time", {
  set.seed(8)
  y <- cumsum(rnorm(100))
  result <- mbh_filter(y, knots = 15L, mstop = 30L, d = 0.02, nu = 0.05)

  expect_equal(result$meta$knots, 15L)
  expect_equal(result$meta$d, 0.02)
  expect_equal(result$meta$mstop, 30L)
  expect_equal(result$meta$nu, 0.05)
  expect_true(is.numeric(result$meta$compute_time))
})
