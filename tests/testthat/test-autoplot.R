# ── Tests: autoplot.macrofilter() ─────────────────────────────────────────────

geoms_of <- function(p) vapply(p$layers, function(l) class(l$geom)[1L],
                               character(1L))

# ── 1. Returns a ggplot ──────────────────────────────────────────────────────

test_that("autoplot returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  y <- ts(cumsum(rnorm(80)) + 10, start = c(2000, 1), frequency = 4)
  p <- ggplot2::autoplot(hp_filter(y))
  expect_s3_class(p, "ggplot")
})

# ── 2. Input validation ──────────────────────────────────────────────────────

test_that("autoplot rejects non-macrofilter input", {
  skip_if_not_installed("ggplot2")
  # call the method directly: dispatch would never route a plain list here
  expect_error(autoplot.macrofilter(list(a = 1)), "macrofilter")
})

# ── 3. Ribbon present iff bands exist ─────────────────────────────────────────

test_that("ribbon layer appears only when bands are present", {
  skip_if_not_installed("ggplot2")
  y <- ts(cumsum(rnorm(80)) + 10, start = c(2000, 1), frequency = 4)

  p_plain <- ggplot2::autoplot(hp_filter(y))
  expect_false("GeomRibbon" %in% geoms_of(p_plain))

  p_band <- ggplot2::autoplot(hp_filter(y, boot_iter = 20L))
  expect_true("GeomRibbon" %in% geoms_of(p_band))
})

# ── 4. X-axis reconstruction ─────────────────────────────────────────────────

test_that("time axis is reconstructed from tsp when available", {
  skip_if_not_installed("ggplot2")

  # ts input: decimal time rebuilt from meta$tsp, starting at 2000
  y_ts <- ts(cumsum(rnorm(80)), start = c(2000, 1), frequency = 4)
  p_ts <- ggplot2::autoplot(hp_filter(y_ts))
  expect_equal(p_ts$data$time[1L], 2000)
  expect_equal(p_ts$labels$x, "Time")          # label hardcoded by design

  # plain numeric: positional index 1..n
  v <- cumsum(rnorm(80))
  p_v <- ggplot2::autoplot(suppressWarnings(hp_filter(v)))
  expect_equal(p_v$data$time, seq_len(80))
})

# ── 5. Subtitle only with confidence bands ───────────────────────────────────

test_that("CI subtitle appears only when bands are present", {
  skip_if_not_installed("ggplot2")
  y <- ts(cumsum(rnorm(80)) + 10, start = c(2000, 1), frequency = 4)

  expect_null(ggplot2::autoplot(hp_filter(y))$labels$subtitle)
  expect_match(
    ggplot2::autoplot(hp_filter(y, boot_iter = 20L))$labels$subtitle,
    "Bootstrap Confidence Intervals"
  )
})
