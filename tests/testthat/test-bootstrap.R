# ── Tests: block-bootstrap confidence bands (all four filters) ────────────────

# Helpers: each entry runs a filter with bootstrap on a quarterly ts and on a
# plain numeric vector. MBH uses a small mstop/boot_iter (we test band
# STRUCTURE, not fit quality) to keep the suite fast.

make_ts  <- function(n = 120L) {
  set.seed(123)
  ts(cumsum(rnorm(n)) + 20, start = c(2000, 1), frequency = 4)
}

boot_runners <- list(
  HP       = function(x, bi) hp_filter(x, boot_iter = bi),
  Hamilton = function(x, bi) hamilton_filter(x, boot_iter = bi),
  bHP      = function(x, bi) bhp_filter(x, boot_iter = bi, iter_max = 5L,
                                        stopping = "fixed"),
  MBH      = function(x, bi) suppressMessages(
    mbh_filter(x, boot_iter = bi, mstop = 50L)
  )
)

# ── 1. No bootstrap -> no bands ──────────────────────────────────────────────

test_that("boot_iter = 0 produces no confidence bands", {
  y <- make_ts()
  for (nm in names(boot_runners)) {
    fit <- boot_runners[[nm]](y, 0L)
    expect_null(fit$trend_lower, info = nm)
    expect_null(fit$trend_upper, info = nm)
  }
})

# ── 2. Bootstrap -> well-formed bands ────────────────────────────────────────

test_that("boot_iter > 0 yields length-correct, ordered, containing bands", {
  y <- make_ts()
  n <- length(y)
  for (nm in names(boot_runners)) {
    fit <- boot_runners[[nm]](y, 20L)

    expect_false(is.null(fit$trend_lower), info = nm)
    expect_false(is.null(fit$trend_upper), info = nm)
    expect_length(fit$trend_lower, n)
    expect_length(fit$trend_upper, n)

    ok <- !is.na(fit$trend_lower) & !is.na(fit$trend_upper)
    expect_true(all(fit$trend_lower[ok] <= fit$trend_upper[ok]), info = nm)

    # the point estimate must stay inside its own band (centred by construction)
    inside <- fit$trend[ok] >= fit$trend_lower[ok] - 1e-6 &
              fit$trend[ok] <= fit$trend_upper[ok] + 1e-6
    expect_true(all(inside), info = nm)
  }
})

# ── 3. S3 contract ───────────────────────────────────────────────────────────

test_that("all filters return class c('macrofilter', 'list')", {
  y <- make_ts()
  for (nm in names(boot_runners)) {
    fit <- boot_runners[[nm]](y, 10L)
    expect_identical(class(fit), c("macrofilter", "list"), info = nm)
  }
})

test_that("meta stores temporal identity for ts input", {
  y <- make_ts()
  for (nm in names(boot_runners)) {
    fit <- boot_runners[[nm]](y, 10L)
    expect_identical(fit$meta$ts_class, "ts", info = nm)
    expect_equal(fit$meta$tsp[3L], 4, info = nm)
  }
})

# ── 4. Hamilton lead-in convention ───────────────────────────────────────────

test_that("Hamilton bands are NA exactly over the h + p - 1 lead-in", {
  y <- make_ts()
  fit <- hamilton_filter(y, boot_iter = 20L)   # auto h = 8 (quarterly), p = 4
  lead <- fit$meta$h + fit$meta$p - 1L

  expect_equal(sum(is.na(fit$trend_lower)), lead)
  expect_true(all(is.na(fit$trend_lower[seq_len(lead)])))
  expect_false(anyNA(fit$trend_lower[(lead + 1L):length(y)]))
})

# ── 5. Works on plain numeric input too ──────────────────────────────────────

test_that("bootstrap bands work on plain numeric input", {
  set.seed(7)
  v <- cumsum(rnorm(120)) + 20
  fit <- suppressWarnings(hp_filter(v, boot_iter = 20L))
  expect_length(fit$trend_lower, length(v))
  expect_null(fit$meta$tsp)                     # no temporal identity
})
