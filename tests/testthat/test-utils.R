test_that("enu <-> lonlat round-trips within a centimetre", {
  x <- c(0, 1500, -800, 12000); y <- c(0, -300, 7000, 12000)
  ll <- enu_to_lonlat(x, y, -108.5, 43.2)
  back <- lonlat_to_enu(ll$lon, ll$lat, -108.5, 43.2)
  expect_equal(back$x, x, tolerance = 1e-6)
  expect_equal(back$y, y, tolerance = 1e-6)
  # one degree of latitude near 43N is ~111.1 km
  expect_equal(enu_to_lonlat(0, 111100, 0, 43.2)$lat, 44.2, tolerance = 1e-3)
})

test_that("fourth difference annihilates cubics and scales x^4 correctly", {
  i <- 1:50
  expect_equal(fourth_difference(2 + 3 * i - i^2 + 0.5 * i^3)[3:48], rep(0, 46), tolerance = 1e-8)
  d4 <- fourth_difference(i^4, normalize = FALSE)
  expect_equal(d4[3:48], rep(24, 46), tolerance = 1e-8)
  expect_equal(fourth_difference(i^4)[10], 24 / 16)
  expect_true(all(is.na(fourth_difference(i^4)[c(1, 2, 49, 50)])))
  expect_true(all(is.na(fourth_difference(1:4))))
})

test_that(".runs collapses logical vectors into intervals", {
  r <- magqc:::.runs(c(FALSE, TRUE, TRUE, FALSE, NA, TRUE))
  expect_equal(r$start, c(2, 6))
  expect_equal(r$end, c(3, 6))
  expect_equal(nrow(magqc:::.runs(c(FALSE, FALSE))), 0)
})

# The spike detector's "nothing else is flagged" property is asserted over
# many noise draws, not one seed, so a regression cannot hide behind (or be
# faked by) a single lucky realisation (#14).
steep_curve <- function(seed, spike = 0, sd = 0.02, at = 200) {
  set.seed(seed)
  x <- 200 * sin(seq(0, 6, length.out = 400)) + rnorm(400, 0, sd)
  x[at] <- x[at] + spike
  x
}

test_that("spike detector finds the spike and only the spike across 25 noise draws", {
  hits <- lapply(1:25, function(s) magqc:::.spike_detect(steep_curve(s, spike = 5), k = 21, nsigma = 6))
  found <- vapply(hits, function(h) identical(which(h$outlier), 200L), logical(1))
  expect_true(all(found), info = paste("seeds missing or over-flagging:", paste(which(!found), collapse = " ")))
  amps <- vapply(hits, function(h) h$amplitude[200], numeric(1))
  expect_lt(max(abs(amps - 5)), 0.15)          # estimate noise is ~1.4 sigma = 0.03 nT
  cleaned <- vapply(hits, function(h) h$cleaned[200], numeric(1))
  expect_lt(max(abs(cleaned - 200 * sin(6 * 199 / 399))), 0.15)
})

test_that("spike detector has no false positives on 25 clean draws (10,000 samples)", {
  fp <- vapply(1:25, function(s) sum(magqc:::.spike_detect(steep_curve(s), k = 21, nsigma = 6)$outlier), integer(1))
  expect_equal(sum(fp), 0L)
})

test_that("spike detector sensitivity floor is where it was measured, not where naive sigma says", {
  # Naively a spike A gives a 4th difference of 6A against sqrt(70)*sd of
  # noise, so 0.25 nT at sd 0.02 would be ~9 sigma. In practice the running
  # MAD is dilated by a running max (so noisy-segment edges are not read as
  # spikes), which raises the floor: measured over 40 draws, 25x sd is always
  # found, 12x sd about 40% of the time, and <= 8x sd never. Pin the two ends.
  detected <- function(amp) vapply(1:20, function(s) 200L %in% which(magqc:::.spike_detect(steep_curve(s, spike = amp))$outlier), logical(1))
  expect_true(all(detected(0.50)))     # 25x noise sd
  expect_false(any(detected(0.15)))    #  8x noise sd
})

test_that("a missing value is filled before differencing and never flagged itself", {
  x <- steep_curve(1, spike = 5)
  x[100] <- NA
  h <- magqc:::.spike_detect(x, k = 21, nsigma = 6)
  expect_equal(which(h$outlier), 200L)
  expect_true(is.na(h$cleaned[100]))
})

test_that("mean bearing is taken on the circle", {
  ang_diff <- function(a, b) abs(((a - b + 180) %% 360) - 180)
  expect_lt(ang_diff(magqc:::.mean_bearing(c(359, 1)), 0), 1e-8)
  expect_lt(ang_diff(magqc:::.mean_bearing(c(170, 190)), 180), 1e-8)
  expect_equal(magqc:::.compass(magqc:::.mean_bearing(c(358.5, 0.5, 1.2))), "N")
})

test_that("tls line offsets are perpendicular distances", {
  x <- seq(0, 100, by = 1); y <- 2 * x + 5
  y[50] <- y[50] + sqrt(5) * 3  # 3 m perpendicular offset from a slope-2 line
  f <- magqc:::.tls_line(x, y)
  expect_equal(abs(f$offset[50]), 3, tolerance = 0.05)
  expect_lt(max(abs(f$offset[-50])), 0.05)
})
