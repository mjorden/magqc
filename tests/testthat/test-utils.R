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

test_that("spike detector finds a spike on a steep smooth curve and nothing else", {
  set.seed(1)
  x <- 200 * sin(seq(0, 6, length.out = 400)) + rnorm(400, 0, 0.02)
  x[200] <- x[200] + 5
  h <- magqc:::.spike_detect(x, k = 21, nsigma = 6)
  expect_equal(which(h$outlier), 200)
  expect_equal(h$amplitude[200], 5, tolerance = 0.05)
  expect_lt(abs(h$cleaned[200] - 200 * sin(6 * 199 / 399)), 0.1)
  # a series with a missing value is filled before differencing, never flagged there
  x[100] <- NA
  h2 <- magqc:::.spike_detect(x, k = 21, nsigma = 6)
  expect_equal(which(h2$outlier), 200)
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
