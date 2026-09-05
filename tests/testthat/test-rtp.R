# A grid of a point dipole's anomaly, built directly (no survey), for the
# physics checks
dipole_grid <- function(inclination, declination, n = 96, cell = 25, h = 200, x0 = 0, y0 = 0) {
  gx <- (seq_len(n) - n / 2) * cell; gy <- (seq_len(n) - n / 2) * cell
  z <- outer(gx, gy, function(x, y) magqc:::.dipole_anomaly(x, y, x0, y0, h, inclination, declination))
  structure(list(x = gx, y = gy, z = z, cell = cell, channel = "mag_res", method = "mincurv",
                 blank_distance = Inf, n_data = n * n, n_excluded = 0L, n_nodes = n * n, n_blank = 0L,
                 rms_fit = 0, iterations = 0L, bounds = NULL), class = "magqc_grid")
}
interior <- function(z, frac = 0.6) {
  n <- nrow(z); m <- ncol(z); i <- round(n * (1 - frac) / 2):round(n * (1 + frac) / 2)
  j <- round(m * (1 - frac) / 2):round(m * (1 + frac) / 2)
  z[i, j]
}
peak_at <- function(g) { k <- which(g$z == max(g$z, na.rm = TRUE), arr.ind = TRUE)[1, ]; c(g$x[k[1]], g$y[k[2]]) }

test_that("RTP of an inclined-field dipole anomaly is the pole anomaly of the same dipole", {
  obs <- dipole_grid(60, 20)
  pole <- dipole_grid(90, 0)
  r <- reduce_to_pole(obs, 60, 20)
  expect_s3_class(r, "magqc_grid"); expect_equal(r$channel, "rtp")
  a <- interior(r$z); b <- interior(pole$z)
  expect_gt(stats::cor(as.vector(a), as.vector(b)), 0.995)
  expect_lt(sqrt(mean((a - b)^2)) / diff(range(b)), 0.03)
  # the observed anomaly peaks off the source; the reduced one sits on it
  expect_gt(sqrt(sum(peak_at(obs)^2)), 2 * obs$cell)
  expect_lte(sqrt(sum(peak_at(r)^2)), obs$cell * 1.5)
  # ...and is symmetric about it (the source sits at index n/2, so mirror
  # about that node rather than the matrix centre)
  zc <- r$z; c0 <- nrow(zc) / 2; k <- -25:25
  fwd <- zc[c0 + k, c0 + k]; bwd <- zc[c0 - k, c0 - k]
  expect_lt(max(abs(fwd - bwd)) / diff(range(zc)), 0.03)
})

test_that("RTP is the identity in a vertical field and bounded at low latitude", {
  pole <- dipole_grid(90, 0)
  r <- reduce_to_pole(pole, 90, 0)
  expect_lt(max(abs(interior(r$z) - interior(pole$z))) / diff(range(pole$z)), 0.02)
  # the declination is irrelevant at the pole
  r2 <- reduce_to_pole(pole, 90, 137)
  expect_equal(r2$z, r$z, tolerance = 1e-8)
  # 5 degrees inclination: the raw operator would amplify some wavenumbers
  # ~100x; the amplitude inclination keeps the result the size of the anomaly
  low <- dipole_grid(5, 0)
  rl <- reduce_to_pole(low, 5, 0)
  expect_equal(attr(rl, "amp_inclination"), 20)
  expect_lt(diff(range(rl$z)), 5 * diff(range(dipole_grid(90, 0)$z)))
  expect_true(all(is.finite(rl$z)))
})

test_that("blanks are infilled for the transform and restored; inputs are validated", {
  g <- dipole_grid(60, 20)
  g$z[1:5, ] <- NA; g$z[40:44, 60:64] <- NA
  r <- reduce_to_pole(g, 60, 20)
  expect_true(all(is.na(r$z[1:5, ])) && all(is.na(r$z[40:44, 60:64])))
  expect_true(all(is.finite(r$z[!is.na(g$z)])))
  full <- reduce_to_pole(dipole_grid(60, 20), 60, 20)
  expect_lt(max(abs(r$z[50:80, 20:50] - full$z[50:80, 20:50])) / diff(range(full$z)), 0.05)
  expect_error(reduce_to_pole(1, 60, 20), "magqc_grid")
  expect_error(reduce_to_pole(g, NA, 20), "finite")
  expect_output(print(r), "reduced to the pole")
  expect_equal(magqc:::.fftfreq(8, 1), c(0, 1, 2, 3, -4, -3, -2, -1) / 8)
  expect_equal(magqc:::.fftfreq(7, 2), c(0, 1, 2, 3, -3, -2, -1) / 14)
})

test_that("run_qc reduces the residual grid to the pole and the map and panel carry it", {
  res <- run_qc(sim_survey(seed = 4, defects = FALSE))
  expect_s3_class(res$rtp, "magqc_grid")
  expect_equal(attr(res$rtp, "inclination"), res$igrf$I)
  expect_equal(res$rtp$source_channel, "mag_res")
  expect_equal(dim(res$rtp$z), dim(res$grid$z))
  expect_equal(is.na(res$rtp$z), is.na(res$grid$z))
  m <- qc_map(res)
  ctl <- Filter(function(cl) cl$method == "addLayersControl", m$x$calls)[[1]]
  expect_true(all(c("Gridded field", "Reduced to pole") %in% unlist(ctl$args[[2]])))
  expect_true(any(vapply(m$x$calls, function(cl) cl$method == "hideGroup", logical(1))))
  expect_equal(length(m$jsHooks$render), 2)
  p <- plotly::plotly_build(plot_grid(res, "rtp"))
  expect_equal(p$x$data[[1]]$type, "heatmap")
  expect_null(run_qc(sim_survey(seed = 4), rtp = FALSE)$rtp)
  expect_null(run_qc(sim_survey(seed = 4), igrf = FALSE)$rtp)
})
