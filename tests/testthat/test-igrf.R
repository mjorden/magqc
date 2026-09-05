sim <- sim_survey(seed = 6, defects = FALSE)

test_that("igrf_field gives a plausible main field for central Wyoming in 2024", {
  f <- igrf_field(-108.5, 43.2, as.POSIXct("2024-06-14 08:00", tz = "UTC"), altitude = 1600)
  expect_equal(nrow(f), 1)
  expect_gt(f$F, 50000); expect_lt(f$F, 55000)
  expect_gt(f$I, 60); expect_lt(f$I, 75)        # steep, positive down
  expect_gt(f$D, 5); expect_lt(f$D, 12)          # east
  expect_equal(attr(f, "epoch"), 2024.45, tolerance = 1e-3)
  expect_match(attr(f, "model"), "^igrf ")
  # the lattice interpolation reproduces a direct evaluation to a fraction of a nT
  s <- sim$data
  k <- c(1, 5000, 12000, nrow(s))
  fl <- igrf_field(s$lon, s$lat, s$time, altitude = 60)
  for (i in k) {
    d <- igrf::igrf(field = "main", year = attr(fl, "epoch"), type = "spheroid", altitude = 0.06,
                    latitude = s$lat[i], longitude = s$lon[i])
    expect_lt(abs(fl$F[i] - d$F), 0.1)
    expect_lt(abs(fl$I[i] - d$I), 0.01)
  }
  expect_error(igrf_field(NA, NA, Sys.time()), "No finite")
})

test_that("run_qc removes the IGRF and grids the residual", {
  res <- run_qc(sim)
  expect_true(all(c("mag_igrf", "mag_res") %in% names(res$survey)))
  expect_equal(res$survey$mag_res, res$survey$mag_lev - res$survey$mag_igrf)
  expect_equal(res$igrf$from, "mag_lev")
  expect_gt(res$igrf$F, 50000)
  expect_equal(res$igrf$altitude, mean(res$survey$gps_alt, na.rm = TRUE))   # GPS height, not clearance
  expect_equal(res$grid$channel, "mag_res")
  # the residual grid is the levelled-field grid minus a smooth surface of a
  # few nT per km (the main field's gradient across the block)
  g_lev <- grid_field(res, channel = "mag_lev")
  d <- res$grid$z - g_lev$z
  expect_lt(diff(range(d, na.rm = TRUE)), 60)
  expect_lt(abs(mean(d, na.rm = TRUE) + res$igrf$F), 40)
  # the map legend and the panel say what they show
  expect_match(magqc:::.channel_label("mag_res"), "residual")
  legends <- Filter(function(cl) cl$method == "addLegend", qc_map(res)$x$calls)
  expect_true(any(grepl("residual", unlist(lapply(legends, function(cl) as.character(unlist(cl$args)))))))
})

test_that("no longitude/latitude: IGRF is skipped and the levelled field is gridded", {
  s <- sim$data; s$lon <- NA_real_; s$lat <- NA_real_
  res <- run_qc(s, spec = sim$spec, base = sim$base)
  expect_null(res$igrf)
  expect_false("mag_res" %in% names(res$survey))
  expect_equal(res$grid$channel, "mag_lev")
  off <- run_qc(sim, igrf = FALSE)
  expect_null(off$igrf)
  expect_equal(off$grid$channel, "mag_lev")
})
