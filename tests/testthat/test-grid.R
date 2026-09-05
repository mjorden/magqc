sim <- sim_survey(seed = 4, defects = FALSE)
res <- run_qc(sim)

test_that("run_qc grids the residual field by default", {
  g <- res$grid
  expect_s3_class(g, "magqc_grid")
  expect_equal(g$channel, "mag_res")
  expect_equal(run_qc(sim, igrf = FALSE)$grid$channel, "mag_lev")
  expect_equal(g$cell, res$spec$line_spacing / 4)
  expect_equal(dim(g$z), c(length(g$x), length(g$y)))
  expect_true(all(diff(g$x) == g$cell) && all(diff(g$y) == g$cell))
  # covers the data extent to within a cell; the margin is inside the
  # blanking distance, so at most a few end-of-line nodes are blanked
  expect_lte(g$x[1], min(res$survey$x)); expect_gte(max(g$x), max(res$survey$x))
  expect_gt(g$x[1], min(res$survey$x) - g$cell); expect_lt(max(g$x), max(res$survey$x) + g$cell)
  expect_lt(g$n_blank, 0.01 * g$n_nodes)
  expect_equal(g$iterations, 0L)   # direct solve
  # a tight blanking distance leaves only a corridor around each line
  g2 <- grid_field(res, blank_distance = 40)   # one cell each side of a line
  expect_gt(g2$n_blank, 0.15 * g2$n_nodes)
  expect_true(is.na(g2$z[1, 1]))
  expect_lt(g$rms_fit, 2)          # nT, against a 2,000 nT anomaly sampled every 6 m
  expect_true(is.list(g$bounds) && g$bounds$north > g$bounds$south && g$bounds$east > g$bounds$west)
  expect_output(print(g), "minimum|mincurv")
  expect_null(run_qc(sim, grid = FALSE)$grid)
})

test_that("flagged spikes are left out of the grid", {
  rd <- run_qc(sim_survey(seed = 42))
  expect_equal(rd$grid$n_excluded, 25)
  v <- rd$survey[[rd$grid$channel]]
  expect_equal(rd$grid$n_data, sum(is.finite(v)) - 25)
  # the grid at a spike location is the smooth field, not the spike
  sp <- rd$flags[rd$flags$check == "spikes", ][1, ]
  g <- rd$grid
  at <- magqc:::.bilinear_at(g$z, (sp$x - g$x[1]) / g$cell + 1, (sp$y - g$y[1]) / g$cell + 1)
  expect_lt(abs(at - v[sp$i_start]), abs(sp$value))
  with_spikes <- grid_field(rd, exclude = integer(0))
  expect_equal(with_spikes$n_excluded, 0)
  expect_output(print(rd$grid), "25 flagged spikes excluded")
})

test_that("minimum curvature reproduces a plane and honours the data", {
  s <- as.data.frame(sim$data)
  s$mag <- 50000 + 0.01 * s$x - 0.02 * s$y
  s$mag_raw <- s$mag
  sv <- as_survey(s)
  g <- grid_field(sv, channel = "mag", cell = 50, blank_distance = 200)
  plane <- outer(g$x, g$y, function(x, y) 50000 + 0.01 * x - 0.02 * y)
  err <- g$z - plane
  # the local-plane binning recovers the plane to ~0.15 nT at worst (cells at
  # line ends holding one or two samples) against a 0.5 nT bias for the mean
  expect_lt(max(abs(err), na.rm = TRUE), 0.25)
  expect_lt(mean(abs(err), na.rm = TRUE), 0.05)
  expect_lt(g$rms_fit, 0.05)
  # minimum curvature overshoots a little between lines at a steep anomaly
  # (the smoothest surface rings); a few percent of the range is normal
  a <- res$grid
  rng <- range(res$survey[[a$channel]], na.rm = TRUE)
  expect_gte(min(a$z, na.rm = TRUE), rng[1] - 0.1 * diff(rng))
  expect_lte(max(a$z, na.rm = TRUE), rng[2] + 0.1 * diff(rng))
})

test_that("idw gridding and argument validation work", {
  g <- grid_field(res, method = "idw", cell = 100)
  expect_equal(g$method, "idw")
  expect_lt(g$rms_fit, 30)
  expect_error(grid_field(sim$data), "`cell` is required")
  expect_error(grid_field(res, channel = "nope"), "not in the survey")
  expect_error(grid_field(1), "magqc_survey")
  expect_error(grid_field(res, cell = 0), "positive")
})

test_that("helpers: bilinear, resample, dilate", {
  z <- matrix(1:6, 2, 3)
  expect_equal(magqc:::.bilinear_at(z, c(1, 2, 1.5), c(1, 3, 2)), c(1, 6, 3.5))
  expect_true(is.na(magqc:::.bilinear_at(z, 0.5, 1)))
  r <- magqc:::.resample(z, 3, 5)
  expect_equal(dim(r), c(3, 5)); expect_equal(r[1, 1], 1); expect_equal(r[3, 5], 6)
  m <- matrix(FALSE, 5, 5); m[3, 3] <- TRUE
  expect_equal(sum(magqc:::.dilate(m, 1)), 9)
  expect_equal(sum(magqc:::.dilate(m, 2)), 25)
})

test_that("the grid reaches the map, the projected fallback and the panel", {
  m <- qc_map(res)
  expect_s3_class(m, "leaflet")
  ctl <- Filter(function(cl) cl$method == "addLayersControl", m$x$calls)[[1]]
  expect_true("Gridded field" %in% unlist(ctl$args[[2]]))
  expect_true(any(grepl("imageOverlay", unlist(lapply(m$jsHooks$render, `[[`, "code")))))
  legends <- Filter(function(cl) cl$method == "addLegend", m$x$calls)
  expect_true(any(vapply(legends, function(cl) identical(cl$args[[1]]$group, "Gridded field"), logical(1))))

  s <- res$survey; s$lon <- NA_real_; s$lat <- NA_real_
  r2 <- run_qc(s, spec = res$spec, base = res$base)
  expect_null(r2$grid$bounds)
  pm <- plotly::plotly_build(qc_map(r2))
  expect_true("heatmap" %in% vapply(pm$x$data, function(t) t$type, character(1)))

  p <- plotly::plotly_build(plot_grid(res))
  expect_equal(p$x$data[[1]]$type, "heatmap")
  expect_null(plot_grid(run_qc(sim, grid = FALSE)))
})
