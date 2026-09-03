test_that("map and diagnostic widgets build", {
  res <- run_qc(sim_survey(seed = 5))
  expect_s3_class(qc_map(res), "leaflet")
  expect_s3_class(plot_fourth_difference(res), "plotly")
  expect_s3_class(plot_clearance(res), "plotly")
  expect_s3_class(plot_diurnal(res), "plotly")
  expect_s3_class(plot_crossovers(res), "plotly")

  # without lon/lat the map falls back to a projected plotly scatter
  d <- res$survey; d$lon <- NULL; d$lat <- NULL
  r2 <- run_qc(as_survey(d), spec = res$spec, base = res$base)
  expect_false(r2$stats$has_lonlat)
  expect_s3_class(qc_map(r2), "plotly")
})

test_that("qc_report renders a self-contained HTML file", {
  skip_if_not(rmarkdown::pandoc_available("2.0"), "pandoc not available")
  res <- run_qc(sim_survey(seed = 5))
  out <- withr::local_tempfile(fileext = ".html")
  qc_report(res, out, title = "Test block")
  expect_true(file.exists(out))
  html <- readLines(out, warn = FALSE)
  expect_true(any(grepl("Test block", html)))
  expect_true(any(grepl("Scorecard", html)))
  expect_true(any(grepl("leaflet", html)))
  expect_false(dir.exists(sub("\\.html$", "_files", out)))
})
