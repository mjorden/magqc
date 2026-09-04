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
  # size budget (#15): with the full plotly bundle the demo report is < 6 MB,
  # with the partial bundle ~3 MB; either way well under the 8.6 MB it was
  expect_lt(file.size(out), 6e6)
})

test_that("diagnostic traces are one per line without per-point labels (#15)", {
  res <- run_qc(sim_survey(seed = 5))
  b <- plotly::plotly_build(plot_fourth_difference(res))$x$data
  lines <- vapply(b, function(t) t$mode == "lines" && length(t$y) > 0, logical(1))
  expect_equal(sum(lines), nrow(res$lines))
  # exactly one legend entry for the series, named for the series, not a line
  shown <- Filter(function(t) isTRUE(t$showlegend) && t$mode == "lines", b)
  expect_equal(vapply(shown, function(t) t$name, character(1)), "fourth difference")
  expect_true(all(vapply(b[lines], function(t) t$type == "scatter", logical(1))))
  expect_equal(sum(vapply(b[lines], function(t) !is.null(t$text), logical(1))), 1)   # only the legend trace
  # no NA gap rows: the traces hold about one point per finite sample
  # (plotly_build drops most of the NA fourth differences at segment ends;
  # a synthetic gap row per line would add nrow(res$lines) points)
  total <- sum(vapply(b[lines], function(t) length(t$y), integer(1)))
  expect_lte(total, nrow(res$survey))
  expect_lt(abs(total - sum(!is.na(res$fourth_difference))), nrow(res$lines))
  # times are strings (numeric epoch ms would render in the viewer's local zone)
  expect_true(all(vapply(b[lines], function(t) is.character(t$x) && grepl("^\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}\\.\\d$", t$x[1]), logical(1))))
  cl <- plotly::plotly_build(plot_clearance(res))$x$data
  n_pts <- sum(vapply(cl[vapply(cl, function(t) t$mode == "lines", logical(1))], function(t) length(t$y), integer(1)))
  expect_lt(n_pts, nrow(res$survey) / 2)                       # decimated
  expect_gt(n_pts, nrow(res$survey) / 4)
  withr::local_options(magqc.partial_bundle = FALSE)
  expect_identical(magqc:::.slim(plot_diurnal(res))$dependencies, plot_diurnal(res)$dependencies)
})
