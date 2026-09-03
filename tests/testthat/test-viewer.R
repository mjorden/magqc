res <- run_qc(sim_survey(seed = 11))

test_that("line profile, flag table and registry build", {
  p <- plot_line_profile(res, "L1050")
  expect_s3_class(p, "plotly")
  expect_error(plot_line_profile(res, "nope"), "not in the survey")
  ft <- flag_table(res)
  expect_equal(nrow(ft), nrow(res$flags))
  expect_named(ft, c("check", "line", "start", "end", "fid", "value", "limit", "units", "finding"))
  expect_true(all(ft$check %in% check_registry()$label))
  clean <- run_qc(sim_survey(seed = 11, defects = FALSE))
  expect_equal(nrow(flag_table(clean)), 0)
})

test_that("map polylines carry escaped labels and layer ids", {
  s <- res$survey
  s$line[s$line == "L1000"] <- "L<b>1000"
  r2 <- run_qc(s, spec = res$spec, base = res$base)
  calls <- qc_map(r2)$x$calls
  polylines <- Filter(function(cl) cl$method == "addPolylines", calls)
  strings <- unlist(lapply(polylines, function(cl) rapply(cl$args, as.character, how = "unlist")))
  # leaflet escapes the label itself; the raw name must survive only as layerId,
  # and nothing may be double-escaped
  expect_equal(sum(strings == "L&lt;b&gt;1000"), 1)
  expect_equal(sum(strings == "L<b>1000"), 1)
  expect_false(any(grepl("&amp;lt;", strings, fixed = TRUE)))
  # popups are not escaped by leaflet, so .popup() must escape line names itself
  all_strings <- unlist(lapply(calls, function(cl) rapply(cl$args, as.character, how = "unlist")))
  popups <- all_strings[startsWith(all_strings, "<b>")]
  expect_gt(length(popups), 0)
  expect_false(any(grepl("L<b>1000", popups, fixed = TRUE)))
})

test_that("the viewer app parses and serves a result", {
  skip_if_not_installed("shiny"); skip_if_not_installed("bslib")
  app_dir <- system.file("shiny", package = "magqc")
  skip_if(!nzchar(app_dir), "app not installed")
  expect_silent(parse(file.path(app_dir, "app.R")))
  expect_error(run_viewer(result = 1), "magqc_result")

  shiny::shinyOptions(magqc_result = res)
  app <- shiny::shinyAppDir(app_dir)
  shiny::testServer(app, {
    session$setInputs(run = 0, source = "sim", seed = 11, n_lines = 14, defects = TRUE, line = "L1050")
    expect_s3_class(res(), "magqc_result")
    expect_equal(res()$stats$n_samples, !!res$stats$n_samples)
    expect_true(grepl("%", output$vb_pct))
    expect_type(output$profile, "character")
    expect_type(output$scorecard, "character")
  })
  shiny::shinyOptions(magqc_result = NULL)
})
