sim <- sim_survey(seed = 42)

test_that("levelling removes the simulated heading bias", {
  res <- run_qc(sim)
  lev <- res$levelling
  expect_s3_class(lev, "magqc_levelling")
  expect_lt(lev$rms[["after"]], 0.5)
  expect_lt(lev$rms[["after"]], lev$rms[["before"]] / 3)
  expect_true("mag_lev" %in% names(res$survey))
  # ties are untouched, traverses shifted by their coefficients
  tie <- res$survey$line_type == "tie"
  expect_equal(res$survey$mag_lev[tie], res$survey$mag[tie])
  expect_false(isTRUE(all.equal(res$survey$mag_lev[!tie], res$survey$mag[!tie])))
  # southbound lines carried the +2.5 nT bias, so their offsets are ~2.5 nT
  # above the northbound ones
  co <- lev$coefficients
  south <- as.integer(sub("L", "", co$line)) %% 20 == 10
  expect_equal(mean(co$c0[south]) - mean(co$c0[!south]), 2.5, tolerance = 0.3)
  # the dropout line only has two crossings, so it gets a constant
  expect_equal(co$order[co$line == "L1070"], "constant")
  expect_true(all(co$order[co$line != "L1070"] == "linear"))
  # heading error is still reported on pre-levelling misfits
  expect_equal(res$scorecard$status[res$scorecard$check == "heading_error"], "fail")
  expect_equal(res$scorecard$metric[res$scorecard$check == "levelling_residual"],
               unname(lev$rms[["after"]]))
  expect_output(print(lev), "before")
})

test_that("a linear correction recovers an injected per-line drift; a constant does not", {
  clean <- sim_survey(seed = 5, defects = FALSE)
  s <- clean$data
  idx <- which(s$line == "L1040")
  along <- cumsum(c(0, sqrt(diff(s$x[idx])^2 + diff(s$y[idx])^2)))
  s$mag_raw[idx] <- s$mag_raw[idx] + 0.002 * along           # 14 nT over the line
  lin <- run_qc(s, spec = clean$spec, base = clean$base, levelling_order = "linear")
  con <- run_qc(s, spec = clean$spec, base = clean$base, levelling_order = "constant")
  after <- function(r) { xo <- r$levelling$crossovers; max(abs(xo$misfit_after[xo$traverse == "L1040"])) }
  expect_lt(after(lin), 0.5)
  expect_gt(after(con), 3)   # ties at 1, 3 and 5 km: a constant leaves +/-4 nT at the outer two
  expect_equal(lin$levelling$coefficients$c1[lin$levelling$coefficients$line == "L1040"], 0.002, tolerance = 0.1)
  # the drift line fails the residual check under constant levelling only
  expect_equal(con$scorecard$status[con$scorecard$check == "levelling_residual"], "fail")
  expect_equal(lin$scorecard$status[lin$scorecard$check == "levelling_residual"], "pass")
  expect_true(any(grepl("L1040", con$flags$line[con$flags$check == "levelling_residual"])))
})

test_that("a clean survey needs no correction", {
  res <- run_qc(sim_survey(seed = 7, defects = FALSE))
  co <- res$levelling$coefficients
  expect_lt(max(abs(co$c0)), 0.3)
  expect_lt(max(abs(co$c1)) * 7000, 0.5)
  expect_lt(res$levelling$rms[["after"]], 0.3)
  expect_equal(res$scorecard$status[res$scorecard$check == "levelling_residual"], "pass")
})

test_that("levelling is skipped cleanly without ties and can be turned off", {
  s <- sim$data[sim$data$line_type == "traverse", ]
  s <- as_survey(as.data.frame(s))
  res <- run_qc(s, spec = sim$spec, base = sim$base)
  expect_null(res$levelling)
  expect_equal(res$scorecard$status[res$scorecard$check == "levelling_residual"], "n/a")
  expect_error(level_ties(s), "No traverse/tie crossovers")
  expect_error(level_ties(1), "magqc_survey")
  off <- run_qc(sim, levelling = FALSE)
  expect_null(off$levelling)
  expect_false("mag_lev" %in% names(off$survey))
  lev <- level_ties(off)
  expect_equal(lev$rms[["before"]], off$scorecard$metric[off$scorecard$check == "crossovers"])
})

test_that("crossover and profile plots show the levelled state", {
  res <- run_qc(sim)
  trace_names <- function(p) vapply(plotly::plotly_build(p)$x$data,
                                    function(t) if (is.null(t$name)) "" else as.character(t$name), character(1))
  expect_true(any(grepl("after levelling", trace_names(plot_crossovers(res)))))
  expect_true("levelled" %in% trace_names(plot_line_profile(res, "L1010")))
})
