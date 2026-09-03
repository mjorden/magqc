# The simulated survey records every defect it injects, with the time window
# each one occupies, so each check can be held to two standards: it must flag
# the interval it was designed to catch, and it must not flag lines that have
# no such defect.

sim <- sim_survey(seed = 42)
res <- run_qc(sim)
truth <- sim$truth
flags <- res$flags

overlaps <- function(f, tr) {
  any(f$line == tr$line & f$time_start <= tr$time_end & f$time_end >= tr$time_start, na.rm = TRUE)
}
truth_for <- function(chk) truth[truth$check == chk, ]
flags_for <- function(chk) flags[flags$check == chk, ]

test_that("a clean survey passes every check", {
  clean <- run_qc(sim_survey(seed = 7, defects = FALSE))
  expect_equal(nrow(clean$flags), 0)
  expect_true(all(clean$scorecard$status == "pass"))
  expect_equal(clean$stats$accepted_pct, 100)
})

test_that("data gaps are caught only where injected", {
  tr <- truth_for("sample_interval"); f <- flags_for("sample_interval")
  expect_true(overlaps(f, tr))
  expect_setequal(unique(f$line), tr$line)
  expect_gt(f$value[1], 10)
})

test_that("sparse sampling is caught as spacing, not as a gap", {
  tr <- truth_for("along_line_spacing"); f <- flags_for("along_line_spacing")
  expect_true(overlaps(f, tr))
  expect_setequal(unique(f$line), tr$line)
  expect_false(tr$line %in% flags_for("sample_interval")$line)
})

test_that("departure from line is caught on the excursion line only", {
  tr <- truth_for("line_deviation"); f <- flags_for("line_deviation")
  expect_true(overlaps(f, tr))
  expect_setequal(unique(f$line), tr$line)
  expect_gt(max(f$value), 25)
})

test_that("a mispositioned line shows up in both adjacent separations", {
  f <- flags_for("line_separation")
  expect_equal(nrow(f), 2)
  expect_true(all(grepl("L1110", f$line)))
  expect_equal(sort(round(f$value, -1)), c(140, 260))
})

test_that("clearance excursion is caught only on its line", {
  tr <- truth_for("clearance"); f <- flags_for("clearance")
  expect_true(overlaps(f, tr))
  expect_setequal(unique(f$line), tr$line)
  expect_true(all(f$value > 0))
})

test_that("every injected spike is flagged and nothing else is", {
  tr <- truth_for("spikes"); f <- flags_for("spikes")
  hit <- vapply(seq_len(nrow(tr)), function(k) overlaps(f, tr[k, ]), logical(1))
  expect_true(all(hit))
  expect_equal(nrow(f), nrow(tr))
})

test_that("fourth-difference noise is caught on the noisy segment only", {
  tr <- truth_for("fourth_difference"); f <- flags_for("fourth_difference")
  expect_true(overlaps(f, tr))
  expect_setequal(unique(f$line), tr$line)
  expect_true(all(f$time_start >= tr$time_start - 1 & f$time_end <= tr$time_end + 1))
})

test_that("diurnal storm is caught and mapped onto the samples flown through it", {
  tr <- truth_for("diurnal"); f <- flags_for("diurnal")
  expect_gt(nrow(f), 0)
  expect_true(all(f$time_start >= tr$time_start - 120 & f$time_end <= tr$time_end + 120))
  expect_true(any(!is.na(f$i_start)))
  expect_true(all(res$diurnal$chord_dev[!res$diurnal$exceed] <= res$spec$max_diurnal_dev, na.rm = TRUE))
})

test_that("heading bias is detected from crossovers", {
  f <- flags_for("heading_error")
  expect_equal(nrow(f), 1)
  expect_equal(f$value, 2.5, tolerance = 0.25)
  # 14 traverses x 3 ties, less the one crossing that falls inside L1070's dropout
  expect_equal(nrow(res$crossovers), 14 * 3 - 1)
  missing <- setdiff(paste("L1070", c("T9000", "T9010", "T9020")),
                     paste(res$crossovers$traverse, res$crossovers$tie))
  expect_length(missing, 1)
  expect_lt(res$scorecard$metric[res$scorecard$check == "crossovers"], res$spec$max_crossover_rms)
})

test_that("scorecard and stats are consistent with the flags", {
  sc <- res$scorecard
  expect_equal(sum(sc$n_flags), nrow(flags))
  expect_true(all(sc$status[sc$n_flags > 0] == "fail"))
  expect_equal(res$stats$checks_failed, sum(sc$status == "fail"))
  expect_lt(res$stats$accepted_pct, 100)
  expect_equal(sum(res$lines$length_km), res$stats$line_km)
  expect_output(print(res), "checks passed")
})

test_that("checks accept a custom spec and a flight plan", {
  # with despiking disabled the 5-40 nT spikes reach the fourth difference, so
  # its limit has to clear 6 * 40 / 16 = 15 nT too
  loose <- survey_spec(max_deviation = 80, clearance_tol = 50, fourth_diff_tol = 20,
                       max_heading_error = 5, line_spacing_tol = 0.5, max_gap = 30,
                       max_along_line_sep = 1000, spike_nsigma = 1e6, max_diurnal_dev = 100)
  r2 <- run_qc(sim, spec = loose)
  expect_equal(nrow(r2$flags), 0)

  plan <- data.frame(line = "L1090", x0 = 9 * 200, y0 = 0, x1 = 9 * 200, y1 = 7000)
  f <- check_line_deviation(res$survey, res$spec, plan = plan)
  expect_true(any(grepl("planned", f$description)))
  expect_gt(max(f$value), 40)
})
