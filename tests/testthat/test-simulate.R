# The simulator must stay truthful across its parameter space: every truth
# row names a line that exists and a window that exists, and each injected
# defect is caught on that line (#3, #4).

truth_is_consistent <- function(sim) {
  tr <- sim$truth
  d <- sim$data
  located <- !is.na(tr$line)
  expect_true(all(tr$line[located] %in% d$line))
  expect_false(anyNA(tr$time_start[located]))
  expect_false(anyNA(tr$time_end[located]))
  expect_true(all(tr$time_end >= tr$time_start, na.rm = TRUE))
  # the storm window has no line but does have a time window
  expect_false(anyNA(tr$time_start[tr$check == "diurnal"]))
}

each_defect_is_caught <- function(sim) {
  res <- run_qc(sim)
  tr <- sim$truth[!is.na(sim$truth$line) & sim$truth$check != "spikes", ]
  for (k in seq_len(nrow(tr))) {
    f <- res$flags[res$flags$check == tr$check[k], ]
    hit <- if (tr$check[k] == "line_separation") any(grepl(tr$line[k], f$line, fixed = TRUE)) else
      any(f$line == tr$line[k] & f$time_start <= tr$time_end[k] + 1 & f$time_end >= tr$time_start[k] - 1, na.rm = TRUE)
    expect_true(hit, info = sprintf("%s on %s not caught", tr$defect[k], tr$line[k]))
  }
  expect_gt(nrow(res$flags[res$flags$check == "spikes", ]), 0)
  expect_equal(nrow(res$flags[res$flags$check == "heading_error", ]), 1)
  res
}

test_that("default geometry keeps the documented defect targets", {
  sim <- sim_survey(seed = 42)
  tr <- sim$truth
  target <- function(chk) unique(tr$line[tr$check == chk])
  expect_equal(target("fourth_difference"), "L1050")
  expect_equal(target("line_deviation"), "L1090")
  expect_equal(target("line_separation"), "L1110")
  expect_equal(target("clearance"), "L1030")
  expect_equal(target("sample_interval"), "L1070")
  expect_equal(target("along_line_spacing"), "L1010")
  expect_equal(sum(tr$check == "spikes"), 25)
  truth_is_consistent(sim)
})

test_that("small blocks and short lines simulate and are caught (#3, #4)", {
  # (3 lines x 3 km would leave a single tie whose only southbound crossing
  # sits inside the dropout, making heading error legitimately n/a; use 5 km)
  for (args in list(list(n_lines = 5), list(n_lines = 3, line_length = 5000),
                    list(line_length = 4000, sample_rate = 5), list(n_lines = 20, line_spacing = 100))) {
    expect_no_warning(sim <- do.call(sim_survey, c(args, list(seed = 3))))
    expect_s3_class(sim$data, "magqc_survey")
    truth_is_consistent(sim)
    n_expected <- if (is.null(args$n_lines)) 14 else args$n_lines
    expect_equal(length(unique(sim$data$line[sim$data$line_type == "traverse"])), n_expected)
    each_defect_is_caught(sim)
  }
})

test_that("defects that cannot be placed are skipped with a warning and no truth row", {
  w <- capture_warnings(sim <- sim_survey(n_lines = 2, line_length = 150, seed = 1))
  expect_true(any(grepl("skipping defect 'data dropout'", w)))
  expect_true(any(grepl("skipping defect 'noisy segment'", w)))
  truth_is_consistent(sim)
  expect_false("sample_interval" %in% sim$truth$check)
  expect_true(all(sim$truth$check %in% c("heading_error", "diurnal", "spikes", "line_separation")))
  expect_error(sim_survey(n_lines = 1), "at least 2")
  expect_error(sim_survey(line_length = 0), "positive")
})

test_that(".pick_lines keeps defects on distinct lines when it can", {
  ln <- sprintf("L%d", 1:14)
  p <- magqc:::.pick_lines(ln, c(a = 0.43, b = 0.71, c = 0.86, d = 0.29, e = 0.57, f = 0.14))
  expect_equal(unname(p), c("L6", "L10", "L12", "L4", "L8", "L2"))
  p3 <- magqc:::.pick_lines(ln[1:3], c(a = 0.43, b = 0.71, c = 0.86, d = 0.29))
  expect_equal(length(unique(p3[1:3])), 3)
  expect_true(p3[["d"]] %in% ln[1:3])
})
