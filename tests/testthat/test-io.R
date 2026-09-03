test_that("as_survey validates and infers", {
  df <- data.frame(line = c("L10", "L10", "T90"), time = as.POSIXct(1:3, origin = "1970-01-01", tz = "UTC"),
                   x = 1:3, y = 1:3, mag_raw = 1:3)
  s <- as_survey(df, origin = c(lon = -100, lat = 40))
  expect_s3_class(s, "magqc_survey")
  expect_equal(s$line_type, c("traverse", "traverse", "tie"))
  expect_true(all(is.finite(s$lon)))
  expect_equal(s$mag, s$mag_raw)
  expect_error(as_survey(df[, -1]), "line")
  expect_error(as_survey(transform(df, mag_raw = NULL)), "mag_raw")
  expect_error(as_survey(transform(df, time = 1:3)), "POSIXct")
})

test_that("read_xyz parses Geosoft-style files", {
  f <- withr::local_tempfile(fileext = ".xyz")
  writeLines(c(
    "/ Survey block A",
    "/ X Y TIME MAG RALT",
    "Line 1000",
    "  100 200 3600.0 54000.1 60",
    "  100 206 3600.1 54000.2 61",
    "  100 212 3600.2 *       59",
    "Tie 9000",
    "  0 500 7200.0 54010.0 60",
    "  6 500 7200.1 54010.5 60"), f)
  s <- read_xyz(f, col_map = c(x = "X", y = "Y", mag_raw = "MAG", radar_alt = "RALT", time = "TIME"),
                time_origin = as.POSIXct("2024-01-01", tz = "UTC"))
  expect_equal(nrow(s), 5)
  expect_equal(unique(s$line), c("1000", "9000"))
  expect_equal(s$line_type, c(rep("traverse", 3), rep("tie", 2)))
  expect_true(is.na(s$mag_raw[3]))
  expect_equal(format(s$time[1], "%H:%M:%S"), "01:00:00")
  expect_error(read_xyz(f, col_map = c(x = "NOPE")), "not in file")
})
