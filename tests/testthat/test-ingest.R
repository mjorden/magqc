nmea <- function(body) sprintf("$%s*%s", body, magqc:::.nmea_xor(body))

test_that("read_nmea parses GGA/RMC/ZDA, checks checksums and handles dates", {
  f <- withr::local_tempfile(fileext = ".nmea")
  lines <- c(
    nmea("GPRMC,235958.00,A,4312.0000,N,10830.0000,W,116.6,359.0,140624,,,A"),
    nmea("GPGGA,235958.00,4312.0000,N,10830.0000,W,1,09,0.9,1580.0,M,-17.0,M,,"),
    nmea("GPGGA,235959.00,4312.0300,N,10830.0000,W,1,09,0.9,1581.0,M,-17.0,M,,"),
    "$GPGGA,235959.50,4312.0450,N,10830.0000,W,1,09,0.9,1581.5,M,-17.0,M,,*00",   # bad checksum
    nmea("GNGGA,000000.00,4312.0600,N,10830.0000,W,1,10,0.8,1582.0,M,-17.0,M,,"),  # midnight passed
    nmea("GPGSV,3,1,11,03,03,111,00,04,15,270,00"),                                 # ignored
    nmea("GPGGA,000001.00,4312.0900,N,10830.0000,W,1,10,0.8,1583.0,M,-17.0,M,,"))
  writeLines(lines, f)
  g <- read_nmea(f)
  expect_equal(nrow(g), 4)
  expect_equal(attr(g, "n_bad_checksum"), 1)
  expect_equal(attr(g, "n_sentences"), 7)
  expect_equal(format(g$time[1], "%Y-%m-%d %H:%M:%S"), "2024-06-14 23:59:58")
  expect_equal(format(g$time[3], "%Y-%m-%d %H:%M:%S"), "2024-06-15 00:00:00")   # rolled over
  expect_equal(g$lat[1], 43.2); expect_equal(g$lon[1], -108.5)
  expect_equal(g$alt, c(1580, 1581, 1582, 1583))
  expect_equal(g$n_sat, c(9L, 9L, 10L, 10L))
  expect_equal(g$speed[1], 116.6 * 0.514444, tolerance = 1e-6)   # RMC within a second
  expect_equal(g$course[1], 359)
  expect_true(is.na(g$speed[3]))                                   # no RMC near the later fixes

  # no date sentences: needs `date`
  writeLines(lines[c(2, 3)], f)
  expect_error(read_nmea(f), "pass `date`")
  g2 <- read_nmea(f, date = "2024-06-14")
  expect_equal(format(g2$time[1], "%Y-%m-%d"), "2024-06-14")
  # ZDA supplies the date too
  writeLines(c(nmea("GPZDA,235958.00,14,06,2024,00,00"), lines[c(2, 3)]), f)
  expect_equal(format(read_nmea(f)$time[1], "%Y-%m-%d"), "2024-06-14")
  writeLines("$GPGSV,3,1,11,03,03,111,00", f)
  expect_error(read_nmea(f), "No GGA")
})

test_that("read_mag_log handles the presets and each time style", {
  f <- withr::local_tempfile(fileext = ".mag")
  writeLines(c("# header", "$ 52341.123,1234,08:00:00.100", "$ 52341.456,1230,08:00:00.200",
               "junk line", "$ 52341.789,1231,08:00:00.300"), f)
  m <- read_mag_log(f, "geometrics", date = "2024-06-14")
  expect_equal(nrow(m), 3)
  expect_equal(m$mag_raw, c(52341.123, 52341.456, 52341.789))
  expect_equal(m$signal, c(1234, 1230, 1231))
  expect_equal(as.numeric(diff(m$time)), c(0.1, 0.1), tolerance = 1e-4)
  secs_since <- function(t, day) as.numeric(t) - as.numeric(as.POSIXct(day, tz = "UTC"))
  expect_equal(secs_since(m$time[1], "2024-06-14"), 8 * 3600 + 0.1, tolerance = 1e-6)
  expect_error(read_mag_log(f, "geometrics"), "pass `date`")
  expect_output(print(mag_format("geometrics")), "hms")

  writeLines(c("52341.12 99 235959.9", "52341.34 99 000000.0", "52341.56 98 000000.1"), f)
  g <- read_mag_log(f, mag_format("gem"), date = "2024-06-14")
  expect_equal(format(g$time[2], "%Y-%m-%d %H:%M:%S"), "2024-06-15 00:00:00")   # midnight

  writeLines(c("t,field", "0.0,52000.1", "0.1,52000.2"), f)
  s <- read_mag_log(f, mag_format(columns = c(time = 1, mag = 2), time = "seconds", skip = 1), date = "2024-06-14")
  expect_equal(as.numeric(s$time - as.POSIXct("2024-06-14", tz = "UTC"), units = "secs"), c(0, 0.1), tolerance = 1e-6)
  e <- read_mag_log(f, mag_format(columns = c(time = 1, mag = 2), time = "elapsed", skip = 1),
                    date = as.POSIXct("2024-06-14 08:00:00", tz = "UTC"))
  expect_equal(secs_since(e$time[2], "2024-06-14"), 8 * 3600 + 0.1, tolerance = 1e-6)
  i <- read_mag_log(f, mag_format(columns = c(mag = 2), time = "index", rate = 10, skip = 1),
                    date = as.POSIXct("2024-06-14 08:00:00", tz = "UTC"))
  expect_equal(secs_since(i$time[2], "2024-06-14"), 8 * 3600 + 0.1, tolerance = 1e-6)
  writeLines(c("2024-06-14 08:00:00.0 52000.1", "2024-06-14 08:00:00.1 52000.2"), f)
  dtf <- read_mag_log(f, mag_format(sep = " (?=[0-9]+\\.[0-9]+$)", columns = c(time = 1, mag = 2),
                                    time = "datetime", time_format = "%Y-%m-%d %H:%M:%OS"))
  expect_equal(nrow(dtf), 2)
  expect_error(mag_format(columns = c(time = 1)), "`mag`")
  expect_error(mag_format(columns = c(mag = 1)), "`time`")
  expect_error(mag_format(columns = c(mag = 1), time = "index"), "`rate`")
})

test_that("merge_streams interpolates positions, honours lag and drops outages", {
  t0 <- as.POSIXct("2024-06-14 08:00:00", tz = "UTC")
  gps <- tibble::tibble(time = t0 + seq(0, 20, by = 1), lat = 43.2 + seq(0, 20, by = 1) * 1e-4,
                        lon = -108.5, alt = 1600)
  mag <- tibble::tibble(time = t0 + seq(0, 20, by = 0.1), mag_raw = 52000)
  m <- merge_streams(mag, gps)
  expect_equal(nrow(m), 201)
  expect_equal(m$lat[6], 43.2 + 0.5e-4)                      # halfway between fixes
  expect_equal(attr(m, "n_dropped"), 0)
  expect_true(all(c("fid", "x", "y", "gps_alt") %in% names(m)))
  # a positive lag: the magnetometer stamps are later, so the sample at
  # mag time 1.0 s belongs at GPS time 0.5 s
  ml <- merge_streams(mag, gps, lag = 0.5)
  expect_equal(ml$lat[ml$time == t0 + 1], 43.2 + 0.5e-4)
  # an outage: fixes missing for 5 s, samples inside it dropped
  gps2 <- gps[!(gps$time > t0 + 5 & gps$time < t0 + 11), ]
  m2 <- merge_streams(mag, gps2, max_gap = 1)
  expect_gt(attr(m2, "n_dropped"), 30)
  expect_error(merge_streams(mag, gps[1, ]), "Fewer than two")
  expect_error(merge_streams(mag, gps, lag = 1e6), "within `max_gap`")
})

test_that("assign_lines finds the lines of a racetrack pattern and names them by position", {
  sim <- sim_survey(seed = 8, defects = FALSE, n_lines = 6)
  files <- sim_raw_files(sim, withr::local_tempdir(), lag = 0)
  raw <- merge_streams(read_mag_log(files$mag, files$format, date = files$date), read_nmea(files$gps),
                       origin = files$origin)
  seg <- assign_lines(raw)
  s <- attr(seg, "segments")
  expect_equal(nrow(s), 6 + 3)
  expect_equal(sort(s$line[s$line_type == "traverse"]), sprintf("L%d", 1000 + (0:5) * 10))
  expect_equal(sort(s$line[s$line_type == "tie"]), sprintf("T%d", 9000 + (0:2) * 10))
  expect_gt(attr(seg, "n_turn_samples"), 0)
  expect_false(anyNA(seg$line))
  # names follow position: L1000 is the westernmost traverse, T9000 the southernmost tie
  expect_lt(mean(seg$x[seg$line == "L1000"]), mean(seg$x[seg$line == "L1010"]))
  expect_lt(mean(seg$y[seg$line == "T9000"]), mean(seg$y[seg$line == "T9010"]))
  # a plan renames the runs it matches
  plan <- data.frame(line = c("A1", "A2"), x0 = c(0, 200), y0 = 0, x1 = c(0, 200), y1 = 7000)
  sp <- attr(assign_lines(raw, plan = plan), "segments")
  expect_true(all(c("A1", "A2") %in% sp$line))
  expect_false("L1000" %in% sp$line)
  # keeping the turns
  kept <- assign_lines(raw, drop_turns = FALSE)
  expect_equal(sum(is.na(kept$line)), attr(seg, "n_turn_samples"))
  expect_error(assign_lines(raw, min_length = 1e6), "No straight run")
})

test_that("raw files round-trip to the simulated survey and its defects are still caught", {
  sim <- sim_survey(seed = 42)
  files <- sim_raw_files(sim, withr::local_tempdir(), lag = 0.3)
  mag <- read_mag_log(files$mag, files$format, date = files$date)
  gps <- read_nmea(files$gps)
  expect_equal(attr(gps, "n_bad_checksum"), 0)
  raw <- merge_streams(mag, gps, lag = files$lag, origin = files$origin)
  survey <- as_survey(assign_lines(raw))
  expect_equal(sort(unique(survey$line)), sort(unique(sim$data$line)))
  # positions agree with the simulated ones to well under a metre
  key <- function(line, time) paste(line, round(as.numeric(time), 1))
  j <- match(key(survey$line, survey$time), key(sim$data$line, sim$data$time))
  ok <- !is.na(j)
  expect_gt(mean(ok), 0.98)
  expect_lt(stats::quantile(sqrt((survey$x[ok] - sim$data$x[j[ok]])^2 + (survey$y[ok] - sim$data$y[j[ok]])^2), 0.99), 1)
  expect_equal(survey$mag_raw[ok], sim$data$mag_raw[j[ok]])
  # the QC sees the same defects
  res <- run_qc(survey, spec = sim$spec, base = sim$base)
  tr <- sim$truth[!is.na(sim$truth$line), ]
  for (chk in unique(tr$check)) {
    lines_expected <- unique(tr$line[tr$check == chk])
    lines_flagged <- unique(res$flags$line[res$flags$check == chk])
    expect_true(all(vapply(lines_expected, function(l) any(grepl(l, lines_flagged, fixed = TRUE)), logical(1))),
                info = chk)
  }
  expect_equal(nrow(res$flags[res$flags$check == "heading_error", ]), 1)
  expect_error(sim_raw_files(1, tempdir()), "sim_survey")
})
