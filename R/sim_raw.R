#' Write a simulated survey as raw instrument and GPS files
#'
#' Turns a [sim_survey()] into the two files an acquisition system would
#' leave behind - a magnetometer log and an NMEA GPS log - including the
#' turns between lines, which the simulator itself does not fly. The
#' magnetometer log is written in the `"geometrics"` [mag_format()] layout
#' with one extra column (`$ field,signal,hh:mm:ss.sss,radar_alt`); the GPS
#' log has a `GGA` fix at `gps_rate` Hz and an `RMC` (date, speed, course)
#' once a second, all with checksums, its clock `lag` seconds behind the
#' magnetometer's.
#'
#' Reading the pair back with [read_mag_log()], [read_nmea()],
#' [merge_streams()] and [assign_lines()] should reproduce the simulated
#' survey - that round trip is the ingestion test.
#'
#' @param sim A `magqc_sim`.
#' @param dir Directory to write into.
#' @param gps_rate GPS fix rate, Hz.
#' @param lag Magnetometer clock minus GPS clock, seconds.
#' @param turn_time Seconds spent turning between lines.
#' @return A list with `mag` and `gps` (the file paths), `format` (the
#'   [mag_format()] that reads the magnetometer log), `date`, `lag` and
#'   `origin`.
#' @examples
#' \dontrun{
#' files <- sim_raw_files(sim_survey(), tempdir(), lag = 0.3)
#' mag <- read_mag_log(files$mag, files$format, date = files$date)
#' gps <- read_nmea(files$gps)
#' survey <- as_survey(assign_lines(merge_streams(mag, gps, lag = files$lag, origin = files$origin)))
#' }
#' @export
sim_raw_files <- function(sim, dir, gps_rate = 5, lag = 0, turn_time = 90) {
  if (!inherits(sim, "magqc_sim")) stop("`sim` must come from sim_survey().", call. = FALSE)
  d <- as.data.frame(sim$data)
  d <- d[order(d$time), ]
  origin <- sim$origin
  dt <- stats::median(diff(as.numeric(d$time)))

  # ---- the flown path: lines as simulated, turns synthesised between them
  pieces <- list(); prev_end <- NULL
  for (idx in .by_line(d)) {
    ln <- d[idx, ]
    if (!is.null(prev_end)) {
      p0 <- prev_end; p1 <- ln[1, ]
      # cubic Bezier from the end of the last line to the start of this one,
      # leaving and arriving along the line directions
      d0 <- c(p0$x - p0$x_prev, p0$y - p0$y_prev); d0 <- d0 / sqrt(sum(d0^2))
      d1 <- c(ln$x[2] - ln$x[1], ln$y[2] - ln$y[1]); d1 <- d1 / sqrt(sum(d1^2))
      chord <- sqrt((p1$x - p0$x)^2 + (p1$y - p0$y)^2)
      L <- max(chord * 0.6, 150)
      c0 <- c(p0$x, p0$y) + d0 * L; c1 <- c(p1$x, p1$y) - d1 * L
      nt <- max(2L, as.integer(round(turn_time / dt)) - 1L)
      s <- seq(0, 1, length.out = nt + 2)[-c(1, nt + 2)]
      bx <- (1 - s)^3 * p0$x + 3 * (1 - s)^2 * s * c0[1] + 3 * (1 - s) * s^2 * c1[1] + s^3 * p1$x
      by <- (1 - s)^3 * p0$y + 3 * (1 - s)^2 * s * c0[2] + 3 * (1 - s) * s^2 * c1[2] + s^3 * p1$y
      tt <- as.numeric(p0$time) + s * (as.numeric(p1$time) - as.numeric(p0$time))
      pieces[[length(pieces) + 1]] <- data.frame(
        time = tt, x = bx, y = by,
        mag_raw = p0$mag_raw + s * (p1$mag_raw - p0$mag_raw) + stats::rnorm(nt, 0, 0.02),
        radar_alt = p0$radar_alt + s * (p1$radar_alt - p0$radar_alt),
        gps_alt = p0$gps_alt + s * (p1$gps_alt - p0$gps_alt))
    }
    pieces[[length(pieces) + 1]] <- data.frame(
      time = as.numeric(ln$time), x = ln$x, y = ln$y, mag_raw = ln$mag_raw,
      radar_alt = ln$radar_alt, gps_alt = if ("gps_alt" %in% names(ln)) ln$gps_alt else ln$radar_alt + 1500)
    prev_end <- ln[nrow(ln), ]; prev_end$x_prev <- ln$x[nrow(ln) - 1]; prev_end$y_prev <- ln$y[nrow(ln) - 1]
  }
  path <- do.call(rbind, pieces)
  path <- path[order(path$time), ]
  ll <- enu_to_lonlat(path$x, path$y, origin[["lon"]], origin[["lat"]])
  path$lon <- ll$lon; path$lat <- ll$lat
  date <- as.Date(as.POSIXct(path$time[1], origin = "1970-01-01", tz = "UTC"))

  # ---- magnetometer log
  tod <- function(t) .hms_str(t, 3, sep = ":")
  mag_file <- file.path(dir, "flight.mag")
  writeLines(c("# simulated magnetometer log: $ field_nT,signal,time_utc,radar_alt_m",
               sprintf("$ %.3f,%d,%s,%.1f", path$mag_raw, 1200L + sample(0:99, nrow(path), TRUE),
                       tod(path$time), path$radar_alt)), mag_file)

  # ---- NMEA log, its clock `lag` seconds behind
  tg <- seq(min(path$time), max(path$time), by = 1 / gps_rate)
  glon <- stats::approx(path$time, path$lon, xout = tg)$y
  glat <- stats::approx(path$time, path$lat, xout = tg)$y
  galt <- stats::approx(path$time, path$gps_alt, xout = tg)$y
  gx <- stats::approx(path$time, path$x, xout = tg)$y; gy <- stats::approx(path$time, path$y, xout = tg)$y
  spd <- c(sqrt(diff(gx)^2 + diff(gy)^2) * gps_rate, NA); spd[length(spd)] <- spd[length(spd) - 1]
  crs <- .bearing(c(diff(gx), NA), c(diff(gy), NA)); crs[length(crs)] <- crs[length(crs) - 1]
  tgs <- tg - lag
  hms <- function(t) .hms_str(t, 2, sep = "")
  ddmmyy <- format(as.POSIXct(tgs, origin = "1970-01-01", tz = "UTC"), "%d%m%y")
  nmea <- function(body) sprintf("$%s*%s", body, .nmea_xor(body))
  gga <- nmea(sprintf("GPGGA,%s,%s,%s,%s,%s,1,%02d,%.1f,%.1f,M,-17.0,M,,", hms(tgs),
                      .nmea_fmt(abs(glat), 2), ifelse(glat >= 0, "N", "S"),
                      .nmea_fmt(abs(glon), 3), ifelse(glon >= 0, "E", "W"),
                      9L + sample(0:3, length(tg), TRUE), 0.8 + stats::runif(length(tg), 0, 0.4), galt))
  whole <- which(abs(tgs - round(tgs)) < 1e-6 | seq_along(tgs) %% gps_rate == 1)
  rmc <- nmea(sprintf("GPRMC,%s,A,%s,%s,%s,%s,%.1f,%.1f,%s,,,A", hms(tgs[whole]),
                      .nmea_fmt(abs(glat[whole]), 2), ifelse(glat[whole] >= 0, "N", "S"),
                      .nmea_fmt(abs(glon[whole]), 3), ifelse(glon[whole] >= 0, "E", "W"),
                      spd[whole] / 0.514444, crs[whole], ddmmyy[whole]))
  gps_file <- file.path(dir, "flight.nmea")
  order_idx <- order(c(tgs, tgs[whole]), c(rep(1, length(tgs)), rep(0, length(whole))))
  writeLines(c(gga, rmc)[order_idx], gps_file)

  list(mag = mag_file, gps = gps_file,
       format = mag_format(sep = ",", prefix = "^\\$", strip = "^\\$\\s*",
                           columns = c(mag = 1, signal = 2, time = 3, radar_alt = 4),
                           time = "hms", time_format = "%H:%M:%OS"),
       date = date, lag = lag, origin = origin)
}

#' Time of day as hh:mm:ss.sss, rounded
#'
#' `format(x, "%OSn")` truncates fractional seconds, so 08:00:00.1 stored as
#' 0.0999... prints as .099; round the numeric seconds first.
#' @noRd
.hms_str <- function(t, digits, sep = ":") {
  tod <- round(t %% 86400, digits)
  h <- floor(tod / 3600); m <- floor((tod - h * 3600) / 60); s <- tod - h * 3600 - m * 60
  sprintf(paste0("%02d", sep, "%02d", sep, "%0", 3 + digits, ".", digits, "f"), as.integer(h), as.integer(m), s)
}

#' Decimal degrees to NMEA ddmm.mmmm (or dddmm.mmmm)
#' @noRd
.nmea_fmt <- function(deg, dwidth) {
  d <- floor(deg); m <- (deg - d) * 60
  sprintf(paste0("%0", dwidth, "d%07.4f"), as.integer(d), m)
}
