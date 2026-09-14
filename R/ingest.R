# ---- NMEA --------------------------------------------------------------------

#' Read a raw GPS log (NMEA 0183)
#'
#' Parses `GGA` (position, altitude, fix quality), `RMC` (speed, course,
#' date) and `ZDA` (date) sentences from any talker (`GP`, `GN`, `GL`, `GA`,
#' `GB`). Sentences with a checksum are verified and dropped when it fails.
#' One row is returned per `GGA` fix, with the most recent `RMC` speed and
#' course attached when it is within a second.
#'
#' NMEA carries only time of day; the date comes from `RMC`/`ZDA` sentences
#' when the log has them, otherwise from `date`. A time of day that steps
#' back by more than twelve hours is taken as midnight passing.
#'
#' @param file Path to the log.
#' @param date A `Date` (or `"YYYY-MM-DD"`) for logs without date sentences.
#' @return A tibble: `time` (POSIXct, UTC), `lat`, `lon` (decimal degrees),
#'   `alt` (metres above the geoid, from `GGA`), `quality` (GGA fix
#'   quality), `n_sat`, `hdop`, `speed` (m/s, from RMC) and `course`
#'   (degrees, from RMC); attributes `n_sentences`, `n_bad_checksum`.
#' @examples
#' \dontrun{
#' gps <- read_nmea("flight.nmea")
#' }
#' @export
read_nmea <- function(file, date = NULL) {
  lines <- readLines(file, warn = FALSE)
  lines <- trimws(lines)
  n_sentences <- sum(startsWith(lines, "$"))
  keep <- grepl("^\\$G[PNLAB](GGA|RMC|ZDA),", lines)
  lines <- lines[keep]
  if (!length(lines)) stop("No GGA, RMC or ZDA sentences in ", file, ".", call. = FALSE)
  ok <- .nmea_checksum_ok(lines)
  n_bad <- sum(!ok)
  lines <- sub("\\*[0-9A-Fa-f]{2}$", "", lines[ok])
  fields <- strsplit(lines, ",", fixed = TRUE)
  type <- substr(vapply(fields, `[`, "", 1), 4, 6)

  num <- function(x) suppressWarnings(as.numeric(x))
  get <- function(f, k) if (length(f) >= k) f[k] else NA_character_
  tod <- vapply(fields, function(f) .nmea_tod(get(f, 2)), numeric(1))

  # a date per sentence: from ZDA / RMC where present, carried forward (and
  # back to the start), else the argument
  sdate <- rep(NA_character_, length(fields))
  is_rmc <- type == "RMC"; is_zda <- type == "ZDA"
  sdate[is_rmc] <- vapply(fields[is_rmc], function(f) {
    d <- get(f, 10); if (is.na(d) || nchar(d) != 6) NA_character_ else
      sprintf("20%s-%s-%s", substr(d, 5, 6), substr(d, 3, 4), substr(d, 1, 2))
  }, "")
  sdate[is_zda] <- vapply(fields[is_zda], function(f) {
    d <- get(f, 3); m <- get(f, 4); y <- get(f, 5)
    if (any(is.na(c(d, m, y))) || !nzchar(y)) NA_character_ else sprintf("%s-%s-%s", y, m, d)
  }, "")
  if (all(is.na(sdate))) {
    if (is.null(date)) stop("The log has no RMC/ZDA date sentences; pass `date`.", call. = FALSE)
    sdate[] <- as.character(as.Date(date))
  } else {
    sdate <- .fill_down(sdate)
  }
  base <- as.numeric(as.POSIXct(sdate, tz = "UTC"))
  time <- base + tod
  # midnight passing that the date sentences did not capture
  step <- diff(time)
  back <- which(step < -43200) + 1
  for (b in back) time[b:length(time)] <- time[b:length(time)] + 86400

  gga <- type == "GGA"
  fg <- fields[gga]
  lat <- vapply(fg, function(f) .nmea_deg(get(f, 3), get(f, 4)), numeric(1))
  lon <- vapply(fg, function(f) .nmea_deg(get(f, 5), get(f, 6)), numeric(1))
  out <- tibble::tibble(
    time = as.POSIXct(time[gga], origin = "1970-01-01", tz = "UTC"),
    lat = lat, lon = lon,
    alt = num(vapply(fg, function(f) get(f, 10), "")),
    quality = as.integer(num(vapply(fg, function(f) get(f, 7), ""))),
    n_sat = as.integer(num(vapply(fg, function(f) get(f, 8), ""))),
    hdop = num(vapply(fg, function(f) get(f, 9), "")),
    speed = NA_real_, course = NA_real_)

  if (any(is_rmc)) {
    fr <- fields[is_rmc]
    rmc <- data.frame(
      time = time[is_rmc],
      speed = num(vapply(fr, function(f) get(f, 8), "")) * 0.514444,   # knots to m/s
      course = num(vapply(fr, function(f) get(f, 9), "")))
    rmc <- rmc[order(rmc$time), ]
    k <- findInterval(as.numeric(out$time), rmc$time)
    k[k == 0] <- NA
    near <- !is.na(k) & abs(as.numeric(out$time) - rmc$time[pmax(k, 1)]) <= 1
    out$speed[near] <- rmc$speed[k[near]]
    out$course[near] <- rmc$course[k[near]]
  }
  out <- out[is.finite(out$lat) & is.finite(out$lon), ]
  out <- out[order(out$time), ]
  attr(out, "n_sentences") <- n_sentences
  attr(out, "n_bad_checksum") <- n_bad
  out
}

.nmea_checksum_ok <- function(lines) {
  has <- grepl("\\*[0-9A-Fa-f]{2}$", lines)
  ok <- rep(TRUE, length(lines))
  body <- sub("^\\$", "", sub("\\*[0-9A-Fa-f]{2}$", "", lines[has]))
  given <- toupper(sub(".*\\*", "", lines[has]))
  calc <- vapply(body, .nmea_xor, "")
  ok[has] <- calc == given
  ok
}

.nmea_xor <- function(s) {
  vapply(s, function(one) {
    x <- 0L
    for (b in utf8ToInt(one)) x <- bitwXor(x, b)
    sprintf("%02X", x)
  }, "", USE.NAMES = FALSE)
}

#' NMEA time of day (hhmmss.ss) in seconds
#' @noRd
.nmea_tod <- function(s) {
  v <- suppressWarnings(as.numeric(s))
  if (is.na(v)) return(NA_real_)
  h <- floor(v / 10000); m <- floor((v - h * 10000) / 100); sec <- v - h * 10000 - m * 100
  h * 3600 + m * 60 + sec
}

#' NMEA ddmm.mmmm plus hemisphere to decimal degrees
#' @noRd
.nmea_deg <- function(s, hemi) {
  v <- suppressWarnings(as.numeric(s))
  if (is.na(v)) return(NA_real_)
  d <- floor(v / 100); mn <- v - d * 100
  out <- d + mn / 60
  if (!is.na(hemi) && hemi %in% c("S", "W")) -out else out
}

.fill_down <- function(x) {
  if (all(is.na(x))) return(x)
  first <- which(!is.na(x))[1]
  x[seq_len(first - 1)] <- x[first]
  for (i in seq_along(x)[-1]) if (is.na(x[i])) x[i] <- x[i - 1]
  x
}

# ---- magnetometer logs --------------------------------------------------------

#' Describe a raw magnetometer log format
#'
#' Instrument logs are line-oriented text with the field reading, usually a
#' signal-quality value, and a time stamp in some vendor arrangement. A
#' `magqc_mag_format` says how to get from a line to those values:
#' which lines to keep (`prefix`), what to strip before splitting (`strip`),
#' the separator, which column holds what (`columns`), and how the time
#' column is written.
#'
#' Two presets are provided, written from the vendors' documentation rather
#' than against files in hand - verify them against your own logs:
#' `"geometrics"` for a G-88x serial stream captured with a time stamp
#' (`$ 52341.123,1234,08:00:00.100`), and `"gem"` for a GSM-19/GSMP text
#' dump (`52341.12 99 080000.1`, whitespace-separated).
#'
#' @param preset `"geometrics"`, `"gem"`, or `NULL` for a custom format.
#' @param sep Separator, as a regular expression (`","`, `"\\\\s+"`).
#' @param prefix Regular expression a line must match to be a data record
#'   (`NULL` keeps every non-blank, non-comment line).
#' @param strip Regular expression removed from each kept line before it is
#'   split (a leading `$`, say).
#' @param skip Header lines to skip.
#' @param comment Lines starting with this are ignored.
#' @param columns Named integer vector of column positions. `mag` (the field,
#'   nT) is required; `time` is required unless `time = "index"`. Any other
#'   name (`signal`, `radar_alt`, `temperature`, ...) is carried through as
#'   a column of that name.
#' @param time How the time column is written: `"hms"` (time of day,
#'   parsed with `time_format`), `"seconds"` (seconds since midnight),
#'   `"datetime"` (date and time in one field, parsed with `time_format`),
#'   `"elapsed"` (seconds since the start passed to [read_mag_log()]), or
#'   `"index"` (no time column; samples are `1 / rate` apart from the
#'   start).
#' @param time_format `strptime` format for `"hms"` and `"datetime"`.
#' @param rate Sample rate, Hz, for `time = "index"`.
#' @param na Strings that mean missing.
#' @return A `magqc_mag_format`.
#' @examples
#' mag_format("geometrics")
#' mag_format(sep = "\\s+", columns = c(time = 1, mag = 2, radar_alt = 5),
#'            time = "seconds")
#' @export
mag_format <- function(preset = NULL, sep = ",", prefix = NULL, strip = NULL, skip = 0L,
                       comment = "#", columns = c(time = 1, mag = 2),
                       time = c("hms", "seconds", "datetime", "elapsed", "index"),
                       time_format = "%H:%M:%OS", rate = NULL, na = c("", "NA", "*", "*****")) {
  time <- match.arg(time)
  if (!is.null(preset)) {
    preset <- match.arg(preset, c("geometrics", "gem"))
    return(switch(preset,
      geometrics = mag_format(sep = ",", prefix = "^\\$", strip = "^\\$\\s*",
                              columns = c(mag = 1, signal = 2, time = 3), time = "hms",
                              time_format = "%H:%M:%OS"),
      gem = mag_format(sep = "\\s+", prefix = "^\\s*[0-9]", columns = c(mag = 1, signal = 2, time = 3),
                       time = "hms", time_format = "%H%M%OS")))
  }
  if (!"mag" %in% names(columns)) stop("`columns` must name a `mag` column.", call. = FALSE)
  if (time != "index" && !"time" %in% names(columns)) stop("`columns` must name a `time` column unless `time = \"index\"`.", call. = FALSE)
  if (time == "index" && is.null(rate)) stop("`rate` is required for `time = \"index\"`.", call. = FALSE)
  structure(list(sep = sep, prefix = prefix, strip = strip, skip = as.integer(skip), comment = comment,
                 columns = columns, time = time, time_format = time_format, rate = rate, na = na),
            class = "magqc_mag_format")
}

#' @export
print.magqc_mag_format <- function(x, ...) {
  cat("<magqc raw magnetometer log format>\n")
  cat(sprintf("  columns: %s\n", paste(sprintf("%s = %d", names(x$columns), x$columns), collapse = ", ")))
  cat(sprintf("  sep '%s'%s%s; time: %s%s\n", x$sep,
              if (!is.null(x$prefix)) sprintf(", keep lines matching '%s'", x$prefix) else "",
              if (!is.null(x$strip)) sprintf(", strip '%s'", x$strip) else "",
              x$time, if (x$time %in% c("hms", "datetime")) sprintf(" (%s)", x$time_format) else
                if (x$time == "index") sprintf(" at %g Hz", x$rate) else ""))
  invisible(x)
}

#' Read a raw magnetometer log
#'
#' @param file Path to the log.
#' @param format A [mag_format()] (or a preset name).
#' @param date For `time = "hms"` and `"seconds"`, the date of the first
#'   record (`Date` or `"YYYY-MM-DD"`); midnight passing is detected. For
#'   `"elapsed"` and `"index"`, a POSIXct start time.
#' @return A tibble with `time` (POSIXct, UTC), `mag_raw` (nT) and any other
#'   columns the format names; attributes `n_lines`, `n_records`.
#' @examples
#' \dontrun{
#' mag <- read_mag_log("flight.mag", mag_format("geometrics"), date = "2024-06-14")
#' }
#' @export
read_mag_log <- function(file, format = mag_format("geometrics"), date = NULL) {
  if (is.character(format)) format <- mag_format(format)
  if (!inherits(format, "magqc_mag_format")) stop("`format` must come from mag_format().", call. = FALSE)
  f <- format
  lines <- readLines(file, warn = FALSE)
  n_lines <- length(lines)
  if (f$skip > 0) lines <- lines[-seq_len(min(f$skip, length(lines)))]
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  if (nzchar(f$comment)) lines <- lines[!startsWith(lines, f$comment)]
  if (!is.null(f$prefix)) lines <- lines[grepl(f$prefix, lines)]
  if (!is.null(f$strip)) lines <- sub(f$strip, "", lines)
  if (!length(lines)) stop("No data records in ", file, " with this format.", call. = FALSE)
  parts <- strsplit(lines, f$sep, perl = TRUE)
  need <- max(f$columns)
  parts <- parts[lengths(parts) >= need]
  if (!length(parts)) stop("No record in ", file, " has ", need, " columns.", call. = FALSE)
  col <- function(k) { v <- trimws(vapply(parts, `[`, "", k)); v[v %in% f$na] <- NA; v }
  out <- tibble::tibble(.rows = length(parts))
  for (nm in names(f$columns)) {
    v <- col(f$columns[[nm]])
    out[[if (nm == "mag") "mag_raw" else nm]] <- if (nm == "time" && f$time %in% c("hms", "datetime")) v else
      suppressWarnings(as.numeric(v))
  }
  out$time <- .mag_time(out$time, f, date, n = nrow(out))
  keep <- is.finite(out$mag_raw) & !is.na(out$time)
  out <- out[keep, ]
  out <- out[, c("time", "mag_raw", setdiff(names(out), c("time", "mag_raw")))]
  attr(out, "n_lines") <- n_lines
  attr(out, "n_records") <- nrow(out)
  out
}

.mag_time <- function(v, f, date, n) {
  need_date <- function() {
    if (is.null(date)) stop("This format's time is a time of day; pass `date`.", call. = FALSE)
    as.numeric(as.POSIXct(as.character(as.Date(date)), tz = "UTC"))
  }
  t <- switch(f$time,
    hms = {
      d0 <- need_date()
      tod <- as.numeric(as.POSIXct(paste("2000-01-01", v), format = paste("%Y-%m-%d", f$time_format), tz = "UTC")) -
        as.numeric(as.POSIXct("2000-01-01", tz = "UTC"))
      d0 + tod
    },
    seconds = need_date() + as.numeric(v),
    datetime = as.numeric(as.POSIXct(v, format = f$time_format, tz = "UTC")),
    elapsed = {
      if (is.null(date)) stop("`time = \"elapsed\"` needs the start time in `date`.", call. = FALSE)
      as.numeric(as.POSIXct(date, tz = "UTC")) + as.numeric(v)
    },
    index = {
      if (is.null(date)) stop("`time = \"index\"` needs the start time in `date`.", call. = FALSE)
      as.numeric(as.POSIXct(date, tz = "UTC")) + (seq_len(n) - 1) / f$rate
    })
  if (f$time %in% c("hms", "seconds")) {
    back <- which(diff(t) < -43200) + 1
    for (b in back) t[b:length(t)] <- t[b:length(t)] + 86400
  }
  as.POSIXct(t, origin = "1970-01-01", tz = "UTC")
}
