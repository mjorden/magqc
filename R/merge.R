#' Merge a magnetometer stream with a GPS stream
#'
#' Interpolates the GPS position linearly onto each magnetometer sample's
#' time, after shifting the GPS clock by `lag`. Magnetometer samples with no
#' GPS fix within `max_gap` seconds are dropped (an outage is not a place to
#' invent positions). Local east/north metres are added from `origin`, so
#' the result is ready for [assign_lines()] and then [as_survey()].
#'
#' @param mag Output of [read_mag_log()] (or any data frame with `time` and
#'   `mag_raw`).
#' @param gps Output of [read_nmea()] (or any data frame with `time`, `lat`,
#'   `lon` and optionally `alt`).
#' @param lag Magnetometer clock minus GPS clock, seconds: a positive lag
#'   means the magnetometer stamps are later than the GPS stamps for the
#'   same instant.
#' @param max_gap Longest distance in time, seconds, from a magnetometer
#'   sample to its nearest GPS fix for the sample to be kept.
#' @param origin `c(lon, lat)` of the local grid origin; the centroid of the
#'   fixes when `NULL`.
#' @return A tibble with `fid`, `time`, `lon`, `lat`, `gps_alt`, `x`, `y`,
#'   `mag_raw` and any other magnetometer columns; attributes `n_dropped`
#'   and `origin`.
#' @examples
#' \dontrun{
#' raw <- merge_streams(read_mag_log("flight.mag", date = "2024-06-14"),
#'                      read_nmea("flight.nmea"), lag = 0.3)
#' survey <- as_survey(assign_lines(raw))
#' }
#' @export
merge_streams <- function(mag, gps, lag = 0, max_gap = 2, origin = NULL) {
  for (nm in c("time", "mag_raw")) if (!nm %in% names(mag)) stop("`mag` needs a `", nm, "` column.", call. = FALSE)
  for (nm in c("time", "lat", "lon")) if (!nm %in% names(gps)) stop("`gps` needs a `", nm, "` column.", call. = FALSE)
  if (!inherits(mag$time, "POSIXct") || !inherits(gps$time, "POSIXct")) stop("`time` columns must be POSIXct.", call. = FALSE)
  gps <- gps[is.finite(gps$lat) & is.finite(gps$lon), ]
  gps <- gps[order(gps$time), ]
  gps <- gps[!duplicated(gps$time), ]
  if (nrow(gps) < 2) stop("Fewer than two usable GPS fixes.", call. = FALSE)
  mag <- mag[order(mag$time), ]

  tg <- as.numeric(gps$time) + lag
  tm <- as.numeric(mag$time)
  k <- findInterval(tm, tg)
  gap <- pmin(abs(tm - tg[pmax(k, 1)]), abs(tm - tg[pmin(k + 1, length(tg))]))
  keep <- tm >= tg[1] & tm <= tg[length(tg)] & gap <= max_gap
  n_dropped <- sum(!keep)
  if (!any(keep)) stop("No magnetometer sample has a GPS fix within `max_gap`; check `lag` and the dates.", call. = FALSE)
  out <- mag[keep, ]
  tmk <- tm[keep]
  out$lon <- stats::approx(tg, gps$lon, xout = tmk, ties = mean)$y
  out$lat <- stats::approx(tg, gps$lat, xout = tmk, ties = mean)$y
  out$gps_alt <- if ("alt" %in% names(gps) && any(is.finite(gps$alt))) stats::approx(tg, gps$alt, xout = tmk, ties = mean)$y else NA_real_
  origin <- origin %||% c(lon = mean(out$lon), lat = mean(out$lat))
  xy <- lonlat_to_enu(out$lon, out$lat, origin[["lon"]], origin[["lat"]])
  out$x <- xy$x; out$y <- xy$y
  out$fid <- seq_len(nrow(out))
  out <- out[, c("fid", "time", "lon", "lat", "gps_alt", "x", "y", "mag_raw",
                 setdiff(names(out), c("fid", "time", "lon", "lat", "gps_alt", "x", "y", "mag_raw")))]
  attr(out, "n_dropped") <- n_dropped
  attr(out, "origin") <- origin
  tibble::as_tibble(out)
}

#' Segment a continuous stream into flight lines
#'
#' A raw stream has no line numbers: the aircraft flies a line, turns, flies
#' the next. Straight flight is found from the heading: positions are
#' smoothed over `window` seconds, the turn rate is the smoothed heading's
#' rate of change, and runs where it stays below `max_turn_rate` (and no
#' time gap longer than `max_gap` intervenes) are candidate lines; runs
#' shorter than `min_length` metres are turns, and consecutive collinear
#' runs are rejoined so a data gap stays inside its line. The traverse
#' axis is the heading of the longest run (or `azimuth`); runs within
#' `max_off_axis` degrees of it (either way) are traverses, within it of the
#' perpendicular ties, and anything else - a transit between blocks - is
#' dropped. Traverses are numbered west to east or south to north across
#' the block (`L1000`, `L1010`, ...), ties along it (`T9000`, ...). A
#' flight-plan table renames runs to the plan line whose axis passes within
#' `tolerance` metres of the run's centroid.
#'
#' @param df Output of [merge_streams()] (needs `time`, `x`, `y`).
#' @param azimuth Traverse azimuth, degrees; estimated when `NULL`.
#' @param window Smoothing window for heading, seconds.
#' @param max_turn_rate Turn rate above which flight is a turn, degrees per
#'   second. A racetrack turn (180 degrees in about 90 s) is 2 deg/s;
#'   straight flight with GPS wander stays under 0.3.
#' @param min_length Shortest run kept as a line, metres.
#' @param max_gap A time gap longer than this, seconds, ends a run.
#' @param max_spread A run whose heading wanders more than this, degrees,
#'   from its mean is a turn however slow it was flown.
#' @param max_off_axis A run is a traverse when its heading is within this
#'   many degrees of the traverse axis (either way), a tie when within it of
#'   the perpendicular, and otherwise a transit between blocks, which is
#'   dropped.
#' @param plan Optional flight plan: `line`, `x0`, `y0`, `x1`, `y1`
#'   (local metres) and optionally `line_type`.
#' @param tolerance Centroid-to-plan-line distance for a plan match, metres.
#' @param drop_turns Drop the turn samples (`TRUE`) or keep them with
#'   `line = NA`.
#' @return `df` with `line` and `line_type` columns (and `dir`, the run's
#'   mean heading); attribute `segments` summarises each run (`line`,
#'   `line_type`, `n`, `length_m`, `heading`).
#' @export
assign_lines <- function(df, azimuth = NULL, window = 3, max_turn_rate = 1, min_length = 500,
                         max_gap = 10, max_spread = 30, max_off_axis = 15, plan = NULL, tolerance = NULL,
                         drop_turns = TRUE) {
  for (nm in c("time", "x", "y")) if (!nm %in% names(df)) stop("`df` needs a `", nm, "` column.", call. = FALSE)
  df <- df[order(df$time), ]
  n <- nrow(df)
  if (n < 10) stop("Too few samples to segment.", call. = FALSE)
  t <- as.numeric(df$time)
  dt <- stats::median(diff(t))
  if (!is.finite(dt) || dt <= 0) stop("Timestamps do not advance.", call. = FALSE)
  k <- max(3L, as.integer(round(window / dt)))
  if (k %% 2 == 0) k <- k + 1L

  sx <- .runmean(df$x, k); sy <- .runmean(df$y, k)
  hdg <- .bearing(c(diff(sx), NA), c(diff(sy), NA))
  hdg[n] <- hdg[n - 1]
  # rate of change of an unwrapped heading, smoothed again
  dh <- diff(hdg); dh <- ((dh + 180) %% 360) - 180
  rate <- c(abs(dh) / diff(t), 0)
  rate <- .runmean(rate, k)
  gap_break <- c(FALSE, diff(t) > max_gap)
  straight <- is.finite(rate) & rate < max_turn_rate

  # runs of straight flight, broken by time gaps
  run <- cumsum(!straight | gap_break | c(TRUE, !straight[-n]))
  run[!straight] <- NA
  ids <- unique(stats::na.omit(run))
  segs <- lapply(ids, function(id) {
    i <- which(run == id)
    len <- sum(sqrt(diff(df$x[i])^2 + diff(df$y[i])^2))
    h <- .mean_bearing(hdg[i])
    spread <- max(abs(((hdg[i] - h + 180) %% 360) - 180), na.rm = TRUE)
    list(idx = i, length = len, heading = h, spread = spread, cx = mean(df$x[i]), cy = mean(df$y[i]))
  })
  segs <- segs[vapply(segs, function(s) s$spread <= max_spread, logical(1))]

  # A data gap or a brief wobble inside a line makes two runs of one line:
  # rejoin consecutive runs that are collinear (same heading, the later
  # one's centroid on the earlier one's axis) so the gap stays inside the
  # line for the data-gap check to find.
  if (length(segs) > 1) {
    merged <- list(segs[[1]])
    for (s in segs[-1]) {
      cur <- merged[[length(merged)]]
      dh <- abs(((s$heading - cur$heading + 180) %% 360) - 180)
      ar <- .deg2rad(cur$heading)
      off <- abs((s$cx - cur$cx) * cos(ar) - (s$cy - cur$cy) * sin(ar))   # across-axis distance
      if (dh <= 10 && off <= 50) {
        i <- c(cur$idx, s$idx)
        merged[[length(merged)]] <- list(
          idx = i, length = cur$length + s$length, heading = .mean_bearing(hdg[i]),
          spread = max(cur$spread, s$spread), cx = mean(df$x[i]), cy = mean(df$y[i]))
      } else merged[[length(merged) + 1]] <- s
    }
    segs <- merged
  }
  segs <- segs[vapply(segs, function(s) s$length >= min_length, logical(1))]
  if (!length(segs)) stop("No straight run of at least ", min_length, " m found; check `max_turn_rate` and `min_length`.", call. = FALSE)

  # the traverse axis as a direction in [0, 180): which way the longest run
  # happened to be flown must not decide the numbering direction
  axis <- (azimuth %||% segs[[which.max(vapply(segs, function(s) s$length, numeric(1)))]]$heading) %% 180
  off_axis <- function(h) abs(((h - axis + 90) %% 180) - 90)   # 0 along the axis, 90 across
  # straight flight that is neither along the traverse axis nor across it is
  # a transit between blocks, not survey data
  for (s in seq_along(segs)) {
    d <- off_axis(segs[[s]]$heading)
    segs[[s]]$line_type <- if (d <= max_off_axis) "traverse" else if (d >= 90 - max_off_axis) "tie" else "transit"
  }
  segs <- segs[vapply(segs, function(s) s$line_type != "transit", logical(1))]
  if (!length(segs)) stop("No run lies along or across the traverse axis; pass `azimuth`.", call. = FALSE)

  # naming: traverses ordered across the block, ties along it; both
  # orderings run west-to-east or south-to-north
  ar <- .deg2rad(axis); along <- c(sin(ar), cos(ar)); across <- c(cos(ar), -sin(ar))
  # the dominant component decides: a near-north axis orders south to north,
  # a near-east axis west to east (a tiny negative minor component must not)
  eastnorth <- function(v) if (v[which.max(abs(v))] < 0) -v else v
  along <- eastnorth(along); across <- eastnorth(across)
  pos <- vapply(segs, function(s) if (s$line_type == "traverse") s$cx * across[1] + s$cy * across[2] else
    s$cx * along[1] + s$cy * along[2], numeric(1))
  for (typ in c("traverse", "tie")) {
    j <- which(vapply(segs, function(s) s$line_type == typ, logical(1)))
    if (!length(j)) next
    o <- j[order(pos[j])]
    for (r in seq_along(o)) segs[[o[r]]]$line <- sprintf("%s%d", if (typ == "traverse") "L" else "T",
                                                        (if (typ == "traverse") 1000 else 9000) + (r - 1) * 10)
  }

  if (!is.null(plan)) {
    tolerance <- tolerance %||% 50
    for (s in seq_along(segs)) {
      d <- .point_line_distance(segs[[s]]$cx, segs[[s]]$cy, plan$x0, plan$y0, plan$x1, plan$y1)
      m <- which.min(d)
      if (length(m) && d[m] <= tolerance) {
        segs[[s]]$line <- as.character(plan$line[m])
        if ("line_type" %in% names(plan)) segs[[s]]$line_type <- as.character(plan$line_type[m])
      }
    }
  }

  df$line <- NA_character_; df$line_type <- NA_character_; df$dir <- NA_character_
  for (s in segs) {
    df$line[s$idx] <- s$line; df$line_type[s$idx] <- s$line_type
    df$dir[s$idx] <- .compass(s$heading)
  }
  summary <- tibble::tibble(
    line = vapply(segs, `[[`, "", "line"), line_type = vapply(segs, `[[`, "", "line_type"),
    n = vapply(segs, function(s) length(s$idx), integer(1)),
    length_m = vapply(segs, `[[`, numeric(1), "length"),
    heading = vapply(segs, `[[`, numeric(1), "heading"))
  n_turn <- sum(is.na(df$line))
  if (drop_turns) df <- df[!is.na(df$line), ]
  attr(df, "segments") <- summary
  attr(df, "n_turn_samples") <- n_turn
  df
}

#' Running mean with edge shrinkage
#' @noRd
.runmean <- function(x, k) {
  n <- length(x); if (n < k) return(x)
  cs <- cumsum(c(0, x)); h <- (k - 1) %/% 2
  lo <- pmax(seq_len(n) - h, 1); hi <- pmin(seq_len(n) + h, n)
  (cs[hi + 1] - cs[lo]) / (hi - lo + 1)
}

#' Perpendicular distance from a point to the infinite line through two points
#' @noRd
.point_line_distance <- function(px, py, x0, y0, x1, y1) {
  dx <- x1 - x0; dy <- y1 - y0
  len <- sqrt(dx^2 + dy^2)
  d <- abs(dy * (px - x0) - dx * (py - y0)) / len
  d[!is.finite(d)] <- Inf
  d
}
