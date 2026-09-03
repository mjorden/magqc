# Every check returns the same flag table so the report, the map and the
# scorecard can treat them interchangeably. One row is one flagged interval:
# `i_start`/`i_end` index into the survey (`.i`), `value` is the statistic that
# tripped the check and `threshold` is the limit it was compared against.
# Survey-level findings (line separation, heading error) have NA `i_start`.
#
# A check also attaches a `metric` attribute - the single summary statistic
# the scorecard reports, with its label and limit.

.flag_cols <- function() {
  tibble::tibble(
    check = character(), line = character(),
    i_start = integer(), i_end = integer(),
    fid_start = numeric(), fid_end = numeric(),
    time_start = as.POSIXct(character(), tz = "UTC"),
    time_end = as.POSIXct(character(), tz = "UTC"),
    x = numeric(), y = numeric(), lon = numeric(), lat = numeric(),
    value = numeric(), threshold = numeric(), units = character(),
    description = character())
}

.flag_runs <- function(survey, idx, runs, check, value, threshold, units, description) {
  if (!nrow(runs)) return(.flag_cols())
  s <- idx[runs$start]; e <- idx[runs$end]
  tibble::tibble(
    check = check, line = survey$line[s],
    i_start = as.integer(s), i_end = as.integer(e),
    fid_start = survey$fid[s], fid_end = survey$fid[e],
    time_start = survey$time[s], time_end = survey$time[e],
    x = (survey$x[s] + survey$x[e]) / 2, y = (survey$y[s] + survey$y[e]) / 2,
    lon = (survey$lon[s] + survey$lon[e]) / 2, lat = (survey$lat[s] + survey$lat[e]) / 2,
    value = value, threshold = threshold, units = units,
    description = description)
}

.bind_flags <- function(x) {
  x <- Filter(Negate(is.null), x)
  if (!length(x)) return(.flag_cols())
  dplyr::bind_rows(x)
}

.with_metric <- function(flags, label, value, threshold, note = NULL) {
  attr(flags, "metric") <- list(label = label, value = value,
                                threshold = threshold, note = note)
  flags
}

.run_extreme <- function(v, runs) {
  vapply(seq_len(nrow(runs)), function(r) {
    s <- v[runs$start[r]:runs$end[r]]
    s <- s[is.finite(s)]
    if (!length(s)) return(NA_real_)
    s[which.max(abs(s))]
  }, numeric(1))
}

.by_line <- function(survey) split(survey$.i, survey$line)

#' Continuous segments: a line, broken wherever the record gaps
#'
#' The noise tests difference consecutive samples, so a data gap - across
#' which the field has moved by hundreds of metres of gradient - would read as
#' a giant spike. Gaps are the gap check's business; the noise checks see
#' each continuous stretch on its own.
#' @noRd
.by_segment <- function(survey, max_gap) {
  out <- list()
  for (idx in .by_line(survey)) {
    if (length(idx) < 2) { out[[length(out) + 1]] <- idx; next }
    brk <- c(0, cumsum(diff(as.numeric(survey$time[idx])) > max_gap))
    for (seg in split(idx, brk)) out[[length(out) + 1]] <- seg
  }
  out
}

# ---- data gaps ---------------------------------------------------------------

#' QC checks
#'
#' Each `check_*()` evaluates one acceptance criterion from a
#' [survey_spec()] against a [as_survey()] table and returns a tibble of
#' flagged intervals (see Details). [run_qc()] runs them all; call them
#' individually to re-test one criterion with a different threshold.
#'
#' @details
#' The returned tibble has one row per flagged interval with columns `check`,
#' `line`, `i_start`, `i_end` (row indices into the survey), `fid_start`,
#' `fid_end`, `time_start`, `time_end`, `x`, `y`, `lon`, `lat` (interval
#' midpoint), `value`, `threshold`, `units` and `description`. It carries a
#' `metric` attribute holding the summary statistic reported on the scorecard.
#'
#' * `check_sample_interval()` - breaks in the data stream longer than
#'   `max_gap` seconds.
#' * `check_along_line_spacing()` - ground distance between consecutive
#'   samples exceeding `max_along_line_sep`, excluding samples on either side
#'   of a gap (those are reported by the gap check).
#' * `check_line_deviation()` - lateral departure from the planned line, or
#'   from a total-least-squares fit through the flown line when no `plan` is
#'   given. Without a plan a whole line flown parallel to but off its planned
#'   position is invisible to this check and is caught by
#'   `check_line_separation()` instead.
#' * `check_line_separation()` - separation between adjacent traverse lines,
#'   measured perpendicular to `line_azimuth`, against `line_spacing` +/-
#'   `line_spacing_tol`.
#' * `check_clearance()` - radar altimeter against `nominal_clearance` +/-
#'   `clearance_tol`.
#' * `check_spikes()` - samples whose fourth difference exceeds `spike_nsigma`
#'   running MADs of the fourth difference. The test is made on the fourth
#'   difference rather than the field itself because a running-median filter
#'   on a smooth, steep anomaly reads ordinary noise as outliers; the fourth
#'   difference removes the geology first.
#' * `check_fourth_difference()` - the fourth-difference noise envelope on
#'   the despiked field (see [fourth_difference()]). Both noise tests are
#'   evaluated on each continuous stretch of a line separately: a data gap
#'   breaks the series, so the field jump across it is not mistaken for a
#'   spike.
#' * `check_diurnal()` - base-station departure from a linear chord over
#'   `diurnal_window` seconds; survey samples acquired inside an exceedance
#'   window are flagged.
#' * `check_crossovers()` - traverse/tie misfit at each intersection (see
#'   [crossovers()]); RMS misfit is the metric, individual intersections above
#'   `max_crossover_abs` are flagged.
#' * `check_heading_error()` - difference in mean crossover misfit between
#'   traverses flown in opposite directions.
#'
#' @param survey A `magqc_survey` from [as_survey()].
#' @param spec A [survey_spec()].
#' @param plan Optional flight plan: a data frame with columns `line`, `x0`,
#'   `y0`, `x1`, `y1` giving the planned start and end of each line.
#' @param base Base-station record: a data frame with `time` (POSIXct) and
#'   `mag_base` (nT).
#' @param xo Crossover table from [crossovers()].
#' @name checks
NULL

#' @rdname checks
#' @export
check_sample_interval <- function(survey, spec) {
  out <- list(); longest <- 0
  for (idx in .by_line(survey)) {
    if (length(idx) < 2) next
    gap <- c(NA_real_, diff(as.numeric(survey$time[idx])))
    longest <- max(longest, max(gap, na.rm = TRUE))
    bad <- which(gap > spec$max_gap)
    if (!length(bad)) next
    out[[length(out) + 1]] <- .flag_runs(
      survey, idx, data.frame(start = bad - 1L, end = bad), "sample_interval",
      value = gap[bad], threshold = spec$max_gap, units = "s",
      description = sprintf("%.1f s break in data stream", gap[bad]))
  }
  .with_metric(.bind_flags(out), "longest gap (s)", longest, spec$max_gap)
}

# ---- along-line spacing -------------------------------------------------------

#' @rdname checks
#' @export
check_along_line_spacing <- function(survey, spec) {
  out <- list(); worst <- 0
  for (idx in .by_line(survey)) {
    if (length(idx) < 2) next
    step <- .step_distance(survey$x[idx], survey$y[idx])
    dt <- c(NA_real_, diff(as.numeric(survey$time[idx])))
    ok <- is.finite(step) & is.finite(dt) & dt <= spec$max_gap
    worst <- max(worst, step[ok], na.rm = TRUE)
    bad <- ok & step > spec$max_along_line_sep
    runs <- .runs(bad)
    if (!nrow(runs)) next
    val <- .run_extreme(step, runs)
    runs$start <- pmax(1L, runs$start - 1L)
    out[[length(out) + 1]] <- .flag_runs(
      survey, idx, runs, "along_line_spacing",
      value = val, threshold = spec$max_along_line_sep, units = "m",
      description = sprintf("samples up to %.1f m apart (limit %.1f m)",
                            val, spec$max_along_line_sep))
  }
  .with_metric(.bind_flags(out), "max sample spacing (m)", worst, spec$max_along_line_sep)
}

# ---- departure from line ------------------------------------------------------

#' @rdname checks
#' @export
check_line_deviation <- function(survey, spec, plan = NULL) {
  out <- list(); all_off <- numeric(0)
  for (idx in .by_line(survey)) {
    if (length(idx) < 3) next
    ln <- survey$line[idx[1]]
    x <- survey$x[idx]; y <- survey$y[idx]
    planned <- !is.null(plan) && ln %in% plan$line
    if (planned) {
      p <- plan[plan$line == ln, ][1, ]
      dvec <- c(p$x1 - p$x0, p$y1 - p$y0); dvec <- dvec / sqrt(sum(dvec^2))
      nrm <- c(-dvec[2], dvec[1])
      off <- (x - p$x0) * nrm[1] + (y - p$y0) * nrm[2]
    } else {
      off <- .tls_line(x, y)$offset
    }
    all_off <- c(all_off, abs(off))
    runs <- .runs(abs(off) > spec$max_deviation)
    if (!nrow(runs)) next
    val <- .run_extreme(off, runs)
    out[[length(out) + 1]] <- .flag_runs(
      survey, idx, runs, "line_deviation",
      value = abs(val), threshold = spec$max_deviation, units = "m",
      description = sprintf("%.0f m off %s line", abs(val),
                            if (planned) "planned" else "best-fit"))
  }
  p95 <- if (length(all_off)) unname(stats::quantile(all_off, 0.95)) else NA_real_
  .with_metric(.bind_flags(out), "95th pct departure (m)", p95, spec$max_deviation,
               note = if (is.null(plan)) "measured against a best-fit line; no flight plan supplied")
}

# ---- line separation ----------------------------------------------------------

#' @rdname checks
#' @export
check_line_separation <- function(survey, spec) {
  tr <- survey[survey$line_type == "traverse", ]
  empty <- .with_metric(.flag_cols(), "max separation error (%)", NA_real_,
                        spec$line_spacing_tol * 100)
  if (!nrow(tr)) return(empty)
  cent <- dplyr::summarise(dplyr::group_by(tr, .data$line),
                           cx = mean(.data$x), cy = mean(.data$y),
                           lon = mean(.data$lon), lat = mean(.data$lat),
                           .groups = "drop")
  if (nrow(cent) < 2) return(empty)
  az <- .deg2rad(spec$line_azimuth)
  nrm <- c(cos(az), -sin(az))
  cent$pos <- cent$cx * nrm[1] + cent$cy * nrm[2]
  cent <- cent[order(cent$pos), ]
  sep <- diff(cent$pos)
  err_pct <- (sep - spec$line_spacing) / spec$line_spacing * 100
  bad <- which(abs(err_pct) > spec$line_spacing_tol * 100)
  flags <- if (length(bad)) {
    a <- bad; b <- bad + 1
    tibble::tibble(
      check = "line_separation",
      line = paste(cent$line[a], cent$line[b], sep = " / "),
      i_start = NA_integer_, i_end = NA_integer_,
      fid_start = NA_real_, fid_end = NA_real_,
      time_start = as.POSIXct(NA, tz = "UTC"), time_end = as.POSIXct(NA, tz = "UTC"),
      x = (cent$cx[a] + cent$cx[b]) / 2, y = (cent$cy[a] + cent$cy[b]) / 2,
      lon = (cent$lon[a] + cent$lon[b]) / 2, lat = (cent$lat[a] + cent$lat[b]) / 2,
      value = sep[bad], threshold = spec$line_spacing, units = "m",
      description = sprintf("%s to %s: %.0f m apart (nominal %g m, %+.0f%%)",
                            cent$line[a], cent$line[b], sep[bad],
                            spec$line_spacing, err_pct[bad]))
  } else .flag_cols()
  .with_metric(flags, "max separation error (%)", max(abs(err_pct)),
               spec$line_spacing_tol * 100)
}

# ---- terrain clearance --------------------------------------------------------

#' @rdname checks
#' @export
check_clearance <- function(survey, spec) {
  if (all(is.na(survey$radar_alt))) {
    return(.with_metric(.flag_cols(), "95th pct |clearance error| (m)", NA_real_,
                        spec$clearance_tol, note = "no radar altimeter column"))
  }
  out <- list(); all_dev <- numeric(0)
  for (idx in .by_line(survey)) {
    dev <- survey$radar_alt[idx] - spec$nominal_clearance
    all_dev <- c(all_dev, abs(dev[is.finite(dev)]))
    runs <- .runs(is.finite(dev) & abs(dev) > spec$clearance_tol)
    if (!nrow(runs)) next
    val <- .run_extreme(dev, runs)
    out[[length(out) + 1]] <- .flag_runs(
      survey, idx, runs, "clearance",
      value = val, threshold = spec$clearance_tol, units = "m",
      description = sprintf("%.0f m %s nominal clearance", abs(val),
                            ifelse(val > 0, "above", "below")))
  }
  p95 <- if (length(all_dev)) unname(stats::quantile(all_dev, 0.95)) else NA_real_
  .with_metric(.bind_flags(out), "95th pct |clearance error| (m)", p95, spec$clearance_tol)
}

# ---- spikes ----------------------------------------------------------------------

#' @rdname checks
#' @export
check_spikes <- function(survey, spec) {
  out <- list(); n_spikes <- 0L
  for (idx in .by_segment(survey, spec$max_gap)) {
    h <- .spike_detect(survey$mag[idx], k = spec$spike_window, nsigma = spec$spike_nsigma)
    runs <- .runs(h$outlier)
    if (!nrow(runs)) next
    n_spikes <- n_spikes + sum(h$outlier)
    val <- .run_extreme(h$amplitude, runs)
    z <- .run_extreme(h$z, runs)
    out[[length(out) + 1]] <- .flag_runs(
      survey, idx, runs, "spikes",
      value = abs(val), threshold = spec$spike_nsigma * abs(val / z), units = "nT",
      description = sprintf("%.1f nT spike (%.0f sigma)", abs(val), abs(z)))
  }
  .with_metric(.bind_flags(out), "spike count", n_spikes, 0)
}

#' Subtract the estimated amplitude of every detected spike
#' @noRd
.despike <- function(x, spec) {
  .spike_detect(x, k = spec$spike_window, nsigma = spec$spike_nsigma)$cleaned
}

# ---- fourth difference -------------------------------------------------------

#' @rdname checks
#' @export
check_fourth_difference <- function(survey, spec) {
  out <- list(); worst <- 0
  series <- rep(NA_real_, nrow(survey))
  for (idx in .by_segment(survey, spec$max_gap)) {
    d4 <- fourth_difference(.despike(survey$mag[idx], spec), normalize = spec$fourth_diff_normalize)
    series[idx] <- d4
    if (any(is.finite(d4))) worst <- max(worst, abs(d4), na.rm = TRUE)
    runs <- .runs(is.finite(d4) & abs(d4) > spec$fourth_diff_tol)
    if (!nrow(runs)) next
    val <- .run_extreme(d4, runs)
    out[[length(out) + 1]] <- .flag_runs(
      survey, idx, runs, "fourth_difference",
      value = abs(val), threshold = spec$fourth_diff_tol, units = "nT",
      description = sprintf("fourth difference %.3f nT (limit %g nT)",
                            abs(val), spec$fourth_diff_tol))
  }
  flags <- .with_metric(.bind_flags(out), "max |4th difference| (nT)", worst, spec$fourth_diff_tol)
  attr(flags, "series") <- series
  flags
}

# ---- diurnal ----------------------------------------------------------------------

#' @rdname checks
#' @export
check_diurnal <- function(base, spec, survey = NULL) {
  stopifnot(all(c("time", "mag_base") %in% names(base)))
  base <- base[order(base$time), ]
  t <- as.numeric(base$time); b <- base$mag_base; n <- length(b)
  empty <- .with_metric(.flag_cols(), "max chord departure (nT)", NA_real_,
                        spec$max_diurnal_dev, note = "base-station record too short")
  if (n < 3) return(empty)
  rate <- stats::median(diff(t))
  w <- max(2L, as.integer(round(spec$diurnal_window / rate)))
  if (n <= w) return(empty)

  dev <- rep(NA_real_, n)
  for (i in seq_len(n - w)) {
    j <- i + w
    if (!is.finite(b[i]) || !is.finite(b[j])) next
    frac <- (t[i:j] - t[i]) / (t[j] - t[i])
    chord <- b[i] + frac * (b[j] - b[i])
    dev[i] <- max(abs(b[i:j] - chord), na.rm = TRUE)
  }
  exceed <- which(is.finite(dev) & dev > spec$max_diurnal_dev)
  affected <- logical(n)
  for (i in exceed) affected[i:(i + w)] <- TRUE
  windows <- .runs(affected)

  stats <- tibble::tibble(time = base$time, mag_base = b, chord_dev = dev, exceed = affected)
  out <- list()
  for (r in seq_len(nrow(windows))) {
    ts <- base$time[windows$start[r]]; te <- base$time[windows$end[r]]
    peak <- max(dev[windows$start[r]:windows$end[r]], na.rm = TRUE)
    desc <- sprintf("acquired during diurnal excursion of %.1f nT over %g s (%s to %s UTC)",
                    peak, spec$diurnal_window,
                    format(ts, "%H:%M:%S"), format(te, "%H:%M:%S"))
    hit <- FALSE
    if (!is.null(survey)) {
      sel <- survey$time >= ts & survey$time <= te
      for (idx in .by_line(survey)) {
        runs <- .runs(sel[idx])
        if (!nrow(runs)) next
        hit <- TRUE
        out[[length(out) + 1]] <- .flag_runs(
          survey, idx, runs, "diurnal", value = peak,
          threshold = spec$max_diurnal_dev, units = "nT", description = desc)
      }
    }
    if (!hit) {
      out[[length(out) + 1]] <- tibble::tibble(
        check = "diurnal", line = NA_character_,
        i_start = NA_integer_, i_end = NA_integer_,
        fid_start = NA_real_, fid_end = NA_real_,
        time_start = ts, time_end = te,
        x = NA_real_, y = NA_real_, lon = NA_real_, lat = NA_real_,
        value = peak, threshold = spec$max_diurnal_dev, units = "nT",
        description = sub("^acquired during", "no samples acquired during", desc))
    }
  }
  flags <- .with_metric(.bind_flags(out), "max chord departure (nT)",
                        if (any(is.finite(dev))) max(dev, na.rm = TRUE) else NA_real_,
                        spec$max_diurnal_dev)
  attr(flags, "base_stats") <- stats
  flags
}

# ---- crossovers ------------------------------------------------------------------

#' @rdname checks
#' @export
check_crossovers <- function(xo, survey, spec) {
  if (!nrow(xo)) {
    return(.with_metric(.flag_cols(), "RMS crossover misfit (nT)", NA_real_,
                        spec$max_crossover_rms, note = "no traverse/tie intersections found"))
  }
  rms <- sqrt(mean(xo$misfit^2, na.rm = TRUE))
  bad <- which(is.finite(xo$misfit) & abs(xo$misfit) > spec$max_crossover_abs)
  flags <- if (length(bad)) {
    tibble::tibble(
      check = "crossovers",
      line = paste(xo$traverse[bad], xo$tie[bad], sep = " x "),
      i_start = xo$i_traverse[bad], i_end = xo$i_traverse[bad],
      fid_start = survey$fid[xo$i_traverse[bad]], fid_end = survey$fid[xo$i_traverse[bad]],
      time_start = xo$time[bad], time_end = xo$time[bad],
      x = xo$x[bad], y = xo$y[bad], lon = xo$lon[bad], lat = xo$lat[bad],
      value = abs(xo$misfit[bad]), threshold = spec$max_crossover_abs, units = "nT",
      description = sprintf("%s crosses %s with %+.1f nT misfit", xo$traverse[bad],
                            xo$tie[bad], xo$misfit[bad]))
  } else .flag_cols()
  .with_metric(flags, "RMS crossover misfit (nT)", rms, spec$max_crossover_rms)
}

# ---- heading error --------------------------------------------------------------

.compass <- function(bearing) {
  c("N", "NE", "E", "SE", "S", "SW", "W", "NW")[1 + (round(bearing / 45) %% 8)]
}

#' @rdname checks
#' @export
check_heading_error <- function(xo, survey, spec) {
  empty <- .with_metric(.flag_cols(), "heading error (nT)", NA_real_,
                        spec$max_heading_error,
                        note = "needs crossovers on traverses flown in both directions")
  if (nrow(xo) < 2) return(empty)
  tr <- survey[survey$line_type == "traverse", ]
  hdg <- dplyr::summarise(dplyr::group_by(tr, .data$line),
                          bearing = .bearing(dplyr::last(.data$x) - dplyr::first(.data$x),
                                             dplyr::last(.data$y) - dplyr::first(.data$y)),
                          .groups = "drop")
  hdg$forward <- abs(((hdg$bearing - spec$line_azimuth + 180) %% 360) - 180) <= 90
  xo$forward <- hdg$forward[match(xo$traverse, hdg$line)]
  by <- dplyr::summarise(dplyr::group_by(xo, .data$forward),
                         n = dplyr::n(), mean_misfit = mean(.data$misfit, na.rm = TRUE),
                         .groups = "drop")
  if (nrow(by) < 2) return(empty)
  fwd <- by[by$forward, ]; rev <- by[!by$forward, ]
  he <- fwd$mean_misfit - rev$mean_misfit
  fwd_lbl <- .compass(mean(hdg$bearing[hdg$forward]))
  rev_lbl <- .compass(mean(hdg$bearing[!hdg$forward]))
  by$heading <- ifelse(by$forward, fwd_lbl, rev_lbl)
  flags <- if (abs(he) > spec$max_heading_error) {
    tibble::tibble(
      check = "heading_error", line = NA_character_,
      i_start = NA_integer_, i_end = NA_integer_,
      fid_start = NA_real_, fid_end = NA_real_,
      time_start = as.POSIXct(NA, tz = "UTC"), time_end = as.POSIXct(NA, tz = "UTC"),
      x = mean(tr$x), y = mean(tr$y), lon = mean(tr$lon), lat = mean(tr$lat),
      value = abs(he), threshold = spec$max_heading_error, units = "nT",
      description = sprintf(
        "mean crossover misfit differs by %.2f nT between %s-bound (%+.2f nT, n=%d) and %s-bound (%+.2f nT, n=%d) traverses",
        abs(he), fwd_lbl, fwd$mean_misfit, fwd$n, rev_lbl, rev$mean_misfit, rev$n))
  } else .flag_cols()
  flags <- .with_metric(flags, "heading error (nT)", abs(he), spec$max_heading_error)
  attr(flags, "by_heading") <- by
  flags
}
