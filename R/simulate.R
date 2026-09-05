#' Simulate an airborne magnetic survey with known defects
#'
#' Generates a block of traverse and tie lines flown over a synthetic
#' terrain and magnetic field, with a base-station diurnal record, and - when
#' `defects = TRUE` - injects a catalogue of the failures a real acceptance
#' review is looking for. Every injected defect is recorded in the returned
#' `truth` table, which is what the package tests assert against and what
#' lets a reader of the report check that each flagged interval is a real
#' problem rather than a false positive.
#'
#' Injected defects (all on traverse lines unless stated):
#' * isolated spikes of 5-40 nT on random samples
#' * a stretch where sensor noise rises to 0.25 nT (fourth difference)
#' * a stretch flown 45 m off the planned line
#' * one whole line flown 60 m off position (line separation)
#' * a stretch flown 35 m too high (terrain clearance)
#' * a 12 s dropout in the data stream (gap)
#' * a stretch with every second sample missing (along-line spacing)
#' * a magnetic-storm excursion on the base station (diurnal)
#' * a +2.5 nT bias on southbound lines (heading error, via crossovers)
#'
#' Defect targets are chosen from the generated geometry (a traverse by its
#' position in the block, a window by fraction of the line or by duration),
#' so any `n_lines >= 2` and any line length work; a defect that cannot be
#' placed on a very short line is skipped with a warning and gets no `truth`
#' row. With the default arguments the targets are L1050 (noise), L1090
#' (off-plan), L1110 (mispositioned), L1030 (clearance), L1070 (dropout) and
#' L1010 (sparse).
#'
#' @param n_lines Number of traverse lines (at least 2).
#' @param line_spacing Traverse spacing, m.
#' @param line_length Traverse length, m.
#' @param tie_spacing Tie-line spacing, m.
#' @param sample_rate Magnetometer sample rate, Hz.
#' @param ground_speed Ground speed, m/s.
#' @param nominal_clearance Terrain clearance, m.
#' @param origin `c(lon, lat)` of the local grid origin.
#' @param start_time POSIXct start of the flight.
#' @param inclination,declination Ambient field direction the anomalies are
#'   induced in, degrees. The defaults match the IGRF at the default origin
#'   (central Wyoming, 2024), so the anomalies are asymmetric the way real
#'   ones are at that latitude and [reduce_to_pole()] straightens them.
#' @param defects Inject defects?
#' @param seed Random seed.
#'
#' @return A list of class `magqc_sim` with elements `data` (an
#'   [as_survey()] tibble), `base` (base-station tibble: `time`,
#'   `mag_base`), `truth` (tibble of injected defects with the time window
#'   each one occupies), `spec` (a [survey_spec()] matching the flight
#'   parameters) and `origin`.
#' @examples
#' sim <- sim_survey(seed = 1)
#' nrow(sim$data)
#' sim$truth
#' @export
sim_survey <- function(n_lines = 14, line_spacing = 200, line_length = 7000,
                       tie_spacing = 2000, sample_rate = 10, ground_speed = 60,
                       nominal_clearance = 60,
                       origin = c(lon = -108.5, lat = 43.2),
                       start_time = as.POSIXct("2024-06-14 08:00:00", tz = "UTC"),
                       inclination = 68, declination = 10,
                       defects = TRUE, seed = 42) {
  if (!is.numeric(n_lines) || n_lines < 2) stop("`n_lines` must be at least 2.", call. = FALSE)
  if (line_length <= 0 || tie_spacing <= 0 || line_spacing <= 0 || sample_rate <= 0 || ground_speed <= 0) {
    stop("`line_length`, `line_spacing`, `tie_spacing`, `sample_rate` and `ground_speed` must be positive.",
         call. = FALSE)
  }
  set.seed(seed)
  dt <- 1 / sample_rate
  step <- ground_speed * dt
  turn_time <- 90

  # ---- flight geometry -----------------------------------------------------
  lines <- list()
  along_tr <- seq(0, line_length, by = step)
  for (i in seq_len(n_lines) - 1) {
    x0 <- i * line_spacing
    northbound <- i %% 2 == 0
    y <- if (northbound) along_tr else rev(along_tr)
    lines[[length(lines) + 1]] <- tibble::tibble(
      line = sprintf("L%d", 1000 + i * 10), line_type = "traverse",
      x = x0 + .smooth_noise(length(y), sd = 4, span = 200),
      y = y, dir = if (northbound) "N" else "S")
  }
  x_min <- -100; x_max <- (n_lines - 1) * line_spacing + 100
  along_tie <- seq(x_min, x_max, by = step)
  # ties 1 km in from each end at the tie spacing; a block too short for that
  # gets a single tie across the middle
  tie_y <- if (line_length >= 2000) seq(1000, line_length - 1000, by = tie_spacing) else line_length / 2
  for (j in seq_along(tie_y)) {
    eastbound <- j %% 2 == 1
    x <- if (eastbound) along_tie else rev(along_tie)
    lines[[length(lines) + 1]] <- tibble::tibble(
      line = sprintf("T%d", 9000 + (j - 1) * 10), line_type = "tie",
      x = x, y = tie_y[j] + .smooth_noise(length(x), sd = 4, span = 200),
      dir = if (eastbound) "E" else "W")
  }

  # ---- time base ------------------------------------------------------------
  t_cursor <- 0
  for (k in seq_along(lines)) {
    n <- nrow(lines[[k]])
    lines[[k]]$t <- t_cursor + seq(0, by = dt, length.out = n)
    t_cursor <- t_cursor + (n - 1) * dt + turn_time
  }
  d <- dplyr::bind_rows(lines)

  # ---- terrain and altitude ---------------------------------------------
  dem <- function(x, y) 1500 + 40 * sin(x / 900) + 30 * cos(y / 1300) + 15 * sin((x + y) / 500)
  d$dem <- dem(d$x, d$y)
  d$radar_alt <- nominal_clearance + .smooth_noise(nrow(d), sd = 3, span = 150)

  # ---- magnetic field -------------------------------------------------------
  # Five buried dipoles magnetised along the ambient field (the total-field
  # anomaly of each is the real inclined-field shape, not a symmetric bump),
  # observed at the sensor height above ground. `amp` is the peak a vertical
  # field would give (2 m / h^3).
  sources <- tibble::tibble(
    xs = c(600, 1500, 2100, 1000, 2400),
    ys = c(1500, 3800, 5600, 6200, 2600),
    depth = c(250, 400, 180, 600, 320),
    amp = c(180, 320, 90, 220, 140))
  anom <- rep(0, nrow(d))
  for (s in seq_len(nrow(sources))) {
    h <- sources$depth[s] + nominal_clearance
    anom <- anom + .dipole_anomaly(d$x, d$y, sources$xs[s], sources$ys[s], h,
                                   inclination, declination, moment = sources$amp[s] * h^3 / 2)
  }
  regional <- 54000 + 0.0025 * d$y + 0.0010 * d$x
  noise_sd <- rep(0.015, nrow(d))

  # ---- base station ---------------------------------------------------------
  t_base <- seq(0, max(d$t) + 60, by = 1)
  diurnal <- 25 * sin(2 * pi * (t_base + 8 * 3600) / 86400 - 1.2) +
    .smooth_noise(length(t_base), sd = 0.4, span = 121)
  base <- tibble::tibble(time = start_time + t_base, mag_base = 54000 + diurnal)

  # ---- defects --------------------------------------------------------------
  truth <- tibble::tibble(defect = character(), check = character(),
                          line = character(),
                          time_start = as.POSIXct(character(), tz = "UTC"),
                          time_end = as.POSIXct(character(), tz = "UTC"),
                          detail = character())
  add_truth <- function(defect, check, line, rows, detail) {
    ts <- if (length(rows)) start_time + min(d$t[rows]) else as.POSIXct(NA, tz = "UTC")
    te <- if (length(rows)) start_time + max(d$t[rows]) else as.POSIXct(NA, tz = "UTC")
    truth <<- dplyr::bind_rows(truth, tibble::tibble(
      defect = defect, check = check, line = line,
      time_start = ts, time_end = te, detail = detail))
  }
  heading_bias <- rep(0, nrow(d))
  drop_rows <- integer(0)
  spike <- rep(0, nrow(d))

  if (isTRUE(defects)) {
    # Defect targets are derived from the geometry that was actually
    # generated - a traverse chosen by its position in the block, a window
    # chosen as a fraction of that line's samples or as a duration - so the
    # simulator stays truthful for any n_lines / line_length / sample_rate
    # (#3, #4). With the defaults these resolve to the same lines and nearly
    # the same samples as the original hard-coded targets.
    tr_lines <- unique(d$line[d$line_type == "traverse"])
    targets <- .pick_lines(tr_lines, c(noisy = 0.43, off_plan = 0.71, misposition = 0.86,
                                       clearance = 0.29, dropout = 0.57, sparse = 0.14))
    defect_rows <- function(defect, line, f0 = 0, f1 = 1, min_samples = 10L) {
      rows <- which(d$line == line)
      n <- length(rows)
      if (n < min_samples) {
        warning(sprintf("sim_survey(): skipping defect '%s' - line %s has %d samples, needs at least %d",
                        defect, line, n, min_samples), call. = FALSE)
        return(integer(0))
      }
      rows[max(1L, round(f0 * n)):min(n, max(1L, round(f1 * n)))]
    }

    # heading error: southbound lines read high
    heading_bias[d$dir == "S"] <- 2.5
    add_truth("heading bias", "heading_error", NA_character_, integer(0),
              "+2.5 nT on southbound traverses")

    # noisy segment: sensor noise raised over the middle of the line
    ln <- targets[["noisy"]]
    r <- defect_rows("noisy segment", ln, 0.34, 0.60, min_samples = 60L)
    if (length(r)) {
      noise_sd[r] <- 0.25
      add_truth("noisy segment", "fourth_difference", ln, r, "sensor noise 0.25 nT")
    }

    # off-plan stretch (smooth 45 m excursion)
    ln <- targets[["off_plan"]]
    r <- defect_rows("off-plan excursion", ln, 0.26, 0.69, min_samples = 60L)
    if (length(r)) {
      d$x[r] <- d$x[r] + 45 * sin(seq(0, pi, length.out = length(r)))
      add_truth("off-plan excursion", "line_deviation", ln, r, "45 m east of plan")
    }

    # whole line mispositioned
    ln <- targets[["misposition"]]
    r <- defect_rows("line mispositioned", ln)
    if (length(r)) {
      d$x[r] <- d$x[r] + 60
      add_truth("line mispositioned", "line_separation", ln, r, "flown 60 m east")
    }

    # clearance excursion
    ln <- targets[["clearance"]]
    r <- defect_rows("clearance excursion", ln, 0.17, 0.43, min_samples = 60L)
    if (length(r)) {
      d$radar_alt[r] <- d$radar_alt[r] + 35 * sin(seq(0, pi, length.out = length(r)))
      add_truth("clearance excursion", "clearance", ln, r, "up to 35 m above nominal")
    }

    # dropout: 12 s of records removed from the middle of the line
    ln <- targets[["dropout"]]
    n_drop <- as.integer(round(12 * sample_rate))
    r <- defect_rows("data dropout", ln, 0.51, 0.51, min_samples = 4L * n_drop)
    if (length(r)) {
      r <- r[1] + seq_len(n_drop) - 1L
      drop_rows <- c(drop_rows, r)
      add_truth("data dropout", "sample_interval", ln, r, "12 s of records missing")
    }

    # sparse sampling: every second sample missing over ~10 s
    ln <- targets[["sparse"]]
    n_sparse <- as.integer(round(10 * sample_rate))
    r <- defect_rows("sparse sampling", ln, 0.69, 0.69, min_samples = 4L * n_sparse)
    if (length(r)) {
      r <- r[1] + seq_len(n_sparse) - 1L
      drop_rows <- c(drop_rows, r[seq(2, length(r), by = 2)])
      add_truth("sparse sampling", "along_line_spacing", ln, r, "every second sample missing")
    }

    # spikes
    tr <- setdiff(which(d$line_type == "traverse"), drop_rows)
    n_spikes <- min(25L, length(tr))
    sp <- sort(sample(tr, n_spikes))
    spike[sp] <- sample(c(-1, 1), n_spikes, TRUE) * stats::runif(n_spikes, 5, 40)
    for (k in sp) add_truth("spike", "spikes", d$line[k], k, sprintf("%.1f nT", spike[k]))

    # diurnal storm on the base station, six minutes starting 40% of the way
    # through the flight
    storm_t0 <- round(0.4 * max(d$t)); storm_len <- 6 * 60
    st <- which(t_base >= storm_t0 & t_base <= storm_t0 + storm_len)
    tt <- (t_base[st] - t_base[st[1]]) / 60
    base$mag_base[st] <- base$mag_base[st] +
      12 * sin(2 * pi * tt / 2) * exp(-((tt - 3) / 1.5)^2)
    add_truth("diurnal storm", "diurnal", NA_character_, integer(0),
              sprintf("+/-12 nT, 2 min period, %s to %s UTC",
                      format(start_time + storm_t0, "%H:%M"),
                      format(start_time + storm_t0 + storm_len, "%H:%M")))
    truth$time_start[truth$check == "diurnal"] <- start_time + storm_t0
    truth$time_end[truth$check == "diurnal"] <- start_time + storm_t0 + storm_len
  }

  d$gps_alt <- d$dem + d$radar_alt
  d$mag_raw <- regional + anom + heading_bias + spike +
    stats::rnorm(nrow(d), 0, noise_sd) +
    stats::approx(t_base, base$mag_base - 54000, xout = d$t)$y
  d$time <- start_time + d$t

  # never let an NA reach the negative subscript (#3)
  drop_rows <- unique(drop_rows[!is.na(drop_rows)])
  if (length(drop_rows)) d <- d[-drop_rows, ]
  d$fid <- seq_len(nrow(d))
  d$dir <- NULL; d$t <- NULL; d$dem <- NULL
  ll <- enu_to_lonlat(d$x, d$y, origin[["lon"]], origin[["lat"]])
  d$lon <- ll$lon; d$lat <- ll$lat

  spec <- survey_spec(line_spacing = line_spacing, tie_spacing = tie_spacing,
                      sample_rate = sample_rate, ground_speed = ground_speed,
                      nominal_clearance = nominal_clearance)
  structure(list(data = as_survey(d), base = base, truth = truth, spec = spec,
                 origin = origin),
            class = "magqc_sim")
}

#' Pick distinct traverse lines by relative position in the block
#'
#' `fracs` are positions in (0, 1]; each maps to the nearest line, and when
#' two defects would land on the same line the later one moves to the
#' nearest unused line so the tests can still attribute each flag to one
#' defect. With fewer lines than defects the leftovers share.
#' @noRd
.pick_lines <- function(lines, fracs) {
  n <- length(lines)
  idx <- pmin(n, pmax(1L, as.integer(round(fracs * n))))
  used <- logical(n)
  for (k in seq_along(idx)) {
    if (!used[idx[k]]) { used[idx[k]] <- TRUE; next }
    free <- which(!used)
    if (!length(free)) next
    idx[k] <- free[which.min(abs(free - idx[k]))]
    used[idx[k]] <- TRUE
  }
  stats::setNames(lines[idx], names(fracs))
}

#' Smooth, zero-mean Gaussian noise
#'
#' White noise convolved with a Gaussian kernel whose standard deviation is
#' `span / 4` samples, then rescaled so the output has standard deviation
#' `sd`. Used for GPS wander, altitude hunting and slow diurnal texture.
#'
#' The kernel matters: a boxcar-smoothed series looks smooth but its
#' *increments* are white, so lateral wander of a few metres crossed with an
#' anomaly gradient of ~1 nT/m produced sample-to-sample jitter that tripped
#' the fourth-difference test on a supposedly clean survey. A Gaussian kernel
#' is smooth in every derivative, like real aircraft motion.
#' @noRd
.smooth_noise <- function(n, sd = 1, span = 100) {
  if (n < 2) return(rep(0, n))
  ksd <- max(1, min(span, n - 1) / 4)
  half <- ceiling(3 * ksd)
  w <- stats::dnorm(-half:half, sd = ksd); w <- w / sum(w)
  x <- stats::rnorm(n + 2 * half)
  s <- stats::filter(x, w, sides = 2)
  s <- as.numeric(s)[(half + 1):(half + n)]
  s <- s - mean(s)
  sdev <- stats::sd(s)
  if (!is.finite(sdev) || sdev == 0) return(rep(0, n))
  s * sd / sdev
}

#' @export
print.magqc_sim <- function(x, ...) {
  cat("<magqc simulated survey>\n")
  cat(sprintf("  %d samples on %d lines (%d traverse, %d tie)\n",
              nrow(x$data), length(unique(x$data$line)),
              length(unique(x$data$line[x$data$line_type == "traverse"])),
              length(unique(x$data$line[x$data$line_type == "tie"]))))
  cat(sprintf("  %d injected defects\n", nrow(x$truth)))
  invisible(x)
}
