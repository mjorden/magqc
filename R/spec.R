#' Define a survey acceptance specification
#'
#' Every QC check in this package reads its threshold from a `magqc_spec`
#' object. The defaults describe a typical fixed-wing high-resolution
#' magnetic survey (200 m traverse spacing, 60 m nominal terrain clearance,
#' 10 Hz acquisition); they are a starting point, not a standard. Real
#' acceptance criteria come from the contract for the survey being flown, so
#' override anything that does not match yours.
#'
#' @param line_spacing Nominal separation between adjacent traverse lines, m.
#' @param tie_spacing Nominal separation between tie (control) lines, m.
#' @param line_azimuth Planned traverse heading, degrees clockwise from north.
#' @param tie_azimuth Planned tie-line heading, degrees clockwise from north.
#' @param sample_rate Acquisition rate of the magnetometer, Hz.
#' @param ground_speed Nominal survey ground speed, m/s.
#' @param max_along_line_sep Maximum accepted ground distance between
#'   consecutive samples, m. Defaults to 1.5x the nominal
#'   `ground_speed / sample_rate`, which tolerates a modest tailwind.
#' @param line_spacing_tol Fractional tolerance on `line_spacing`; 0.10 accepts
#'   an actual separation within +/-10% of nominal.
#' @param max_deviation Maximum lateral departure from the planned (or, absent
#'   a plan, the best-fit) line, m.
#' @param nominal_clearance Planned terrain clearance (radar altimeter), m.
#' @param clearance_tol Accepted departure from `nominal_clearance`, m.
#' @param fourth_diff_tol Fourth-difference noise envelope, nT. Interpreted on
#'   the normalised operator - see [fourth_difference()].
#' @param fourth_diff_normalize Whether the fourth-difference operator is
#'   divided by 16 before comparison with `fourth_diff_tol`.
#' @param spike_nsigma Robust z-score of the fourth difference above which a
#'   sample is called a spike.
#' @param spike_window Window length, in samples, for the running MAD that
#'   scales the fourth difference in the spike test. Must be odd.
#' @param max_diurnal_dev Maximum departure of the base-station field from a
#'   linear chord fitted over `diurnal_window`, nT. This is the usual form of
#'   the diurnal specification: it constrains curvature and short-period
#'   activity rather than the total drift over the flight.
#' @param diurnal_window Chord length for the diurnal test, seconds.
#' @param max_crossover_rms Maximum RMS misfit at traverse/tie intersections
#'   before levelling, nT.
#' @param max_crossover_abs Per-intersection misfit above which a single
#'   crossover is flagged, nT.
#' @param max_heading_error Maximum difference in mean crossover misfit between
#'   opposing flight directions, nT.
#' @param max_gap Longest accepted break in the data stream, seconds.
#'
#' @return An object of class `magqc_spec`.
#' @examples
#' spec <- survey_spec(line_spacing = 100, nominal_clearance = 40)
#' spec$max_along_line_sep
#' @export
survey_spec <- function(line_spacing = 200,
                        tie_spacing = 2000,
                        line_azimuth = 0,
                        tie_azimuth = 90,
                        sample_rate = 10,
                        ground_speed = 60,
                        max_along_line_sep = NULL,
                        line_spacing_tol = 0.10,
                        max_deviation = 25,
                        nominal_clearance = 60,
                        clearance_tol = 15,
                        fourth_diff_tol = 0.05,
                        fourth_diff_normalize = TRUE,
                        spike_nsigma = 6,
                        spike_window = 21,
                        max_diurnal_dev = 3.0,
                        diurnal_window = 60,
                        max_crossover_rms = 2.0,
                        max_crossover_abs = 6.0,
                        max_heading_error = 2.0,
                        max_gap = 1.0) {

  .assert_pos <- function(v, nm) {
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0) {
      stop("`", nm, "` must be a single positive number.", call. = FALSE)
    }
  }
  for (nm in c("line_spacing", "tie_spacing", "sample_rate", "ground_speed",
               "line_spacing_tol", "max_deviation", "nominal_clearance",
               "clearance_tol", "fourth_diff_tol", "spike_nsigma",
               "spike_window", "max_diurnal_dev", "diurnal_window",
               "max_crossover_rms", "max_crossover_abs", "max_heading_error",
               "max_gap")) {
    .assert_pos(get(nm), nm)
  }
  if (spike_window %% 2 == 0) {
    stop("`spike_window` must be odd so the running median is centred.", call. = FALSE)
  }
  if (is.null(max_along_line_sep)) {
    max_along_line_sep <- 1.5 * ground_speed / sample_rate
  }
  .assert_pos(max_along_line_sep, "max_along_line_sep")

  structure(
    list(line_spacing = line_spacing,
         tie_spacing = tie_spacing,
         line_azimuth = line_azimuth %% 360,
         tie_azimuth = tie_azimuth %% 360,
         sample_rate = sample_rate,
         ground_speed = ground_speed,
         max_along_line_sep = max_along_line_sep,
         line_spacing_tol = line_spacing_tol,
         max_deviation = max_deviation,
         nominal_clearance = nominal_clearance,
         clearance_tol = clearance_tol,
         fourth_diff_tol = fourth_diff_tol,
         fourth_diff_normalize = fourth_diff_normalize,
         spike_nsigma = spike_nsigma,
         spike_window = spike_window,
         max_diurnal_dev = max_diurnal_dev,
         diurnal_window = diurnal_window,
         max_crossover_rms = max_crossover_rms,
         max_crossover_abs = max_crossover_abs,
         max_heading_error = max_heading_error,
         max_gap = max_gap),
    class = "magqc_spec"
  )
}

#' @export
print.magqc_spec <- function(x, ...) {
  cat("<magqc survey specification>\n")
  fmt <- function(label, value, unit = "") {
    cat(sprintf("  %-24s %s%s\n", label, format(value), unit))
  }
  fmt("traverse spacing", x$line_spacing, " m")
  fmt("tie spacing", x$tie_spacing, " m")
  fmt("traverse azimuth", x$line_azimuth, " deg")
  fmt("sample rate", x$sample_rate, " Hz")
  fmt("max sample spacing", round(x$max_along_line_sep, 2), " m")
  fmt("max line departure", x$max_deviation, " m")
  fmt("terrain clearance", sprintf("%g +/- %g", x$nominal_clearance, x$clearance_tol), " m")
  fmt("4th-difference limit", x$fourth_diff_tol, " nT")
  fmt("diurnal limit", sprintf("%g nT over %g", x$max_diurnal_dev, x$diurnal_window), " s")
  fmt("crossover RMS limit", x$max_crossover_rms, " nT")
  fmt("heading error limit", x$max_heading_error, " nT")
  fmt("max data gap", x$max_gap, " s")
  invisible(x)
}
