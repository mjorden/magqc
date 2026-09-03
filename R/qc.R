.check_registry <- function() {
  tibble::tibble(
    check = c("sample_interval", "along_line_spacing", "line_deviation",
              "line_separation", "clearance", "spikes", "fourth_difference",
              "diurnal", "crossovers", "heading_error", "levelling_residual"),
    label = c("Data gaps", "Along-line sampling", "Departure from line",
              "Line separation", "Terrain clearance", "Spikes",
              "Fourth-difference noise", "Diurnal variation",
              "Crossover misfit", "Heading error", "Post-levelling residual"))
}

#' Run every QC check on a survey
#'
#' Applies the diurnal correction (when a base-station record is supplied),
#' computes the crossover table, runs each `check_*()` from the
#' [checks] family and assembles the flags, a per-check scorecard, per-line
#' statistics and survey totals into a `magqc_result` that [qc_report()]
#' renders.
#'
#' @param survey A `magqc_survey` from [as_survey()], a plain data frame that
#'   [as_survey()] accepts, or a `magqc_sim` from [sim_survey()] (in which case
#'   `base` and `spec` default to the simulation's own).
#' @param spec A [survey_spec()].
#' @param base Optional base-station record (`time`, `mag_base`). When given,
#'   `mag` is set to `mag_raw` minus the base-station departure from its
#'   median level before any field statistic is computed.
#' @param plan Optional flight plan; see [check_line_deviation()].
#' @param levelling Run tie-line levelling ([level_ties()]) after the checks
#'   and report the post-levelling residual on the scorecard? Skipped
#'   silently when there are no crossovers.
#' @param levelling_order Passed to [level_ties()] as `order`.
#' @return A list of class `magqc_result` with elements `survey` (with a
#'   `mag_lev` column when levelling ran), `base`, `spec`, `flags`,
#'   `scorecard`, `lines`, `stats`, `crossovers`, `diurnal` (base-station
#'   statistics), `fourth_difference` (per-sample series), `heading` (misfit
#'   by flight direction) and `levelling` (a `magqc_levelling`, or `NULL`).
#' @examples
#' res <- run_qc(sim_survey(seed = 3))
#' res$scorecard
#' res$levelling
#' @export
run_qc <- function(survey, spec = survey_spec(), base = NULL, plan = NULL,
                   levelling = TRUE, levelling_order = c("linear", "constant")) {
  if (inherits(survey, "magqc_sim")) {
    if (missing(spec)) spec <- survey$spec
    base <- base %||% survey$base
    survey <- survey$data
  }
  if (!inherits(survey, "magqc_survey")) survey <- as_survey(survey)
  if (!inherits(spec, "magqc_spec")) stop("`spec` must come from survey_spec().", call. = FALSE)

  if (!is.null(base)) {
    base <- base[order(base$time), ]
    datum <- stats::median(base$mag_base, na.rm = TRUE)
    # ties = mean: duplicate base timestamps are averaged rather than warned about (#2)
    corr <- stats::approx(as.numeric(base$time), base$mag_base, xout = as.numeric(survey$time),
                          rule = 2, ties = mean)$y
    survey$mag <- survey$mag_raw - (corr - datum)
  }

  xo <- crossovers(survey)
  lev <- NULL
  if (isTRUE(levelling) && nrow(xo)) {
    lev <- level_ties(survey, xo, order = match.arg(levelling_order))
    survey$mag_lev <- lev$mag_lev
  }
  checks <- list(
    sample_interval    = check_sample_interval(survey, spec),
    along_line_spacing = check_along_line_spacing(survey, spec),
    line_deviation     = check_line_deviation(survey, spec, plan),
    line_separation    = check_line_separation(survey, spec),
    clearance          = check_clearance(survey, spec),
    spikes             = check_spikes(survey, spec),
    fourth_difference  = check_fourth_difference(survey, spec),
    diurnal            = if (!is.null(base)) check_diurnal(base, spec, survey) else
      .with_metric(.flag_cols(), "max chord departure (nT)", NA_real_,
                   spec$max_diurnal_dev, note = "no base-station record supplied"),
    crossovers         = check_crossovers(xo, survey, spec),
    heading_error      = check_heading_error(xo, survey, spec),
    levelling_residual = check_levelling_residual(lev, survey, spec))

  flags <- .bind_flags(unname(checks))
  reg <- .check_registry()

  step <- rep(NA_real_, nrow(survey))
  for (idx in .by_line(survey)) step[idx] <- .step_distance(survey$x[idx], survey$y[idx])
  step[!is.finite(step)] <- 0
  flagged <- logical(nrow(survey))
  per_check_samples <- integer(nrow(reg))
  for (k in seq_len(nrow(reg))) {
    f <- checks[[reg$check[k]]]
    m <- logical(nrow(survey))
    for (r in which(!is.na(f$i_start))) m[f$i_start[r]:f$i_end[r]] <- TRUE
    per_check_samples[k] <- sum(m)
    flagged <- flagged | m
  }

  scorecard <- reg
  # index by name, never by position (#9): the registry order is the report
  # order and must not have to match the order the checks were run in
  stopifnot(setequal(names(checks), reg$check))
  metric_of <- function(chk, field) attr(checks[[chk]], "metric")[[field]]
  by_check <- function(f, type) unname(vapply(reg$check, f, type))
  scorecard$metric_label <- by_check(function(chk) metric_of(chk, "label"), character(1))
  scorecard$metric <- by_check(function(chk) as.numeric(metric_of(chk, "value")), numeric(1))
  scorecard$threshold <- by_check(function(chk) as.numeric(metric_of(chk, "threshold")), numeric(1))
  scorecard$note <- by_check(function(chk) metric_of(chk, "note") %||% NA_character_, character(1))
  scorecard$n_flags <- by_check(function(chk) nrow(checks[[chk]]), integer(1))
  scorecard$n_samples <- per_check_samples
  scorecard$pct_samples <- 100 * per_check_samples / nrow(survey)
  scorecard$status <- ifelse(
    is.na(scorecard$metric), "n/a",
    ifelse(scorecard$metric > scorecard$threshold | scorecard$n_flags > 0, "fail", "pass"))

  survey$.step <- step
  survey$.flagged <- flagged
  lines <- dplyr::summarise(
    dplyr::group_by(survey, .data$line, .data$line_type),
    n_samples = dplyr::n(),
    start = min(.data$time), end = max(.data$time),
    length_km = sum(.data$.step) / 1000,
    flagged_km = sum(.data$.step[.data$.flagged]) / 1000,
    n_flags = 0L, .groups = "drop")
  fl_by_line <- table(flags$line[!is.na(flags$i_start)])
  lines$n_flags <- as.integer(fl_by_line[lines$line]); lines$n_flags[is.na(lines$n_flags)] <- 0L
  lines$status <- ifelse(lines$n_flags > 0, "fail", "pass")
  lines <- lines[order(lines$line_type == "tie", lines$line), ]
  survey$.step <- NULL; survey$.flagged <- NULL

  total_km <- sum(lines$length_km)
  stats <- list(
    n_samples = nrow(survey),
    n_lines = nrow(lines),
    n_traverse = sum(lines$line_type == "traverse"),
    n_tie = sum(lines$line_type == "tie"),
    line_km = total_km,
    traverse_km = sum(lines$length_km[lines$line_type == "traverse"]),
    tie_km = sum(lines$length_km[lines$line_type == "tie"]),
    flagged_km = sum(lines$flagged_km),
    accepted_pct = if (total_km > 0) 100 * (1 - sum(lines$flagged_km) / total_km) else NA_real_,
    start = min(survey$time), end = max(survey$time),
    duration_h = as.numeric(difftime(max(survey$time), min(survey$time), units = "hours")),
    n_crossovers = nrow(xo),
    n_flags = nrow(flags),
    checks_passed = sum(scorecard$status == "pass"),
    checks_failed = sum(scorecard$status == "fail"),
    checks_total = sum(scorecard$status != "n/a"),
    has_lonlat = any(is.finite(survey$lon)),
    diurnal_corrected = !is.null(base))

  structure(list(
    survey = survey, base = base, spec = spec, flags = flags,
    scorecard = scorecard, lines = lines, stats = stats, crossovers = xo,
    diurnal = attr(checks$diurnal, "base_stats"),
    fourth_difference = attr(checks$fourth_difference, "series"),
    heading = attr(checks$heading_error, "by_heading"),
    levelling = lev),
    class = "magqc_result")
}

#' @export
print.magqc_result <- function(x, ...) {
  s <- x$stats
  cat("<magqc result>\n")
  cat(sprintf("  %s samples, %d lines, %.1f line-km, %.1f h of flying\n",
              format(s$n_samples, big.mark = ","), s$n_lines, s$line_km, s$duration_h))
  cat(sprintf("  %d of %d checks passed; %.1f line-km flagged (%.1f%% accepted)\n\n",
              s$checks_passed, s$checks_total, s$flagged_km, s$accepted_pct))
  sc <- x$scorecard
  for (k in seq_len(nrow(sc))) {
    cat(sprintf("  %-5s %-24s %-32s %s\n",
                toupper(sc$status[k]), sc$label[k],
                sprintf("%s = %s", sc$metric_label[k],
                        if (is.na(sc$metric[k])) "n/a" else format(signif(sc$metric[k], 3))),
                if (sc$n_flags[k] > 0) sprintf("%d flag(s)", sc$n_flags[k]) else ""))
  }
  invisible(x)
}
