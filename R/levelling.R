#' Tie-line levelling
#'
#' Removes the between-line component of the crossover misfit by fitting a
#' low-order correction to each traverse line and subtracting it, holding
#' the tie lines fixed. With `order = "constant"` each traverse gets a single
#' offset (the mean of its crossover misfits); with `order = "linear"` a
#' traverse with at least `min_crossovers_linear` crossings gets an offset
#' plus a drift term in distance along the line, and the others fall back to
#' a constant. Lines with no crossovers are left untouched.
#'
#' The levelled channel is a separate column, `mag_lev`, so every noise and
#' geometry check keeps operating on the unlevelled, diurnal-corrected field.
#' The heading-error check is likewise computed on pre-levelling misfits -
#' levelling would otherwise absorb exactly the signal it looks for.
#'
#' @param x A `magqc_survey` or a `magqc_result` (whose crossover table is
#'   reused).
#' @param xo Crossover table from [crossovers()]; computed when `NULL`.
#' @param order `"linear"` (default) or `"constant"`.
#' @param min_crossovers_linear Crossings a traverse needs before a drift
#'   term is fitted.
#' @return An object of class `magqc_levelling`: a list with `coefficients`
#'   (one row per traverse: `line`, `order`, `n_crossovers`, `c0` nT, `c1`
#'   nT/m), `crossovers` (the input table plus `along_traverse`,
#'   `misfit_before`, `misfit_after`), `rms` (`c(before, after)`), `mag_lev`
#'   (levelled field aligned with the survey rows) and `order`.
#' @examples
#' res <- run_qc(sim_survey(seed = 2), levelling = FALSE)
#' lev <- level_ties(res)
#' lev$rms
#' @export
level_ties <- function(x, xo = NULL, order = c("linear", "constant"),
                       min_crossovers_linear = 3) {
  order <- match.arg(order)
  if (inherits(x, "magqc_result")) {
    survey <- x$survey
    xo <- xo %||% x$crossovers
  } else if (inherits(x, "magqc_survey")) {
    survey <- x
    xo <- xo %||% crossovers(survey)
  } else {
    stop("`x` must be a magqc_survey or a magqc_result.", call. = FALSE)
  }
  if (!nrow(xo)) {
    stop("No traverse/tie crossovers: levelling needs at least one tie line crossing the traverses.",
         call. = FALSE)
  }

  along <- rep(NA_real_, nrow(survey))
  for (idx in .by_line(survey)) {
    along[idx] <- cumsum(c(0, sqrt(diff(survey$x[idx])^2 + diff(survey$y[idx])^2)))
  }
  xo$along_traverse <- along[xo$i_traverse]
  xo$misfit_before <- xo$misfit

  lines <- unique(survey$line[survey$line_type == "traverse"])
  coef <- lapply(lines, function(ln) {
    j <- which(xo$traverse == ln & is.finite(xo$misfit) & is.finite(xo$along_traverse))
    n <- length(j)
    if (n == 0) return(data.frame(line = ln, order = "none", n_crossovers = 0L, c0 = NA_real_, c1 = NA_real_))
    if (order == "linear" && n >= min_crossovers_linear &&
        diff(range(xo$along_traverse[j])) > 0) {
      fit <- stats::lm.fit(cbind(1, xo$along_traverse[j]), xo$misfit[j])
      b <- fit$coefficients
      if (all(is.finite(b))) {
        return(data.frame(line = ln, order = "linear", n_crossovers = n, c0 = b[[1]], c1 = b[[2]]))
      }
    }
    data.frame(line = ln, order = "constant", n_crossovers = n, c0 = mean(xo$misfit[j]), c1 = 0)
  })
  coef <- tibble::as_tibble(do.call(rbind, coef))

  corr <- rep(0, nrow(survey))
  for (k in which(coef$order != "none")) {
    idx <- which(survey$line == coef$line[k])
    corr[idx] <- coef$c0[k] + coef$c1[k] * along[idx]
  }
  m <- match(xo$traverse, coef$line)
  c0 <- coef$c0[m]; c1 <- coef$c1[m]
  c0[is.na(c0)] <- 0; c1[is.na(c1)] <- 0
  xo$misfit_after <- xo$misfit - (c0 + c1 * xo$along_traverse)

  structure(list(
    coefficients = coef,
    crossovers = xo,
    rms = c(before = sqrt(mean(xo$misfit_before^2, na.rm = TRUE)),
            after = sqrt(mean(xo$misfit_after^2, na.rm = TRUE))),
    mag_lev = survey$mag - corr,
    order = order),
    class = "magqc_levelling")
}

#' @export
print.magqc_levelling <- function(x, ...) {
  cat("<magqc tie-line levelling>\n")
  cat(sprintf("  %s corrections on %d traverses (%d linear, %d constant, %d without crossovers)\n",
              x$order, nrow(x$coefficients), sum(x$coefficients$order == "linear"),
              sum(x$coefficients$order == "constant"), sum(x$coefficients$order == "none")))
  cat(sprintf("  crossover RMS: %.2f nT before, %.2f nT after (%d crossovers)\n",
              x$rms[["before"]], x$rms[["after"]], nrow(x$crossovers)))
  invisible(x)
}

#' @rdname checks
#' @param lev A `magqc_levelling` from [level_ties()], or `NULL` when
#'   levelling was not run.
#' @export
check_levelling_residual <- function(lev, survey, spec) {
  empty <- .with_metric(.flag_cols(), "RMS misfit after levelling (nT)", NA_real_,
                        spec$max_levelled_rms, note = "levelling not run (no crossovers)")
  if (is.null(lev)) return(empty)
  xo <- lev$crossovers
  bad <- which(is.finite(xo$misfit_after) & abs(xo$misfit_after) > spec$max_levelled_abs)
  flags <- if (length(bad)) {
    tibble::tibble(
      check = "levelling_residual",
      line = paste(xo$traverse[bad], xo$tie[bad], sep = " x "),
      i_start = xo$i_traverse[bad], i_end = xo$i_traverse[bad],
      fid_start = survey$fid[xo$i_traverse[bad]], fid_end = survey$fid[xo$i_traverse[bad]],
      time_start = xo$time[bad], time_end = xo$time[bad],
      x = xo$x[bad], y = xo$y[bad], lon = xo$lon[bad], lat = xo$lat[bad],
      value = abs(xo$misfit_after[bad]), threshold = spec$max_levelled_abs, units = "nT",
      description = sprintf("%+.1f nT misfit at %s x %s survives %s levelling (was %+.1f nT)",
                            xo$misfit_after[bad], xo$traverse[bad], xo$tie[bad],
                            lev$coefficients$order[match(xo$traverse[bad], lev$coefficients$line)],
                            xo$misfit_before[bad]))
  } else .flag_cols()
  .with_metric(flags, "RMS misfit after levelling (nT)", unname(lev$rms[["after"]]),
               spec$max_levelled_rms,
               note = sprintf("%.2f nT before levelling", lev$rms[["before"]]))
}
