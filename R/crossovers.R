#' Traverse / tie-line crossover table
#'
#' Finds every intersection between a traverse line and a tie line and
#' interpolates the diurnal-corrected field on each line at the crossing
#' point. The misfit `mag_traverse - mag_tie` is the raw material for the
#' crossover and heading-error checks and, later, for levelling.
#'
#' Each flown line is first reduced to its total-least-squares axis and the
#' two axes intersected to find a candidate crossing. The candidate is then
#' refined by re-fitting each line locally (about ten samples either side of
#' the nearest sample) and intersecting the local axes, so a line that wanders
#' off its overall axis still crosses where it actually crosses. The crossing
#' is kept if it falls inside both lines' extents and within `max_dist` of a
#' real sample on each; the field is then linearly interpolated between the
#' bracketing samples along each line.
#'
#' @param survey A `magqc_survey`.
#' @param max_dist Reject intersections further than this from the nearest
#'   sample on either line, m. Defaults to three sample spacings.
#' @return A tibble with one row per crossover: `traverse`, `tie`, `x`, `y`,
#'   `lon`, `lat`, `time`, `mag_traverse`, `mag_tie`, `misfit`,
#'   `alt_traverse`, `alt_tie`, `i_traverse`, `i_tie`.
#' @export
crossovers <- function(survey, max_dist = NULL) {
  empty <- tibble::tibble(
    traverse = character(), tie = character(), x = numeric(), y = numeric(),
    lon = numeric(), lat = numeric(), time = as.POSIXct(character(), tz = "UTC"),
    mag_traverse = numeric(), mag_tie = numeric(), misfit = numeric(),
    alt_traverse = numeric(), alt_tie = numeric(),
    i_traverse = integer(), i_tie = integer())
  tr_lines <- unique(survey$line[survey$line_type == "traverse"])
  ti_lines <- unique(survey$line[survey$line_type == "tie"])
  if (!length(tr_lines) || !length(ti_lines)) return(empty)

  fits <- lapply(split(survey, survey$line), function(s) {
    if (nrow(s) < 3) return(NULL)
    f <- .tls_line(s$x, s$y)
    list(s = s, fit = f, bbox = c(range(s$x), range(s$y)),
         step = stats::median(.step_distance(s$x, s$y), na.rm = TRUE))
  })
  local_fit <- function(L, k, half = 10) {
    i <- max(1, k - half):min(nrow(L$s), k + half)
    list(i = i, fit = .tls_line(L$s$x[i], L$s$y[i]))
  }
  nearest <- function(L, P) {
    d <- sqrt((L$s$x - P[1])^2 + (L$s$y - P[2])^2)
    k <- which.min(d)
    list(k = k, dist = d[k])
  }
  interp <- function(L, lf, P, col) {
    v <- L$s[[col]][lf$i]
    if (sum(is.finite(v)) < 2) return(NA_real_)
    aP <- (P[1] - lf$fit$cx) * lf$fit$dir[1] + (P[2] - lf$fit$cy) * lf$fit$dir[2]
    stats::approx(lf$fit$along, v, xout = aP, ties = mean, rule = 2)$y
  }
  inside <- function(P, bbox, pad) {
    P[1] >= bbox[1] - pad && P[1] <= bbox[2] + pad &&
      P[2] >= bbox[3] - pad && P[2] <= bbox[4] + pad
  }

  out <- list()
  for (a in tr_lines) for (b in ti_lines) {
    A <- fits[[a]]; B <- fits[[b]]
    if (is.null(A) || is.null(B)) next
    P <- .line_intersection(c(A$fit$cx, A$fit$cy), A$fit$dir, c(B$fit$cx, B$fit$cy), B$fit$dir)
    if (anyNA(P)) next
    md <- max_dist %||% (3 * max(A$step, B$step, na.rm = TRUE))
    if (!inside(P, A$bbox, md) || !inside(P, B$bbox, md)) next
    nA <- nearest(A, P); nB <- nearest(B, P)
    # refine: intersect the local axes around the nearest samples, twice
    for (it in 1:2) {
      fa <- local_fit(A, nA$k); fb <- local_fit(B, nB$k)
      P2 <- .line_intersection(c(fa$fit$cx, fa$fit$cy), fa$fit$dir, c(fb$fit$cx, fb$fit$cy), fb$fit$dir)
      if (anyNA(P2)) break
      P <- P2
      nA <- nearest(A, P); nB <- nearest(B, P)
    }
    if (nA$dist > md || nB$dist > md) next
    fa <- local_fit(A, nA$k); fb <- local_fit(B, nB$k)
    kA <- nA$k; kB <- nB$k
    out[[length(out) + 1]] <- tibble::tibble(
      traverse = a, tie = b, x = P[1], y = P[2],
      lon = A$s$lon[kA], lat = A$s$lat[kA], time = A$s$time[kA],
      mag_traverse = interp(A, fa, P, "mag"), mag_tie = interp(B, fb, P, "mag"),
      alt_traverse = interp(A, fa, P, "radar_alt"), alt_tie = interp(B, fb, P, "radar_alt"),
      i_traverse = A$s$.i[kA], i_tie = B$s$.i[kB])
  }
  if (!length(out)) return(empty)
  xo <- dplyr::bind_rows(out)
  xo$misfit <- xo$mag_traverse - xo$mag_tie
  xo[, names(empty)]
}
