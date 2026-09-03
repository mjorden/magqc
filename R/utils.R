#' @keywords internal
"_PACKAGE"

`%||%` <- function(a, b) if (is.null(a)) b else a

.deg2rad <- function(d) d * pi / 180
.rad2deg <- function(r) r * 180 / pi

#' Convert local east/north metres to WGS84 longitude and latitude
#'
#' A flat-earth conversion using the standard series expansions for metres per
#' degree of latitude and longitude on the WGS84 ellipsoid. Accurate to well
#' under a metre over a survey block of a few tens of kilometres, which is the
#' only scale this package uses it at. For anything larger, project the data
#' properly before loading it.
#'
#' @param x,y East and north offsets from the origin, in metres.
#' @param lon0,lat0 Origin longitude and latitude, in decimal degrees.
#' @return A list with elements `lon` and `lat`.
#' @export
enu_to_lonlat <- function(x, y, lon0, lat0) {
  m_lat <- .m_per_deg_lat(lat0)
  lat <- lat0 + y / m_lat
  m_lon <- .m_per_deg_lon((lat0 + lat) / 2)
  list(lon = lon0 + x / m_lon, lat = lat)
}

#' @rdname enu_to_lonlat
#' @param lon,lat Longitude and latitude in decimal degrees.
#' @export
lonlat_to_enu <- function(lon, lat, lon0, lat0) {
  m_lat <- .m_per_deg_lat(lat0)
  m_lon <- .m_per_deg_lon((lat0 + lat) / 2)
  list(x = (lon - lon0) * m_lon, y = (lat - lat0) * m_lat)
}

.m_per_deg_lat <- function(lat_deg) {
  p <- .deg2rad(lat_deg)
  111132.92 - 559.82 * cos(2 * p) + 1.175 * cos(4 * p) - 0.0023 * cos(6 * p)
}

.m_per_deg_lon <- function(lat_deg) {
  p <- .deg2rad(lat_deg)
  111412.84 * cos(p) - 93.5 * cos(3 * p) + 0.118 * cos(5 * p)
}

#' Great-circle-free planar distance between consecutive points
#' @noRd
.step_distance <- function(x, y) {
  c(NA_real_, sqrt(diff(x)^2 + diff(y)^2))
}

#' Bearing in degrees clockwise from north
#' @noRd
.bearing <- function(dx, dy) (.rad2deg(atan2(dx, dy))) %% 360

#' Running median with automatic window clamping
#'
#' `stats::runmed()` requires an odd window no longer than the series. This
#' wrapper clamps the window and returns the input unchanged when the series is
#' too short to filter, so callers never have to special-case short lines.
#' @noRd
.runmed_safe <- function(x, k) {
  n <- length(x)
  if (n < 3) return(x)
  k <- min(k, if (n %% 2 == 0) n - 1 else n)
  if (k %% 2 == 0) k <- k - 1
  if (k < 3) return(x)
  # "constant" carries the first/last full-window median out to the ends.
  # The default "median" shrinks the window at the edges, so a spike on the
  # last sample becomes its own scale estimate and can never be detected.
  stats::runmed(x, k, endrule = "constant")
}

#' Running maximum over a centred window
#' @noRd
.run_max <- function(x, k) {
  n <- length(x); h <- (k - 1) %/% 2
  vapply(seq_len(n), function(i) {
    v <- x[max(1, i - h):min(n, i + h)]
    if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
  }, numeric(1))
}

#' Spike detection on the fourth difference
#'
#' A running-median (Hampel) test on the raw field fails on airborne data: on
#' a monotone stretch of a smooth anomaly the running median *is* the centre
#' sample, the residual collapses to zero and ordinary noise reads as many
#' sigma. The fourth-difference operator annihilates anything up to a cubic,
#' so applying the robust test to it isolates the sensor noise from the
#' geology. A single-sample spike of amplitude `A` leaves a `+A, -4A, +6A,
#' -4A, +A` pattern in the fourth difference; each cluster of samples that
#' trips the threshold is reduced to the sample with the largest fourth
#' difference, and `A` is estimated as one sixth of it.
#'
#' The scale is a running MAD of the fourth difference, dilated by a running
#' maximum so the quiet side of a noisy stretch is judged against the noisy
#' sigma too - otherwise the first samples of a noisy stretch read as spikes
#' against the quiet scale next door. That dilation is also why the practical
#' detection floor is higher than the naive `6A / (sqrt(70) sigma)` z-score
#' suggests: with the defaults (`k = 21`, `nsigma = 6`) a spike of 25x the
#' noise sd is always found, ~12x about 40% of the time, and 8x or less
#' never (measured in tests/testthat/test-utils.R).
#'
#' @return A list: `outlier` (logical), `amplitude` and `z` (NA except at
#'   outliers), `cleaned` (the series with each spike estimate subtracted).
#' @noRd
.spike_detect <- function(x, k = 21, nsigma = 6) {
  n <- length(x)
  none <- list(outlier = logical(n), amplitude = rep(NA_real_, n),
               z = rep(NA_real_, n), cleaned = x)
  ok <- is.finite(x)
  if (sum(ok) < 7) return(none)
  xf <- x
  if (!all(ok)) xf <- stats::approx(which(ok), x[ok], xout = seq_len(n), rule = 2)$y

  # The centred operator is undefined on the first and last two samples, so
  # extend the series by linear extrapolation and difference the extension.
  # At the edges this degrades gracefully to a lower-order test (a spike on
  # the last sample shows up as its second difference) rather than no test.
  ext <- c(3 * xf[1] - 2 * xf[2], 2 * xf[1] - xf[2], xf,
           2 * xf[n] - xf[n - 1], 3 * xf[n] - 2 * xf[n - 1])
  d4 <- fourth_difference(ext, normalize = FALSE)[3:(n + 2)]
  a <- abs(d4)
  scale <- .run_max(1.4826 * .runmed_safe(a, k), k)
  floor_s <- 0.5 * 1.4826 * stats::median(a)
  if (!is.finite(floor_s) || floor_s <= 0) floor_s <- .Machine$double.eps
  scale <- pmax(scale, floor_s)
  z4 <- d4 / scale
  cand <- ok & is.finite(z4) & abs(z4) > nsigma

  out <- none
  runs <- .runs(cand)
  for (r in seq_len(nrow(runs))) {
    i <- runs$start[r]:runs$end[r]
    j <- i[which.max(abs(d4[i]))]
    out$outlier[j] <- TRUE
    out$amplitude[j] <- if (j > 2 && j < n - 1) d4[j] / 6 else
      if (j <= 2) x[j] - (2 * xf[j + 1] - xf[j + 2]) else x[j] - (2 * xf[j - 1] - xf[j - 2])
    out$z[j] <- z4[j]
    out$cleaned[j] <- x[j] - out$amplitude[j]
  }
  out
}

#' Fourth difference of a series
#'
#' The classic airborne magnetic noise-envelope statistic:
#' `(x[i-2] - 4 x[i-1] + 6 x[i] - 4 x[i+1] + x[i+2])`. Specifications differ on
#' whether the operator is normalised by its coefficient sum of 16; when
#' `normalize = TRUE` (the default here) the result is on the same scale as the
#' underlying noise amplitude, which is the convention the default
#' `fourth_diff_tol` in [survey_spec()] assumes.
#'
#' @param x Numeric vector of field readings, in survey order.
#' @param normalize Divide by 16. See details.
#' @return A numeric vector the same length as `x`, `NA` in the first and last
#'   two positions.
#' @export
fourth_difference <- function(x, normalize = TRUE) {
  n <- length(x)
  out <- rep(NA_real_, n)
  if (n < 5) return(out)
  i <- 3:(n - 2)
  out[i] <- x[i - 2] - 4 * x[i - 1] + 6 * x[i] - 4 * x[i + 1] + x[i + 2]
  if (normalize) out <- out / 16
  out
}

#' Total-least-squares fit of a straight line through a point cloud
#'
#' Returns the centroid, the unit direction vector of the principal axis, and
#' the signed perpendicular offset of every point from that axis. Used to
#' measure how straight a flown line is when no flight plan is supplied.
#' @noRd
.tls_line <- function(x, y) {
  cx <- mean(x); cy <- mean(y)
  m <- cbind(x - cx, y - cy)
  s <- svd(m, nu = 0, nv = 2)
  dir <- s$v[, 1]
  nrm <- c(-dir[2], dir[1])
  list(cx = cx, cy = cy, dir = dir, normal = nrm,
       offset = as.vector(m %*% nrm),
       along  = as.vector(m %*% dir))
}

#' Intersection of two infinite lines given as point + direction
#' @noRd
.line_intersection <- function(p1, d1, p2, d2) {
  den <- d1[1] * d2[2] - d1[2] * d2[1]
  if (abs(den) < 1e-12) return(c(NA_real_, NA_real_))
  t <- ((p2[1] - p1[1]) * d2[2] - (p2[2] - p1[2]) * d2[1]) / den
  c(p1[1] + t * d1[1], p1[2] + t * d1[2])
}

#' Longest run of TRUE values, expressed as index ranges
#'
#' Collapses a logical vector into the contiguous `TRUE` intervals it contains,
#' which is how every check turns per-sample tests into flagged intervals.
#' @noRd
.runs <- function(flag) {
  flag[is.na(flag)] <- FALSE
  if (!any(flag)) return(data.frame(start = integer(0), end = integer(0)))
  r <- rle(flag)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  keep <- r$values
  data.frame(start = starts[keep], end = ends[keep])
}
