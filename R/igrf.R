#' IGRF main field at survey points
#'
#' Evaluates the International Geomagnetic Reference Field (via the `igrf`
#' package) on a small longitude/latitude lattice over the points, at one
#' altitude and epoch, and interpolates total intensity, inclination and
#' declination bilinearly to each point. Over a survey block the main field
#' varies by a few nT per kilometre and is smooth, so the lattice is exact
#' to well under a tenth of a nanotesla and thousands of samples cost the
#' same as twenty-five model evaluations.
#'
#' @param lon,lat Longitude and latitude, decimal degrees.
#' @param time A POSIXct (or a vector; the mean sets the epoch).
#' @param altitude Height above the WGS84 spheroid, metres. The main field
#'   falls by about 0.025 nT per metre, so a rough value is enough.
#' @param n Lattice size (`n` by `n` evaluations).
#' @return A data frame with one row per point: `F` (total intensity, nT),
#'   `I` (inclination, degrees, positive down), `D` (declination, degrees,
#'   positive east); attributes `epoch` (decimal year) and `model`
#'   (the `igrf` package version).
#' @examples
#' igrf_field(-108.5, 43.2, as.POSIXct("2024-06-14", tz = "UTC"), altitude = 1600)
#' @export
igrf_field <- function(lon, lat, time, altitude = 0, n = 5L) {
  rlang::check_installed("igrf", reason = "to evaluate the IGRF main field")
  ok <- is.finite(lon) & is.finite(lat)
  if (!any(ok)) stop("No finite longitude/latitude to evaluate the IGRF at.", call. = FALSE)
  t0 <- mean(as.numeric(time), na.rm = TRUE)
  lt <- as.POSIXlt(t0, origin = "1970-01-01", tz = "UTC")
  epoch <- lt$year + 1900 + (lt$yday + (lt$hour * 3600 + lt$min * 60 + lt$sec) / 86400) / 365.25
  alt_km <- altitude / 1000

  pad <- function(r) if (diff(r) < 1e-4) r + c(-1, 1) * 5e-3 else r
  lon_r <- pad(range(lon[ok])); lat_r <- pad(range(lat[ok]))
  gl <- seq(lon_r[1], lon_r[2], length.out = n); gb <- seq(lat_r[1], lat_r[2], length.out = n)
  F_ <- I_ <- D_ <- matrix(NA_real_, n, n)
  for (i in seq_len(n)) for (j in seq_len(n)) {
    r <- igrf::igrf(field = "main", year = epoch, type = "spheroid", altitude = alt_km,
                    latitude = gb[j], longitude = gl[i])
    F_[i, j] <- r$F; I_[i, j] <- r$I; D_[i, j] <- r$D
  }
  # clamp to the lattice: a point on its edge can land a rounding error outside
  fi <- pmin(pmax((lon - gl[1]) / (gl[2] - gl[1]) + 1, 1), n)
  fj <- pmin(pmax((lat - gb[1]) / (gb[2] - gb[1]) + 1, 1), n)
  out <- data.frame(F = .bilinear_at(F_, fi, fj), I = .bilinear_at(I_, fi, fj), D = .bilinear_at(D_, fi, fj))
  attr(out, "epoch") <- epoch
  attr(out, "model") <- paste0("igrf ", as.character(utils::packageVersion("igrf")))
  out
}

#' Remove the IGRF from a survey
#'
#' Adds `mag_igrf` (the main field at each sample) and `mag_res` (the
#' residual: the levelled field when present, otherwise `mag`, minus the
#' main field) to the survey, and returns the summary the report prints.
#' @return A list: the `survey` with the two columns and `igrf`, a list of
#'   `model`, `epoch`, `altitude` (m), and `F`, `I`, `D` at the block centre.
#' @noRd
.remove_igrf <- function(survey, altitude = NULL) {
  ok <- is.finite(survey$lon) & is.finite(survey$lat)
  if (!any(ok)) return(list(survey = survey, igrf = NULL))
  if (is.null(altitude)) {
    altitude <- if ("gps_alt" %in% names(survey) && any(is.finite(survey$gps_alt))) {
      mean(survey$gps_alt, na.rm = TRUE)
    } else if (any(is.finite(survey$radar_alt))) {
      mean(survey$radar_alt, na.rm = TRUE)   # terrain clearance stands in for height above the spheroid
    } else 0
  }
  f <- igrf_field(survey$lon, survey$lat, survey$time, altitude = altitude)
  base <- if ("mag_lev" %in% names(survey)) survey$mag_lev else survey$mag
  survey$mag_igrf <- f$F
  survey$mag_res <- base - f$F
  ctr <- which.min((survey$lon - mean(survey$lon[ok]))^2 + (survey$lat - mean(survey$lat[ok]))^2)
  list(survey = survey,
       igrf = list(model = attr(f, "model"), epoch = attr(f, "epoch"), altitude = altitude,
                   F = f$F[ctr], I = f$I[ctr], D = f$D[ctr],
                   from = if ("mag_lev" %in% names(survey)) "mag_lev" else "mag"))
}

#' Human label for a field channel
#' @noRd
.channel_label <- function(channel) {
  switch(channel,
         rtp = "residual field reduced to the pole",
         mag_res = "residual field (IGRF removed)",
         mag_lev = "levelled field",
         mag = "field",
         mag_raw = "raw field",
         channel)
}
