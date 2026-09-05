#' Reduction to the pole
#'
#' Transforms a total-field anomaly grid observed in an inclined field into
#' the anomaly the same sources would produce at the magnetic pole, where
#' anomalies sit directly over their sources and are symmetric. The filter
#' is applied in the wavenumber domain (Blakely 1995): with
#' `theta = sin(I) + i cos(I) cos(D - phi)` for a wavevector at azimuth
#' `phi`, the induced-magnetisation operator is `1 / theta^2`. Near the
#' magnetic equator that operator blows up for wavevectors along the
#' declination, so the amplitude uses an *amplitude inclination* `I_a`
#' (MacLeod, Jones & Dai 1993): the phase is taken from the true inclination
#' and the amplitude from `max(|I|, I_a)`, which is exact when `|I| >= I_a`.
#'
#' Blanked nodes are infilled by minimum curvature before filtering and
#' blanked again after; the grid is mirror-extended and cosine-tapered by
#' a quarter of its size on each side so the periodic transform sees no
#' edge step. The grid mean is removed and restored.
#'
#' @param grid A `magqc_grid` (a residual field; see [grid_field()]).
#' @param inclination,declination Ambient field direction, degrees
#'   (inclination positive down, declination positive east).
#' @param amp_inclination Amplitude inclination, degrees; the amplitude of
#'   the operator is evaluated at `max(abs(inclination), amp_inclination)`.
#' @param pad Extension on each side as a fraction of the grid size.
#' @return A `magqc_grid` with `channel = "rtp"` and attributes
#'   `inclination`, `declination`, `amp_inclination`.
#' @references Blakely, R. J. (1995). Potential Theory in Gravity and
#'   Magnetic Applications, 12.3. MacLeod, I. N., Jones, K. and Dai, T. F.
#'   (1993). 3-D analytic signal in the interpretation of total magnetic
#'   field data at low magnetic latitudes. Exploration Geophysics 24,
#'   679-688.
#' @examples
#' res <- run_qc(sim_survey(seed = 4))
#' r <- reduce_to_pole(res$grid, res$igrf$I, res$igrf$D)
#' r
#' @export
reduce_to_pole <- function(grid, inclination, declination, amp_inclination = 20, pad = 0.25) {
  if (!inherits(grid, "magqc_grid")) stop("`grid` must be a magqc_grid.", call. = FALSE)
  for (v in c(inclination, declination, amp_inclination)) {
    if (!is.numeric(v) || length(v) != 1 || !is.finite(v)) stop("Field angles must be single finite numbers.", call. = FALSE)
  }
  z <- grid$z
  blank <- is.na(z)
  if (all(blank)) stop("The grid is entirely blank.", call. = FALSE)
  if (any(blank)) {
    con <- z
    z <- .mincurv(con, max_iter = 2000L, tol = 0.01)   # smooth infill of the blanks
    attr(z, "iterations") <- NULL
  }
  mu <- mean(z)
  z <- z - mu

  nx <- nrow(z); ny <- ncol(z)
  px <- max(8L, as.integer(ceiling(pad * nx))); py <- max(8L, as.integer(ceiling(pad * ny)))
  zp <- .extend_taper(z, px, py)
  nxp <- nrow(zp); nyp <- ncol(zp)

  kx <- 2 * pi * .fftfreq(nxp, grid$cell); ky <- 2 * pi * .fftfreq(nyp, grid$cell)
  KX <- matrix(kx, nxp, nyp); KY <- matrix(ky, nxp, nyp, byrow = TRUE)
  phi <- atan2(KX, KY)                              # wavevector azimuth from north, east positive
  I <- .deg2rad(inclination); D <- .deg2rad(declination)
  Ia <- .deg2rad(max(abs(inclination), amp_inclination))
  theta <- sin(I) + 1i * cos(I) * cos(D - phi)
  theta_a <- sin(Ia) + 1i * cos(Ia) * cos(D - phi)
  op <- Conj(theta)^2 / (Mod(theta)^2 * Mod(theta_a)^2)
  op[1, 1] <- 1                                     # DC: nothing to rotate

  Z <- stats::fft(zp)
  out <- Re(stats::fft(Z * op, inverse = TRUE)) / length(zp)
  out <- out[px + seq_len(nx), py + seq_len(ny)] + mu
  out[blank] <- NA_real_

  g <- grid
  g$z <- out
  g$channel <- "rtp"
  g$source_channel <- grid$channel
  g$rms_fit <- NA_real_
  attr(g, "inclination") <- inclination
  attr(g, "declination") <- declination
  attr(g, "amp_inclination") <- max(abs(inclination), amp_inclination)
  g
}

#' Frequencies of an n-point transform with sample spacing d (cycles per unit)
#' @noRd
.fftfreq <- function(n, d) {
  c(seq(0, floor((n - 1) / 2)), seq(-floor(n / 2), -1)) / (n * d)
}

#' Mirror-extend a matrix by (px, py) on each side and cosine-taper the extension
#' @noRd
.extend_taper <- function(z, px, py) {
  nx <- nrow(z); ny <- ncol(z)
  refl <- function(n, p) {
    # indices of a mirror extension, clamped for grids shorter than the pad
    lo <- pmin(pmax(rev(seq_len(p)) + 1L, 1L), n)
    hi <- pmin(pmax(n - seq_len(p), 1L), n)
    c(lo, seq_len(n), hi)
  }
  zp <- z[refl(nx, px), refl(ny, py), drop = FALSE]
  wx <- c(0.5 * (1 - cos(pi * seq_len(px) / (px + 1))), rep(1, nx), rev(0.5 * (1 - cos(pi * seq_len(px) / (px + 1)))))
  wy <- c(0.5 * (1 - cos(pi * seq_len(py) / (py + 1))), rep(1, ny), rev(0.5 * (1 - cos(pi * seq_len(py) / (py + 1)))))
  zp * outer(wx, wy)
}

#' Total-field anomaly of a point dipole with induced magnetisation
#'
#' The anomaly at points `(x, y)` on the plane `z = 0` of a dipole at
#' `(x0, y0)` and depth `h`, magnetised along the field direction
#' `(inclination, declination)`; `moment` scales it (a vertical field gives
#' a peak of `2 moment / h^3`). Used by the simulator and the RTP tests.
#' @noRd
.dipole_anomaly <- function(x, y, x0, y0, h, inclination, declination, moment = 1e9) {
  I <- .deg2rad(inclination); D <- .deg2rad(declination)
  f <- c(cos(I) * sin(D), cos(I) * cos(D), sin(I))     # east, north, down
  rx <- x - x0; ry <- y - y0; rz <- -h                  # source to observation (down positive)
  r <- sqrt(rx^2 + ry^2 + rz^2)
  fr <- f[1] * rx + f[2] * ry + f[3] * rz
  moment * (3 * fr^2 / r^5 - 1 / r^3)
}
