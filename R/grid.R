#' Grid the field
#'
#' Interpolates a field channel onto a regular grid in the survey's local
#' coordinates. The default method is minimum curvature (Briggs 1974): the
#' samples are binned to the nearest node, and the free nodes are relaxed to
#' the biharmonic solution - the smoothest surface through the data - with a
#' coarse-to-fine start so the relaxation converges in a few hundred sweeps.
#' `"idw"` is an inverse-distance alternative. Nodes farther than
#' `blank_distance` from any sample are blanked, so the map shows the field
#' only where it was flown.
#'
#' Up to 40,000 nodes the minimum-curvature system is solved exactly (sparse
#' LU); larger grids fall back to Jacobi relaxation from a coarse-to-fine
#' start, bounded by `max_iter` and `tol`.
#'
#' @param x A `magqc_result` (cell and blanking default from its spec) or a
#'   `magqc_survey` (then `cell` is required).
#' @param channel Field column to grid. Defaults to `mag_lev` when levelling
#'   ran, otherwise `mag`.
#' @param cell Grid cell size, m. Defaults to a quarter of the traverse
#'   spacing, the usual choice for line data.
#' @param method `"mincurv"` (default) or `"idw"`.
#' @param blank_distance Nodes farther than this (Chebyshev distance, m)
#'   from any sample are `NA`. Defaults to half the traverse spacing, which
#'   keeps every node between adjacent lines and blanks beyond the block.
#' @param exclude Row indices of samples to leave out. For a `magqc_result`
#'   this defaults to the samples flagged as spikes - a known point defect
#'   is not field, and a minimum-curvature surface would otherwise ring
#'   around each one.
#' @param max_iter,tol Relaxation stops when the largest change in a sweep
#'   falls below `tol` nT or after `max_iter` sweeps at the finest level.
#' @return An object of class `magqc_grid`: `x`, `y` (node coordinates, m),
#'   `z` (matrix `length(x)` by `length(y)`, nT, `NA` where blanked),
#'   `cell`, `channel`, `method`, `blank_distance`, `n_data`, `n_excluded`, `n_nodes`,
#'   `n_blank`, `rms_fit` (RMS of sample minus grid interpolated at the
#'   sample), `iterations`, and `bounds` (west/east/south/north in degrees)
#'   when the survey has longitude and latitude.
#' @references Briggs, I. C. (1974). Machine contouring using minimum
#'   curvature. Geophysics, 39(1), 39-48.
#' @examples
#' res <- run_qc(sim_survey(seed = 4), grid = FALSE)
#' g <- grid_field(res)
#' g
#' @export
grid_field <- function(x, channel = NULL, cell = NULL, method = c("mincurv", "idw"),
                       blank_distance = NULL, exclude = NULL, max_iter = 2000L, tol = 0.01) {
  method <- match.arg(method)
  if (inherits(x, "magqc_result")) {
    survey <- x$survey
    cell <- cell %||% x$spec$line_spacing / 4
    blank_distance <- blank_distance %||% x$spec$line_spacing / 2
    exclude <- exclude %||% x$flags$i_start[x$flags$check == "spikes"]
  } else if (inherits(x, "magqc_survey")) {
    survey <- x
    if (is.null(cell)) stop("`cell` is required when gridding a bare survey (a quarter of the line spacing is usual).", call. = FALSE)
    blank_distance <- blank_distance %||% 2 * cell
  } else {
    stop("`x` must be a magqc_survey or a magqc_result.", call. = FALSE)
  }
  channel <- channel %||% if ("mag_lev" %in% names(survey)) "mag_lev" else "mag"
  if (!channel %in% names(survey)) stop("Channel `", channel, "` is not in the survey.", call. = FALSE)
  .grid_survey(survey, channel, cell, method, blank_distance, max_iter, tol, exclude)
}

.grid_survey <- function(survey, channel, cell, method, blank_distance, max_iter, tol, exclude = NULL) {
  v <- survey[[channel]]
  exclude <- exclude[!is.na(exclude)]
  if (length(exclude)) v[exclude] <- NA_real_   # flagged spikes are not field
  ok <- is.finite(v) & is.finite(survey$x) & is.finite(survey$y)
  px <- survey$x[ok]; py <- survey$y[ok]; pv <- v[ok]
  if (length(pv) < 10) stop("Too few finite samples to grid.", call. = FALSE)
  if (!is.numeric(cell) || cell <= 0) stop("`cell` must be positive.", call. = FALSE)

  # nodes from just inside a cell of the data extent: a wider margin only
  # shows the free-edge overshoot of the smoothest surface beyond the last line
  x0 <- floor(min(px) / cell) * cell
  y0 <- floor(min(py) / cell) * cell
  nx <- as.integer(ceiling((max(px) - x0) / cell)) + 1L
  ny <- as.integer(ceiling((max(py) - y0) / cell)) + 1L
  if (nx * ny > 4e6) stop("Grid would have ", nx, " x ", ny, " nodes; use a larger `cell`.", call. = FALSE)
  gx <- x0 + (seq_len(nx) - 1) * cell
  gy <- y0 + (seq_len(ny) - 1) * cell

  # bin samples to the nearest node; the node value is a local plane through
  # the cell's samples evaluated at the node, not their mean (see .bin_plane)
  ix <- as.integer(round((px - x0) / cell)) + 1L
  iy <- as.integer(round((py - y0) / cell)) + 1L
  con <- .bin_plane(px - gx[ix], py - gy[iy], pv, ix + (iy - 1L) * nx, nx * ny)
  dim(con) <- c(nx, ny)

  z <- switch(method,
              mincurv = .mincurv(con, max_iter, tol),
              idw = .idw_grid(con, gx, gy, blank_distance))
  iterations <- attr(z, "iterations") %||% NA_integer_
  attr(z, "iterations") <- NULL

  k <- as.integer(ceiling(blank_distance / cell))
  keep <- .dilate(!is.na(con), k)
  z[!keep] <- NA_real_

  fit <- .bilinear_at(z, (px - x0) / cell + 1, (py - y0) / cell + 1)
  resid <- pv - fit

  bounds <- NULL
  if (all(c("lon", "lat") %in% names(survey)) && any(is.finite(survey$lon) & is.finite(survey$lat))) {
    ll <- survey[is.finite(survey$lon) & is.finite(survey$lat), ]
    if (nrow(ll) >= 3 && stats::sd(ll$x) > 0 && stats::sd(ll$y) > 0) {
      # the local frame is a linear map of lon/lat over a survey block
      flon <- stats::lm(lon ~ x + y, data = ll); flat <- stats::lm(lat ~ x + y, data = ll)
      corners <- data.frame(x = c(gx[1] - cell / 2, gx[nx] + cell / 2), y = c(gy[1] - cell / 2, gy[ny] + cell / 2))
      lon_c <- stats::predict(flon, corners); lat_c <- stats::predict(flat, corners)
      bounds <- list(west = min(lon_c), east = max(lon_c), south = min(lat_c), north = max(lat_c))
    }
  }

  structure(list(
    x = gx, y = gy, z = z, cell = cell, channel = channel, method = method,
    blank_distance = blank_distance, n_data = length(pv), n_excluded = length(exclude), n_nodes = nx * ny,
    n_blank = sum(!keep), rms_fit = sqrt(mean(resid^2, na.rm = TRUE)),
    iterations = iterations, bounds = bounds),
    class = "magqc_grid")
}

#' @export
print.magqc_grid <- function(x, ...) {
  cat("<magqc grid>\n")
  cat(sprintf("  %s: %d x %d nodes at %g m (%s), %d of %d blanked\n", x$channel, length(x$x), length(x$y),
              x$cell, x$method, x$n_blank, x$n_nodes))
  cat(sprintf("  %d samples%s; RMS fit %.2f nT; range %.1f to %.1f nT\n", x$n_data,
              if (x$n_excluded > 0) sprintf(" (%d flagged spikes excluded)", x$n_excluded) else "",
              x$rms_fit, min(x$z, na.rm = TRUE), max(x$z, na.rm = TRUE)))
  invisible(x)
}

#' Node constraints from a local plane per cell
#'
#' Line data sample a cell densely along the line and hardly at all across
#' it. The mean of a cell's samples sits at their centroid, which at a line
#' end or a partly covered cell is up to half a cell from the node - for a
#' 0.02 nT/m gradient and 50 m cells that is a 0.5 nT bias. Instead fit
#' `v = a + b dx + c dy` to each cell's samples (dx, dy relative to the
#' node) with a light ridge on the slopes, so a direction the samples do
#' not span defaults to zero slope and `a` is the value at the node. Solved
#' per cell with Cramer's rule, vectorised over cells.
#' @noRd
.bin_plane <- function(dx, dy, v, key, n_nodes, ridge = 1) {
  S <- function(w) rowsum(w, key)
  n <- S(rep(1, length(v))); sx <- S(dx); sy <- S(dy)
  sxx <- S(dx^2) + ridge * n; syy <- S(dy^2) + ridge * n; sxy <- S(dx * dy)
  sv <- S(v); svx <- S(v * dx); svy <- S(v * dy)
  # [n sx sy; sx sxx sxy; sy sxy syy] [a b c] = [sv svx svy]
  det <- n * (sxx * syy - sxy^2) - sx * (sx * syy - sxy * sy) + sy * (sx * sxy - sxx * sy)
  a <- (sv * (sxx * syy - sxy^2) - sx * (svx * syy - sxy * svy) + sy * (svx * sxy - sxx * svy)) / det
  a[!is.finite(a)] <- (sv / n)[!is.finite(a)]
  out <- rep(NA_real_, n_nodes)
  out[as.integer(rownames(n))] <- a[, 1]
  out
}

# ---- minimum curvature ------------------------------------------------------------

#' Minimum curvature: the constrained biharmonic system
#'
#' `con` is the node-binned data (NA where free). Every free node satisfies
#' the 13-point discrete biharmonic equation
#' `20 u - 8 (adjacent) + 2 (diagonal) + (two away) = 0`; constrained nodes
#' are held at their data value; the two ghost layers outside the grid are
#' linear extrapolations of the edge (a free edge: zero second derivative).
#' Up to `direct_max` nodes the system is solved exactly with a sparse LU
#' (Matrix); beyond that Jacobi relaxation with a coarse-to-fine start takes
#' over, since the LU fill-in of a 2-D biharmonic operator grows too fast.
#' @noRd
.mincurv <- function(con, max_iter, tol, direct_max = 40000L) {
  if (length(con) <= direct_max) .mincurv_direct(con) else .mincurv_relax(con, max_iter, tol)
}

.mincurv_direct <- function(con) {
  nx <- nrow(con); ny <- ncol(con); n <- nx * ny
  # 1-D expansion of a padded coordinate into real nodes: interior -> itself;
  # ghosts -> linear extrapolation of the two edge nodes. The 2-D ghost is the
  # tensor product of the two 1-D expansions.
  expand1 <- function(n1) {
    idx <- matrix(0L, n1 + 4, 2); w <- matrix(0, n1 + 4, 2)
    c1 <- (-1):(n1 + 2)
    inside <- c1 >= 1 & c1 <= n1
    idx[inside, 1] <- c1[inside]; w[inside, 1] <- 1
    lo <- c1 < 1; hi <- c1 > n1
    idx[lo, ] <- rep(c(1L, 2L), each = sum(lo)); w[lo, ] <- cbind(1 - c1[lo] + 1, -(1 - c1[lo]))
    idx[hi, ] <- rep(c(n1, n1 - 1L), each = sum(hi)); w[hi, ] <- cbind(c1[hi] - n1 + 1, -(c1[hi] - n1))
    list(idx = idx, w = w)   # row k corresponds to coordinate c = k - 2
  }
  EI <- expand1(nx); EJ <- expand1(ny)

  free <- which(is.na(con)); fixed <- which(!is.na(con))
  fi <- ((free - 1L) %% nx) + 1L; fj <- ((free - 1L) %/% nx) + 1L
  offsets <- rbind(c(0, 0, 20), c(1, 0, -8), c(-1, 0, -8), c(0, 1, -8), c(0, -1, -8),
                   c(1, 1, 2), c(1, -1, 2), c(-1, 1, 2), c(-1, -1, 2),
                   c(2, 0, 1), c(-2, 0, 1), c(0, 2, 1), c(0, -2, 1))
  rows <- integer(0); cols <- integer(0); vals <- numeric(0)
  for (o in seq_len(nrow(offsets))) {
    I <- fi + offsets[o, 1] + 2L; J <- fj + offsets[o, 2] + 2L   # padded-table rows
    for (a in 1:2) for (b in 1:2) {
      wi <- EI$w[I, a]; wj <- EJ$w[J, b]
      w <- wi * wj * offsets[o, 3]
      keep <- w != 0
      if (!any(keep)) next
      col <- EI$idx[I[keep], a] + (EJ$idx[J[keep], b] - 1L) * nx
      rows <- c(rows, free[keep]); cols <- c(cols, col); vals <- c(vals, w[keep])
    }
  }
  rows <- c(rows, fixed); cols <- c(cols, fixed); vals <- c(vals, rep(1, length(fixed)))
  A <- Matrix::sparseMatrix(i = rows, j = cols, x = vals, dims = c(n, n))
  b <- numeric(n); b[fixed] <- con[fixed]
  u <- as.vector(Matrix::solve(A, b))
  z <- matrix(u, nx, ny)
  attr(z, "iterations") <- 0L
  z
}

#' Jacobi relaxation with a coarse-to-fine start (large grids)
#'
#' Each level bins the constraints to a grid `L` times coarser, starts from
#' the previous level's solution resampled bilinearly (the coarsest from the
#' data mean), and relaxes until the largest change is below `tol` or
#' `max_iter` sweeps.
#' @noRd
.mincurv_relax <- function(con, max_iter, tol) {
  nx <- nrow(con); ny <- ncol(con)
  levels <- 2^(max(0L, floor(log2(min(nx, ny) / 8))):0)
  u <- NULL; total <- 0L
  for (L in levels) {
    cnx <- as.integer(ceiling((nx - 1) / L)) + 1L
    cny <- as.integer(ceiling((ny - 1) / L)) + 1L
    idx <- which(!is.na(con), arr.ind = TRUE)
    ci <- as.integer(floor((idx[, 1] - 1) / L + 0.5)) + 1L
    cj <- as.integer(floor((idx[, 2] - 1) / L + 0.5)) + 1L
    ckey <- ci + (cj - 1L) * cnx
    cmeans <- rowsum(con[idx], ckey) / as.vector(table(ckey)[as.character(sort(unique(ckey)))])
    ccon <- matrix(NA_real_, cnx, cny)
    ccon[as.integer(rownames(cmeans))] <- cmeans[, 1]
    u <- if (is.null(u)) matrix(mean(ccon, na.rm = TRUE), cnx, cny) else .resample(u, cnx, cny)
    r <- .relax(u, ccon, max_iter, tol)
    u <- r$u; total <- total + r$iter
  }
  attr(u, "iterations") <- total
  u
}

#' Jacobi sweeps of the biharmonic stencil, constrained nodes held
#' @noRd
.relax <- function(u, con, max_iter, tol) {
  nx <- nrow(u); ny <- ncol(u)
  fixed <- !is.na(con)
  u[fixed] <- con[fixed]
  if (nx < 3 || ny < 3) return(list(u = u, iter = 0L))
  i <- 3:(nx + 2); j <- 3:(ny + 2)
  it <- 0L
  for (it in seq_len(max_iter)) {
    p <- .pad2(u)
    new <- (8 * (p[i + 1, j] + p[i - 1, j] + p[i, j + 1] + p[i, j - 1]) -
              2 * (p[i + 1, j + 1] + p[i + 1, j - 1] + p[i - 1, j + 1] + p[i - 1, j - 1]) -
              (p[i + 2, j] + p[i - 2, j] + p[i, j + 2] + p[i, j - 2])) / 20
    new[fixed] <- con[fixed]
    delta <- max(abs(new - u))
    u <- new
    if (is.finite(delta) && delta < tol) break
  }
  list(u = u, iter = it)
}

#' Pad a matrix by two on each side with linear extrapolation (free edges)
#' @noRd
.pad2 <- function(u) {
  nx <- nrow(u); ny <- ncol(u)
  p <- matrix(0, nx + 4, ny + 4)
  p[3:(nx + 2), 3:(ny + 2)] <- u
  p[2, 3:(ny + 2)] <- 2 * u[1, ] - u[2, ];        p[1, 3:(ny + 2)] <- 3 * u[1, ] - 2 * u[2, ]
  p[nx + 3, 3:(ny + 2)] <- 2 * u[nx, ] - u[nx - 1, ]; p[nx + 4, 3:(ny + 2)] <- 3 * u[nx, ] - 2 * u[nx - 1, ]
  p[, 2] <- 2 * p[, 3] - p[, 4];                  p[, 1] <- 3 * p[, 3] - 2 * p[, 4]
  p[, ny + 3] <- 2 * p[, ny + 2] - p[, ny + 1];   p[, ny + 4] <- 3 * p[, ny + 2] - 2 * p[, ny + 1]
  p
}

#' Bilinear resample of a matrix to new dimensions (index space)
#' @noRd
.resample <- function(u, nx2, ny2) {
  fi <- if (nx2 > 1) (seq_len(nx2) - 1) * (nrow(u) - 1) / (nx2 - 1) + 1 else 1
  fj <- if (ny2 > 1) (seq_len(ny2) - 1) * (ncol(u) - 1) / (ny2 - 1) + 1 else 1
  matrix(.bilinear_at(u, rep(fi, times = ny2), rep(fj, each = nx2)), nx2, ny2)
}

#' Bilinear interpolation of a matrix at fractional index positions
#'
#' `NA` where any of the four surrounding nodes is `NA` or the position is
#' outside the matrix.
#' @noRd
.bilinear_at <- function(z, fi, fj) {
  nx <- nrow(z); ny <- ncol(z)
  out <- rep(NA_real_, length(fi))
  inside <- is.finite(fi) & is.finite(fj) & fi >= 1 & fi <= nx & fj >= 1 & fj <= ny
  fi <- fi[inside]; fj <- fj[inside]
  i0 <- pmin(floor(fi), nx - 1L); j0 <- pmin(floor(fj), ny - 1L)
  if (nx == 1L) i0 <- rep(1L, length(fi)); if (ny == 1L) j0 <- rep(1L, length(fj))
  t <- fi - i0; s <- fj - j0
  i1 <- pmin(i0 + 1L, nx); j1 <- pmin(j0 + 1L, ny)
  v <- (1 - t) * (1 - s) * z[cbind(i0, j0)] + t * (1 - s) * z[cbind(i1, j0)] +
    (1 - t) * s * z[cbind(i0, j1)] + t * s * z[cbind(i1, j1)]
  out[inside] <- v
  out
}

#' Grow a logical mask by k cells (8-connected)
#' @noRd
.dilate <- function(m, k) {
  nx <- nrow(m); ny <- ncol(m)
  for (step in seq_len(k)) {
    p <- matrix(FALSE, nx + 2, ny + 2); p[2:(nx + 1), 2:(ny + 1)] <- m
    m <- m | p[1:nx, 2:(ny + 1)] | p[3:(nx + 2), 2:(ny + 1)] | p[2:(nx + 1), 1:ny] | p[2:(nx + 1), 3:(ny + 2)] |
      p[1:nx, 1:ny] | p[3:(nx + 2), 1:ny] | p[1:nx, 3:(ny + 2)] | p[3:(nx + 2), 3:(ny + 2)]
  }
  m
}

#' Inverse-distance-squared gridding from the node-binned data
#' @noRd
.idw_grid <- function(con, gx, gy, radius) {
  idx <- which(!is.na(con), arr.ind = TRUE)
  dx <- gx[idx[, 1]]; dy <- gy[idx[, 2]]; dv <- con[idx]
  nx <- length(gx); ny <- length(gy)
  z <- matrix(NA_real_, nx, ny)
  gxx <- rep(gx, times = ny); gyy <- rep(gy, each = nx)
  chunk <- max(1L, as.integer(2e6 / length(dv)))
  for (s in seq(1, length(gxx), by = chunk)) {
    e <- min(s + chunk - 1, length(gxx))
    d2 <- outer(gxx[s:e], dx, "-")^2 + outer(gyy[s:e], dy, "-")^2
    w <- 1 / pmax(d2, (0.1)^2)
    w[d2 > radius^2] <- 0
    num <- w %*% dv; den <- rowSums(w)
    z[s:e] <- ifelse(den > 0, num / den, NA_real_)
  }
  z
}

# ---- colouring -------------------------------------------------------------------

#' Colour lookup for a grid: Spectral, histogram-equalised by default
#' @noRd
.grid_scale <- function(z, stretch = c("equalize", "linear"), n = 256L) {
  stretch <- match.arg(stretch)
  v <- z[is.finite(z)]
  pal <- grDevices::hcl.colors(n, "Spectral", rev = TRUE)
  if (!length(v)) return(list(pal = pal, breaks = numeric(0), index = function(x) rep(NA_integer_, length(x))))
  breaks <- if (stretch == "equalize") stats::quantile(v, seq(0, 1, length.out = n + 1), names = FALSE) else
    seq(min(v), max(v), length.out = n + 1)
  breaks <- unique(breaks)
  if (length(breaks) < 2) breaks <- c(breaks[1] - 0.5, breaks[1] + 0.5)
  index <- function(x) {
    k <- findInterval(x, breaks, rightmost.closed = TRUE, all.inside = TRUE)
    k[!is.finite(x)] <- NA_integer_
    as.integer(round((k - 1) / (length(breaks) - 2 + 1e-9) * (n - 1))) + 1L
  }
  list(pal = pal, breaks = breaks, index = index)
}

#' A grid as a base64 PNG (north up) plus legend swatches
#' @noRd
.grid_png <- function(g, stretch = "equalize") {
  sc <- .grid_scale(g$z, stretch)
  idx <- sc$index(g$z)                             # nx x ny
  rgb <- grDevices::col2rgb(sc$pal) / 255          # 3 x n
  nx <- length(g$x); ny <- length(g$y)
  arr <- array(0, c(ny, nx, 4))                    # rows = y from the top
  ok <- !is.na(idx)
  for (ch in 1:3) {
    m <- matrix(0, nx, ny); m[ok] <- rgb[ch, idx[ok]]
    arr[, , ch] <- t(m)[ny:1, , drop = FALSE]
  }
  a <- matrix(0, nx, ny); a[ok] <- 1
  arr[, , 4] <- t(a)[ny:1, , drop = FALSE]
  raw <- png::writePNG(arr, target = raw())
  v <- g$z[is.finite(g$z)]
  q <- stats::quantile(v, c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 1), names = FALSE)
  list(base64 = base64enc::base64encode(raw),
       legend_colors = sc$pal[sc$index(q)],
       legend_labels = formatC(q, format = "f", digits = 0, big.mark = ","))
}

#' Plotly colourscale for a grid (list of position/colour pairs)
#' @noRd
.grid_colorscale <- function(n = 11L) {
  pal <- grDevices::hcl.colors(n, "Spectral", rev = TRUE)
  lapply(seq_len(n), function(k) list((k - 1) / (n - 1), pal[k]))
}

#' Gridded field
#'
#' Plotly heatmap of a [grid_field()] result in the survey's local frame,
#' equal-aspect, histogram-equalised colours (the colourbar ticks are the
#' field values at those quantiles).
#'
#' @param x A `magqc_result` with a `grid` element, or a `magqc_grid`.
#' @return A plotly widget, or `NULL` when there is no grid.
#' @export
plot_grid <- function(x) {
  g <- if (inherits(x, "magqc_grid")) x else x$grid
  if (is.null(g)) return(NULL)
  sc <- .grid_scale(g$z, "equalize")
  # equalised values in [0, 1] for colour, real values for hover
  zi <- matrix((sc$index(g$z) - 1) / (length(sc$pal) - 1), nrow(g$z), ncol(g$z))
  q <- stats::quantile(g$z[is.finite(g$z)], seq(0, 1, by = 0.2), names = FALSE)
  p <- plotly::plot_ly(x = g$x, y = g$y, z = t(zi), customdata = t(round(g$z, 1)), type = "heatmap",
                       colorscale = .grid_colorscale(), zmin = 0, zmax = 1, hoverongaps = FALSE,
                       hovertemplate = I("%{customdata:.1f} nT<br>x %{x:.0f}, y %{y:.0f}<extra></extra>"),
                       colorbar = list(title = list(text = "nT"), tickvals = seq(0, 1, by = 0.2),
                                       ticktext = formatC(q, format = "f", digits = 0, big.mark = ","),
                                       thickness = 12, len = 0.9))
  p <- .plotly_base(p, "northing (m)", "easting (m)", date_axis = FALSE)
  # equal aspect by shrinking the axis domain, not by widening the range:
  # a north-south block then sits centred at true shape instead of being
  # padded out to +/- several km of nothing
  plotly::layout(p, xaxis = list(constrain = "domain"),
                 yaxis = list(scaleanchor = "x", scaleratio = 1, constrain = "domain"),
                 hovermode = "closest",
                 annotations = list(list(
                   xref = "paper", yref = "paper", x = 0, y = 1.02, xanchor = "left", yanchor = "bottom",
                   showarrow = FALSE, font = list(color = .pal$ink2, size = 12),
                   text = sprintf("%s, %g m cells, %s; RMS fit %.2f nT", g$channel, g$cell,
                                  if (g$method == "mincurv") "minimum curvature" else "inverse distance", g$rms_fit))))
}
