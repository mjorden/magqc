# Report palette. Categorical slots are assigned by entity and never cycled:
# traverse lines are always slot 1, tie lines always slot 2. Status colours
# are reserved for pass/fail and never reused for a series; wherever one
# appears it is paired with a text label so colour never carries the meaning
# alone.
.pal <- list(
  traverse  = "#2a78d6", tie = "#eb6834",
  # categorical slots in fixed order - assigned by position, never recycled (#7)
  series    = c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"),
  good      = "#0ca30c", critical = "#d03b3b",
  ink       = "#0b0b0b", ink2 = "#52514e", muted = "#898781",
  grid      = "#e1e0d9", axis = "#c3c2b7", surface = "#fcfcfb",
  band      = "rgba(137,135,129,0.12)", flag_fill = "rgba(208,59,59,0.16)")

.font <- 'system-ui, -apple-system, "Segoe UI", sans-serif'

.plotly_base <- function(p, ytitle, xtitle = "time (UTC)", date_axis = TRUE) {
  xaxis <- list(title = xtitle, gridcolor = .pal$grid, zeroline = FALSE,
                linecolor = .pal$axis, tickfont = list(color = .pal$muted))
  if (date_axis) xaxis$type <- "date"
  plotly::layout(
    p,
    font = list(family = .font, color = .pal$ink2, size = 12),
    paper_bgcolor = .pal$surface, plot_bgcolor = .pal$surface,
    margin = list(l = 60, r = 20, t = 10, b = 45),
    hovermode = "x unified",
    xaxis = xaxis,
    yaxis = list(title = ytitle, gridcolor = .pal$grid, zeroline = FALSE,
                 linecolor = .pal$axis, tickfont = list(color = .pal$muted)),
    legend = list(orientation = "h", x = 0, y = 1.08, font = list(color = .pal$ink2)))
}

#' Colours for a set of series names: fixed slot order, one per name
#'
#' More than eight distinct names (impossible for compass headings, the only
#' current caller) would exceed the palette; those beyond it get the muted
#' ink rather than a recycled hue that would alias two series (#7).
#' @noRd
.series_colours <- function(names) {
  n <- length(names)
  cols <- c(.pal$series[seq_len(min(n, length(.pal$series)))],
            rep(.pal$muted, max(0, n - length(.pal$series))))
  stats::setNames(cols, names)
}

#' Insert a row of NAs between lines so plotly does not join them across turns
#' @noRd
.gap_between_lines <- function(df) {
  parts <- split(df, factor(df$line, levels = unique(df$line)))
  out <- lapply(parts, function(p) {
    g <- p[1, ]; g[] <- NA; g$line <- p$line[1]
    rbind(p, g)
  })
  do.call(rbind, out)
}

.flag_mask <- function(result, check) {
  f <- result$flags[result$flags$check == check & !is.na(result$flags$i_start), ]
  m <- logical(nrow(result$survey))
  for (r in seq_len(nrow(f))) m[f$i_start[r]:f$i_end[r]] <- TRUE
  m
}

# ---- map -------------------------------------------------------------------

#' Map of flight lines and flagged intervals
#'
#' A leaflet map when the survey has longitude/latitude, otherwise a plotly
#' scatter in local metres. Traverse and tie lines are drawn in their own
#' colours; each check with findings gets a toggleable layer of red intervals
#' (or markers for single-sample and survey-level findings) with a popup
#' describing the finding.
#'
#' @param result A `magqc_result` from [run_qc()].
#' @param decimate Draw every n-th sample of each line; the full resolution is
#'   kept for flagged intervals.
#' @return An htmlwidget.
#' @export
qc_map <- function(result, decimate = 5) {
  s <- result$survey
  flags <- result$flags
  reg <- .check_registry()
  if (isTRUE(result$stats$has_lonlat)) .leaflet_map(s, flags, reg, decimate) else
    .plotly_map(s, flags, reg, decimate)
}

.popup <- function(f) {
  sprintf("<b>%s</b><br/>%s<br/><span style='color:%s'>%s%s%s</span>",
          htmltools::htmlEscape(ifelse(is.na(f$line), "survey-wide", f$line)),
          htmltools::htmlEscape(f$description), .pal$ink2,
          ifelse(is.na(f$time_start), "", format(f$time_start, "%H:%M:%S")),
          ifelse(is.na(f$time_end) | f$time_end == f$time_start, "",
                 paste0(" to ", format(f$time_end, "%H:%M:%S"))),
          ifelse(is.na(f$fid_start), "", sprintf("&nbsp;&nbsp;fid %s%s", f$fid_start,
                 ifelse(f$fid_end == f$fid_start, "", paste0("-", f$fid_end)))))
}

.leaflet_map <- function(s, flags, reg, decimate) {
  m <- leaflet::leaflet(options = leaflet::leafletOptions(preferCanvas = TRUE))
  # OpenStreetMap needs no key; Carto's light basemap now watermarks without one.
  m <- leaflet::addProviderTiles(m, "OpenStreetMap.Mapnik", group = "OpenStreetMap")
  m <- leaflet::addProviderTiles(m, "Esri.WorldImagery", group = "Esri imagery")
  groups <- c("Traverse lines", "Tie lines")
  for (ln in split(s, s$line)) {
    d <- ln[seq(1, nrow(ln), by = decimate), ]
    if (nrow(ln) > 1 && !(nrow(ln) %in% seq(1, nrow(ln), by = decimate))) d <- rbind(d, ln[nrow(ln), ])
    tie <- ln$line_type[1] == "tie"
    # Line names come straight from the XYZ file header. leaflet >= 2.0
    # HTML-escapes `label` itself (escaping here again would double-escape);
    # popups are not escaped by leaflet, which is why .popup() escapes (#1).
    # layerId lets the viewer link a click on a line to its profile.
    m <- leaflet::addPolylines(
      m, lng = d$lon, lat = d$lat, weight = 1.6, opacity = 0.85,
      color = if (tie) .pal$tie else .pal$traverse,
      group = if (tie) "Tie lines" else "Traverse lines",
      label = ln$line[1], layerId = ln$line[1])
  }
  for (k in seq_len(nrow(reg))) {
    f <- flags[flags$check == reg$check[k], ]
    if (!nrow(f)) next
    grp <- sprintf("%s (%d)", reg$label[k], nrow(f))
    groups <- c(groups, grp)
    for (r in seq_len(nrow(f))) {
      pop <- .popup(f[r, ])
      if (!is.na(f$i_start[r]) && f$i_end[r] > f$i_start[r]) {
        seg <- s[f$i_start[r]:f$i_end[r], ]
        m <- leaflet::addPolylines(m, lng = seg$lon, lat = seg$lat, weight = 6,
                                   opacity = 0.9, color = .pal$critical,
                                   group = grp, popup = pop)
      } else if (is.finite(f$lon[r])) {
        m <- leaflet::addCircleMarkers(
          m, lng = f$lon[r], lat = f$lat[r],
          radius = if (is.na(f$i_start[r])) 10 else 6,
          stroke = TRUE, color = .pal$surface, weight = 2,
          fillColor = .pal$critical, fillOpacity = 0.9, group = grp, popup = pop)
      }
    }
  }
  m <- leaflet::addLayersControl(m, baseGroups = c("OpenStreetMap", "Esri imagery"),
                                 overlayGroups = groups,
                                 options = leaflet::layersControlOptions(collapsed = FALSE))
  m <- leaflet::addLegend(m, position = "bottomleft",
                          colors = c(.pal$traverse, .pal$tie, .pal$critical),
                          labels = c("traverse line", "tie line", "flagged interval"),
                          opacity = 0.9)
  m <- leaflet::addScaleBar(m, position = "bottomright")
  ok <- is.finite(s$lon) & is.finite(s$lat)
  pad_lon <- diff(range(s$lon[ok])) * 0.08; pad_lat <- diff(range(s$lat[ok])) * 0.08
  leaflet::fitBounds(m, min(s$lon[ok]) - pad_lon, min(s$lat[ok]) - pad_lat,
                     max(s$lon[ok]) + pad_lon, max(s$lat[ok]) + pad_lat)
}

.plotly_map <- function(s, flags, reg, decimate) {
  d <- s[seq(1, nrow(s), by = decimate), ]
  d <- .gap_between_lines(d)
  p <- plotly::plot_ly()
  for (typ in c("traverse", "tie")) {
    dd <- d[d$line_type == typ | is.na(d$line_type), ]
    p <- plotly::add_trace(p, data = dd, x = ~x, y = ~y, type = "scatter", mode = "lines",
                           line = list(color = if (typ == "tie") .pal$tie else .pal$traverse, width = 1.5),
                           text = ~line, hoverinfo = "text", name = paste(typ, "lines"))
  }
  for (k in seq_len(nrow(reg))) {
    f <- flags[flags$check == reg$check[k], ]
    if (!nrow(f)) next
    idx <- unlist(lapply(which(!is.na(f$i_start)), function(r) c(f$i_start[r]:f$i_end[r], NA)))
    if (length(idx)) {
      seg <- s[idx, ]
      p <- plotly::add_trace(p, data = seg, x = ~x, y = ~y, type = "scatter", mode = "lines+markers",
                             line = list(color = .pal$critical, width = 5),
                             marker = list(color = .pal$critical, size = 6),
                             hoverinfo = "text", text = rep(f$description[!is.na(f$i_start)][1], nrow(seg)),
                             name = sprintf("%s (%d)", reg$label[k], nrow(f)))
    }
    ff <- f[is.na(f$i_start) & is.finite(f$x), ]
    if (nrow(ff)) {
      p <- plotly::add_trace(p, data = ff, x = ~x, y = ~y, type = "scatter", mode = "markers",
                             marker = list(color = .pal$critical, size = 14,
                                           line = list(color = .pal$surface, width = 2)),
                             text = ~description, hoverinfo = "text",
                             name = sprintf("%s (%d)", reg$label[k], nrow(f)))
    }
  }
  p <- .plotly_base(p, "northing (m)", "easting (m)", date_axis = FALSE)
  plotly::layout(p, hovermode = "closest", yaxis = list(scaleanchor = "x", scaleratio = 1))
}

# ---- diagnostics ---------------------------------------------------------------

#' Diagnostic plots
#'
#' Each returns a plotly widget for one section of the report.
#' `plot_fourth_difference()` shows the despiked fourth-difference series
#' against its limit; `plot_clearance()` the radar altimeter against the
#' nominal band; `plot_diurnal()` the base-station record with exceedance
#' windows shaded; `plot_crossovers()` the misfit at every traverse/tie
#' intersection, by flight direction.
#'
#' @param result A `magqc_result` from [run_qc()].
#' @name diagnostics
NULL

#' Times for a plotly date axis, rendered in UTC
#'
#' Not epoch milliseconds: plotly.js shows numeric dates in the viewer's local
#' zone, so a report opened in Denver would put the survey six hours away from
#' its own base-station plot. Strings render as written; a tenth of a second
#' is enough for a 10 Hz survey and saves six characters per sample over
#' plotly's default microsecond format.
#' @noRd
.time_str <- function(time) format(time, "%Y-%m-%d %H:%M:%OS1", tz = "UTC")

#' One trace per line, sized for a self-contained report
#'
#' A single trace over the whole survey needs a per-point line name for the
#' hover and NA rows to break the line at turns, and the ISO timestamps and
#' full-precision doubles dominate the report size (#15). One trace per line
#' carries the name once (`%{fullData.name}` in the hover), needs no gap
#' rows, and rounds `y` to what the units warrant. `decimate` thins smooth
#' series; flagged samples are drawn separately at full resolution.
#' @noRd
.add_line_traces <- function(p, s, y, digits, hover_units, decimate = 1L, legend_name) {
  # The legend shows one entry for the series (plotly hides entries for empty
  # traces, so no stub): the first line's trace is named for the series and
  # carries its line name per point for the hover; every other trace is named
  # for its line, hidden from the legend, and hovers via %{fullData.name}.
  first <- TRUE
  for (idx in .by_line(s)) {
    if (decimate > 1L) {
      keep <- unique(c(seq(1L, length(idx), by = decimate), length(idx)))
      idx <- idx[keep]
    }
    if (first) {
      p <- plotly::add_trace(
        p, x = .time_str(s$time[idx]), y = round(y[idx], digits), type = "scatter", mode = "lines",
        line = list(color = .pal$traverse, width = 1), name = legend_name,
        legendgroup = legend_name, showlegend = TRUE, text = s$line[idx],
        hovertemplate = I(paste0("%{text}: %{y:.", digits, "f} ", hover_units, "<extra></extra>")))
      first <- FALSE
    } else {
      p <- plotly::add_trace(
        p, x = .time_str(s$time[idx]), y = round(y[idx], digits), type = "scatter", mode = "lines",
        line = list(color = .pal$traverse, width = 1), name = s$line[idx[1]],
        legendgroup = legend_name, showlegend = FALSE,
        hovertemplate = I(paste0("%{fullData.name}: %{y:.", digits, "f} ", hover_units, "<extra></extra>")))
    }
  }
  p
}

#' Shrink a plotly widget for the self-contained report (#15)
#'
#' Two things: the built widget carries `attrs` and `visdat` - the R-side
#' inputs plotly used to build the traces, a second copy of every data
#' vector that the JS binding never reads - so they are dropped; and the
#' full plotly.js is 3.7 MB of an otherwise ~2 MB report, while the report's
#' traces are all plain `scatter`, which the ~1 MB "basic" bundle covers.
#' `partial_bundle()` downloads that bundle at render time, so an offline
#' render falls back to the full library silently. Set
#' `options(magqc.partial_bundle = FALSE)` to skip the download attempt (the
#' `attrs` stripping still applies).
#' @noRd
.slim <- function(p) {
  if (is.null(p)) return(p)
  p <- plotly::plotly_build(p)
  p$x$attrs <- NULL
  p$x$visdat <- NULL
  # plotly recycles a scalar hovertemplate to one copy per point (I() does not
  # stop it); a constant vector collapses back to the scalar it came from
  for (i in seq_along(p$x$data)) {
    ht <- p$x$data[[i]]$hovertemplate
    if (length(ht) > 1 && length(unique(ht)) == 1) p$x$data[[i]]$hovertemplate <- ht[[1]]
  }
  if (!isTRUE(getOption("magqc.partial_bundle", TRUE))) return(p)
  tryCatch(suppressWarnings(plotly::partial_bundle(p, type = "basic", local = TRUE, minified = TRUE)),
           error = function(e) p)
}

#' @rdname diagnostics
#' @export
plot_fourth_difference <- function(result) {
  s <- result$survey
  tol <- result$spec$fourth_diff_tol
  d4 <- result$fourth_difference
  bad <- .flag_mask(result, "fourth_difference")
  p <- plotly::plot_ly()
  p <- .add_line_traces(p, s, d4, digits = 4, hover_units = "nT", legend_name = "fourth difference")
  if (any(bad)) {
    p <- plotly::add_trace(p, x = .time_str(s$time[bad]), y = round(d4[bad], 4), type = "scatter", mode = "markers",
                           marker = list(color = .pal$critical, size = 5), name = "exceeds limit",
                           text = s$line[bad], hovertemplate = I("%{text}: %{y:.4f} nT<extra></extra>"))
  }
  p <- .plotly_base(p, "fourth difference (nT)")
  plotly::layout(p, shapes = list(
    list(type = "line", xref = "paper", x0 = 0, x1 = 1, y0 = tol, y1 = tol,
         line = list(color = .pal$muted, dash = "dot", width = 1)),
    list(type = "line", xref = "paper", x0 = 0, x1 = 1, y0 = -tol, y1 = -tol,
         line = list(color = .pal$muted, dash = "dot", width = 1))),
    annotations = list(list(xref = "paper", x = 1, y = tol, text = sprintf("limit %g nT", tol),
                            showarrow = FALSE, xanchor = "right", yanchor = "bottom",
                            font = list(color = .pal$muted, size = 11))))
}

#' @rdname diagnostics
#' @export
plot_clearance <- function(result) {
  s <- result$survey; sp <- result$spec
  bad <- .flag_mask(result, "clearance")
  p <- plotly::plot_ly()
  # the altimeter track is smooth at the sample scale; every third sample is
  # plenty for the overview, and flagged samples are drawn in full below
  p <- .add_line_traces(p, s, s$radar_alt, digits = 1, hover_units = "m", decimate = 3L,
                        legend_name = "radar altimeter")
  if (any(bad)) {
    p <- plotly::add_trace(p, x = .time_str(s$time[bad]), y = round(s$radar_alt[bad], 1), type = "scatter", mode = "markers",
                           marker = list(color = .pal$critical, size = 5), name = "outside tolerance",
                           text = s$line[bad], hovertemplate = I("%{text}: %{y:.1f} m<extra></extra>"))
  }
  p <- .plotly_base(p, "terrain clearance (m)")
  plotly::layout(p, shapes = list(
    list(type = "rect", xref = "paper", x0 = 0, x1 = 1,
         y0 = sp$nominal_clearance - sp$clearance_tol, y1 = sp$nominal_clearance + sp$clearance_tol,
         fillcolor = .pal$band, line = list(width = 0), layer = "below")),
    annotations = list(list(xref = "paper", x = 1, y = sp$nominal_clearance + sp$clearance_tol,
                            text = sprintf("%g \u00b1 %g m", sp$nominal_clearance, sp$clearance_tol),
                            showarrow = FALSE, xanchor = "right", yanchor = "bottom",
                            font = list(color = .pal$muted, size = 11))))
}

#' @rdname diagnostics
#' @export
plot_diurnal <- function(result) {
  d <- result$diurnal
  if (is.null(d)) return(NULL)
  p <- plotly::plot_ly(d, x = ~time)
  p <- plotly::add_trace(p, y = ~round(mag_base, 2), type = "scatter", mode = "lines",
                         line = list(color = .pal$traverse, width = 1.5), name = "base station",
                         customdata = ~round(chord_dev, 2),
                         hovertemplate = "%{y:.1f} nT  (chord departure %{customdata:.2f} nT)<extra></extra>")
  win <- .runs(d$exceed)
  shapes <- lapply(seq_len(nrow(win)), function(r) list(
    type = "rect", xref = "x", yref = "paper",
    x0 = d$time[win$start[r]], x1 = d$time[win$end[r]], y0 = 0, y1 = 1,
    fillcolor = .pal$flag_fill, line = list(width = 0), layer = "below"))
  ann <- if (nrow(win)) list(list(
    xref = "x", yref = "paper", x = d$time[win$start[1]], y = 1,
    text = sprintf("exceeds %g nT / %g s", result$spec$max_diurnal_dev, result$spec$diurnal_window),
    showarrow = FALSE, xanchor = "left", yanchor = "top",
    font = list(color = .pal$critical, size = 11))) else list()
  p <- .plotly_base(p, "base-station field (nT)")
  plotly::layout(p, shapes = shapes, annotations = ann)
}

#' @rdname diagnostics
#' @export
plot_crossovers <- function(result) {
  xo <- result$crossovers
  if (!nrow(xo)) return(NULL)
  sp <- result$spec
  hd <- result$heading
  s <- result$survey
  tr <- s[s$line_type == "traverse", ]
  hdg <- dplyr::summarise(dplyr::group_by(tr, .data$line),
                          bearing = .bearing(dplyr::last(.data$x) - dplyr::first(.data$x),
                                             dplyr::last(.data$y) - dplyr::first(.data$y)),
                          .groups = "drop")
  xo$heading <- paste0(.compass(hdg$bearing[match(xo$traverse, hdg$line)]), "-bound")
  xo <- xo[order(xo$traverse, xo$tie), ]
  xo$id <- seq_len(nrow(xo))
  xo$label <- paste(xo$traverse, "x", xo$tie)
  hs <- sort(unique(xo$heading))
  cols <- .series_colours(hs)
  lev <- result$levelling
  if (!is.null(lev)) {
    key <- paste(lev$crossovers$traverse, lev$crossovers$tie)
    xo$misfit_after <- lev$crossovers$misfit_after[match(paste(xo$traverse, xo$tie), key)]
  }
  p <- plotly::plot_ly()
  for (h in hs) {
    dd <- xo[xo$heading == h, ]
    p <- plotly::add_trace(p, data = dd, x = ~id, y = ~misfit, type = "scatter", mode = "markers",
                           marker = list(color = cols[[h]], size = 9,
                                         line = list(color = .pal$surface, width = 1.5)),
                           name = h, text = ~label,
                           hovertemplate = "%{text}: %{y:+.2f} nT<extra></extra>")
    if (!is.null(lev)) {
      p <- plotly::add_trace(p, data = dd, x = ~id, y = ~misfit_after, type = "scatter", mode = "markers",
                             marker = list(color = cols[[h]], size = 8, symbol = "circle-open",
                                           line = list(width = 1.5)),
                             name = paste(h, "after levelling"), text = ~label,
                             hovertemplate = "%{text}: %{y:+.2f} nT after levelling<extra></extra>")
    }
  }
  # the scorecard is the single source of truth for the statistic (#6)
  sc <- result$scorecard
  rms <- sc$metric[sc$check == "crossovers"]
  if (!length(rms) || is.na(rms)) rms <- sqrt(mean(xo$misfit^2, na.rm = TRUE))
  rms_txt <- if (is.null(lev)) sprintf("RMS %.2f nT (limit %g)", rms, sp$max_crossover_rms) else
    sprintf("RMS %.2f nT before (limit %g), %.2f nT after %s levelling (limit %g)",
            rms, sp$max_crossover_rms, lev$rms[["after"]], lev$order, sp$max_levelled_rms)
  p <- .plotly_base(p, "traverse - tie misfit (nT)", "crossover (ordered by traverse, tie)", date_axis = FALSE)
  plotly::layout(
    p, hovermode = "closest",
    shapes = list(
      list(type = "line", xref = "paper", x0 = 0, x1 = 1, y0 = 0, y1 = 0,
           line = list(color = .pal$axis, width = 1)),
      list(type = "line", xref = "paper", x0 = 0, x1 = 1, y0 = sp$max_crossover_abs,
           y1 = sp$max_crossover_abs, line = list(color = .pal$muted, dash = "dot", width = 1)),
      list(type = "line", xref = "paper", x0 = 0, x1 = 1, y0 = -sp$max_crossover_abs,
           y1 = -sp$max_crossover_abs, line = list(color = .pal$muted, dash = "dot", width = 1))),
    annotations = list(list(
      xref = "paper", yref = "paper", x = 1, y = 1.02, showarrow = FALSE,
      xanchor = "right", yanchor = "bottom", font = list(color = .pal$ink2, size = 12),
      text = rms_txt)))
}

# ---- per-line profile ------------------------------------------------------------

#' Profile of one flight line
#'
#' Three stacked panels against distance along the line: the
#' diurnal-corrected field, its despiked fourth difference with the
#' acceptance limit, and the radar altimeter with its tolerance band. Every
#' interval flagged on the line is shaded across all three panels and
#' labelled with the check that raised it. This is the drill-down view the
#' interactive viewer opens when a line is clicked on the map.
#'
#' @param result A `magqc_result` from [run_qc()].
#' @param line Name of the line to draw.
#' @return A plotly widget.
#' @export
plot_line_profile <- function(result, line) {
  s <- result$survey
  idx <- which(s$line == line)
  if (!length(idx)) stop("Line '", line, "' is not in the survey.", call. = FALSE)
  d <- s[idx, ]
  sp <- result$spec
  along <- cumsum(c(0, sqrt(diff(d$x)^2 + diff(d$y)^2)))
  d4 <- result$fourth_difference[idx]
  reg <- .check_registry()
  f <- result$flags[!is.na(result$flags$line) & result$flags$line == line &
                      !is.na(result$flags$i_start), ]
  hover <- sprintf("fid %s, %s", d$fid, format(d$time, "%H:%M:%S"))

  p1 <- plotly::plot_ly(x = along, y = d$mag, type = "scattergl", mode = "lines",
                        line = list(color = .pal$traverse, width = 1.2), text = hover,
                        hovertemplate = "%{y:.2f} nT<br>%{text}<extra></extra>", name = "field")
  if ("mag_lev" %in% names(d) && any(is.finite(d$mag_lev)) && !isTRUE(all.equal(d$mag_lev, d$mag))) {
    p1 <- plotly::add_trace(p1, x = along, y = d$mag_lev, type = "scattergl", mode = "lines",
                            line = list(color = .pal$tie, width = 1.2, dash = "dash"),
                            hovertemplate = "%{y:.2f} nT levelled<extra></extra>",
                            name = "levelled", inherit = FALSE)
  }
  p2 <- plotly::plot_ly(x = along, y = d4, type = "scattergl", mode = "lines",
                        line = list(color = .pal$traverse, width = 1), text = hover,
                        hovertemplate = "%{y:.4f} nT<br>%{text}<extra></extra>", name = "4th diff")
  for (lim in c(-1, 1) * sp$fourth_diff_tol) {
    p2 <- plotly::add_trace(p2, x = range(along), y = c(lim, lim), type = "scatter", mode = "lines",
                            line = list(color = .pal$muted, dash = "dot", width = 1),
                            hoverinfo = "skip", name = "limit", inherit = FALSE)
  }
  p3 <- plotly::plot_ly(x = along, y = d$radar_alt, type = "scattergl", mode = "lines",
                        line = list(color = .pal$traverse, width = 1.2), text = hover,
                        hovertemplate = "%{y:.1f} m<br>%{text}<extra></extra>", name = "clearance")
  for (lim in sp$nominal_clearance + c(-1, 1) * sp$clearance_tol) {
    p3 <- plotly::add_trace(p3, x = range(along), y = c(lim, lim), type = "scatter", mode = "lines",
                            line = list(color = .pal$muted, dash = "dot", width = 1),
                            hoverinfo = "skip", name = "band", inherit = FALSE)
  }

  min_w <- diff(range(along)) * 0.004
  shapes <- lapply(seq_len(nrow(f)), function(r) {
    x0 <- along[f$i_start[r] - idx[1] + 1]; x1 <- along[f$i_end[r] - idx[1] + 1]
    if (x1 - x0 < min_w) { x0 <- x0 - min_w / 2; x1 <- x1 + min_w / 2 }
    list(type = "rect", xref = "x", yref = "paper", x0 = x0, x1 = x1, y0 = 0, y1 = 1,
         fillcolor = .pal$flag_fill, line = list(width = 0), layer = "below")
  })
  ann <- if (nrow(f) && nrow(f) <= 15) lapply(seq_len(nrow(f)), function(r) list(
    xref = "x", yref = "paper", y = 1.0, yanchor = "bottom", showarrow = FALSE,
    x = (along[f$i_start[r] - idx[1] + 1] + along[f$i_end[r] - idx[1] + 1]) / 2,
    text = reg$label[match(f$check[r], reg$check)],
    font = list(size = 10, color = .pal$critical))) else list()

  ax <- function(title) list(title = title, gridcolor = .pal$grid, zeroline = FALSE,
                             linecolor = .pal$axis, tickfont = list(color = .pal$muted))
  p <- plotly::subplot(p1, p2, p3, nrows = 3, shareX = TRUE, titleY = TRUE,
                       heights = c(0.42, 0.29, 0.29))
  plotly::layout(
    p, showlegend = FALSE, shapes = shapes, annotations = ann,
    font = list(family = .font, color = .pal$ink2, size = 12),
    paper_bgcolor = .pal$surface, plot_bgcolor = .pal$surface,
    margin = list(l = 60, r = 20, t = 30, b = 45), hovermode = "x unified",
    xaxis = ax("distance along line (m)"),
    yaxis = ax("field (nT)"), yaxis2 = ax("4th diff (nT)"), yaxis3 = ax("clearance (m)"))
}

#' Flag table for display
#'
#' The flag tibble from [run_qc()] reshaped for a reader: check labels
#' instead of slugs, times and fiducial ranges as strings, and rounded
#' values. Used by the report and the viewer.
#'
#' @param result A `magqc_result` from [run_qc()].
#' @return A data frame with columns `check`, `line`, `start`, `end`, `fid`,
#'   `value`, `limit`, `units`, `finding`.
#' @export
flag_table <- function(result) {
  f <- result$flags
  lab <- .check_registry()
  if (!nrow(f)) {
    return(data.frame(check = character(), line = character(), start = character(),
                      end = character(), fid = character(), value = numeric(),
                      limit = numeric(), units = character(), finding = character(),
                      stringsAsFactors = FALSE))
  }
  data.frame(
    check = lab$label[match(f$check, lab$check)],
    line = ifelse(is.na(f$line), "survey-wide", f$line),
    start = ifelse(is.na(f$time_start), "", format(f$time_start, "%H:%M:%S")),
    end = ifelse(is.na(f$time_end), "", format(f$time_end, "%H:%M:%S")),
    fid = ifelse(is.na(f$fid_start), "",
                 ifelse(f$fid_end == f$fid_start, as.character(f$fid_start),
                        paste0(f$fid_start, "-", f$fid_end))),
    value = signif(f$value, 3), limit = signif(f$threshold, 3), units = f$units,
    finding = f$description, stringsAsFactors = FALSE)
}

#' Registry of QC checks
#'
#' @return A tibble with the check slug and its display label, in report
#'   order.
#' @export
check_registry <- function() .check_registry()

# ---- report ------------------------------------------------------------------------

#' Render the HTML acceptance report
#'
#' Renders the bundled R Markdown template to a single self-contained HTML
#' file containing the summary statistics, the pass/fail scorecard, the map
#' of flight lines with every flagged interval, the flag table, diagnostic
#' plots and the specification the survey was tested against.
#'
#' @param result A `magqc_result` from [run_qc()].
#' @param output Path of the HTML file to write.
#' @param title Report title.
#' @param open Open the report in the browser afterwards?
#' @return The path to the rendered file, invisibly.
#' @examples
#' \dontrun{
#' res <- run_qc(sim_survey())
#' qc_report(res, "magqc-report.html")
#' }
#' @export
qc_report <- function(result, output = "magqc-report.html",
                      title = "Airborne magnetic survey QA/QC", open = FALSE) {
  stopifnot(inherits(result, "magqc_result"))
  tmpl <- system.file("rmarkdown", "qc_report.Rmd", package = "magqc")
  if (!nzchar(tmpl)) stop("Report template not found; is magqc installed?", call. = FALSE)
  output <- normalizePath(output, mustWork = FALSE)
  # render in a scratch copy so intermediate files never land beside the user's data
  work <- tempfile("magqc-report-")
  dir.create(work)
  file.copy(tmpl, file.path(work, "qc_report.Rmd"))
  rmarkdown::render(
    file.path(work, "qc_report.Rmd"),
    output_file = output, output_dir = dirname(output),
    params = list(result = result, title = title),
    envir = new.env(parent = globalenv()), quiet = TRUE)
  if (isTRUE(open)) utils::browseURL(output)
  invisible(output)
}

#' @importFrom rlang .data
NULL

# DT and knitr are called from the R Markdown template, not from package code;
# referencing them here keeps R CMD check from reporting them as unused.
.template_deps <- function() list(DT::datatable, knitr::knit)
