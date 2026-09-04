# Regenerate the README figures in man/figures with headless Chrome.
#
#   Rscript tools/screenshots.R
#
# Needs chromote, magick and a Chrome/Chromium install. The viewer figure
# needs the app running on http://127.0.0.1:8790 (e.g.
# `run_viewer(launch.browser = FALSE, port = 8790)` in another session);
# it is skipped with a note when that port does not answer.
#
# chromote's element screenshots reset the scroll position before clipping,
# so page elements are captured by growing the viewport to the whole document
# and cropping a single full-page capture with magick.

for (p in c("chromote", "magick", "devtools")) if (!requireNamespace(p, quietly = TRUE)) stop("install ", p)
devtools::load_all(quiet = TRUE)
out_dir <- "man/figures"; dir.create(out_dir, showWarnings = FALSE)
scale <- 1.5

report <- file.path(tempdir(), "demo.html")
qc_report(run_qc(sim_survey()), report, title = "Demo survey — airborne magnetic QA/QC")

b <- chromote::ChromoteSession$new(width = 1280, height = 900)
on.exit(b$close())
js <- function(code) b$Runtime$evaluate(code, returnByValue = TRUE)$result$value
save <- function(img, file) {
  path <- file.path(out_dir, file)
  magick::image_write(img, path, format = "png")
  cat(sprintf("%-30s %5.0f kB\n", file, file.size(path) / 1024))
}

# --- report: overview (viewport-sized), then page elements from one tall capture
b$Page$navigate(paste0("file:///", normalizePath(report, winslash = "/"))); Sys.sleep(5)
tmp <- file.path(tempdir(), "shot.png")
b$screenshot(tmp, cliprect = c(0, 0, 1280, 760), scale = scale, show = FALSE)
save(magick::image_read(tmp), "report-overview.png")

h <- as.integer(js("document.documentElement.scrollHeight")) + 50L
b$Emulation$setDeviceMetricsOverride(width = 1280L, height = h, deviceScaleFactor = 1, mobile = FALSE)
Sys.sleep(8)
js("var p = document.querySelectorAll('.js-plotly-plot'); p[0].id = 'shot-d4'; p[3].id = 'shot-xo'")
b$screenshot(tmp, selector = "html", scale = scale, show = FALSE)
full <- magick::image_read(tmp)
crop <- function(sel, file, pad = 8) {
  r <- unlist(js(sprintf("(function(){var r=document.querySelector('%s').getBoundingClientRect();return [r.left,r.top,r.width,r.height];})()", sel)))
  g <- sprintf("%dx%d+%d+%d", round((r[3] + 2 * pad) * scale), round((r[4] + 2 * pad) * scale),
               round((r[1] - pad) * scale), round((r[2] - pad) * scale))
  save(magick::image_crop(full, g), file)
}
crop(".leaflet-container", "report-map.png")
crop("#shot-xo", "report-crossovers.png")
crop("#shot-d4", "report-fourth-difference.png")

# --- viewer: map tab with a line selected
b$Emulation$clearDeviceMetricsOverride()
up <- tryCatch({ con <- url("http://127.0.0.1:8790/"); on.exit(close(con), add = TRUE); readLines(con, n = 1); TRUE },
               error = function(e) FALSE, warning = function(w) FALSE)
if (up) {
  b$Page$navigate("http://127.0.0.1:8790/"); Sys.sleep(12)
  js("document.querySelector('a.nav-link[data-value=\"Map & line profile\"]').click()"); Sys.sleep(5)
  js("Shiny.setInputValue('map_shape_click', {id: 'L1050', lat: 0, lng: 0, '.nonce': Math.random()})"); Sys.sleep(6)
  b$screenshot(tmp, cliprect = c(0, 0, 1280, 900), scale = scale, show = FALSE)
  save(magick::image_read(tmp), "viewer.png")
} else {
  cat("viewer not running on :8790 - skipped viewer.png\n")
}
