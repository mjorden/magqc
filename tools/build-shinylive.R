# Build the static (webR) viewer site.
#
# shinylive runs the app entirely in the browser, installing CRAN packages
# from the webR repository at load time. magqc itself is not on CRAN, so the
# package's R/ sources are vendored beside app.R and sourced by the app (see
# the `vendored` switch at the top of inst/shiny/app.R). Run from the package
# root:
#
#   Rscript tools/build-shinylive.R            # writes ./site
#   Rscript tools/build-shinylive.R --serve    # ...and serves it locally
#
# Requires the `shinylive` package (>= 0.2) and, on first run, network access
# to download the shinylive assets.

args <- commandArgs(trailingOnly = TRUE)
dest <- "site"
if (!requireNamespace("shinylive", quietly = TRUE)) {
  stop("install.packages('shinylive') first.", call. = FALSE)
}
if (!file.exists("DESCRIPTION") || !file.exists("inst/shiny/app.R")) {
  stop("Run from the magqc package root.", call. = FALSE)
}

pkg_name <- read.dcf("DESCRIPTION", fields = "Package")[[1]]
app_src <- file.path(tempdir(), "viewer-site-src")
unlink(app_src, recursive = TRUE)
dir.create(file.path(app_src, "R"), recursive = TRUE)

# shinylive's dependency scan (run inside webR) installs every package it can
# find a reference to, and this package is not in the webR repository, so the
# vendored app.R must not contain the package name at all. The installed-mode
# branch is dead code in the vendored copy; neutralise its one literal.
pkg_ref <- sprintf("(library|require|requireNamespace|loadNamespace)\\(['\"]?%s|%s::", pkg_name, pkg_name)
app <- readLines("inst/shiny/app.R", warn = FALSE)
app <- sub(sprintf('pkg <- "%s"', pkg_name), "pkg <- NULL", app, fixed = TRUE)
if (any(grepl(pkg_ref, app))) {
  stop("inst/shiny/app.R still references the '", pkg_name,
       "' package after sanitising; shinylive would try to install it.", call. = FALSE)
}
writeLines(app, file.path(app_src, "app.R"))

r_files <- list.files("R", pattern = "[.]R$", full.names = TRUE)
stopifnot(all(file.copy(r_files, file.path(app_src, "R"))))
# The same rule applies to the vendored sources: `pkg::fun` references in
# roxygen examples/comments would be picked up too.
for (f in file.path(app_src, "R", basename(r_files))) {
  src <- readLines(f, warn = FALSE)
  src <- gsub(paste0(pkg_name, ":::?"), "", src, fixed = FALSE)
  # system.file(..., package = "<pkg>") is also a scanner hit; dropping the
  # argument makes it look in base, which returns "" - the right answer for
  # code paths (report template, app dir) that do not exist in the browser.
  src <- gsub(paste0(', package = "', pkg_name, '"'), "", src, fixed = TRUE)
  writeLines(src, f)
}
leftover <- unlist(lapply(file.path(app_src, "R", basename(r_files)),
                          function(f) grep(pkg_ref, readLines(f, warn = FALSE), value = TRUE)))
if (length(leftover)) {
  stop("vendored sources still reference the '", pkg_name, "' package:\n  ",
       paste(leftover, collapse = "\n  "), call. = FALSE)
}

unlink(dest, recursive = TRUE)
shinylive::export(app_src, dest, quiet = FALSE)

# A first visit downloads the webR runtime and the packages (a minute or two);
# shinylive shows only a spinner, so add a note that removes itself once the
# app's iframe has rendered.
index <- file.path(dest, "index.html")
html <- readLines(index, warn = FALSE)
splash <- c(
  '<div id="viewer-splash" style="position:fixed;left:0;right:0;top:58%;text-align:center;',
  'font:14px system-ui,-apple-system,\'Segoe UI\',sans-serif;color:#52514e;z-index:0;pointer-events:none">',
  'Loading R in your browser (webR) &mdash; the first visit downloads the runtime and packages and can take a minute or two.',
  '</div>',
  '<script>(function(){var t=setInterval(function(){var f=document.querySelector("iframe");',
  'var d=f&&f.contentDocument;if(d&&d.querySelector(".bslib-value-box")){var s=document.getElementById("viewer-splash");',
  'if(s)s.remove();clearInterval(t);}},500);})();</script>')
body_end <- grep("</body>", html, fixed = TRUE)
stopifnot(length(body_end) == 1)
html <- append(html, splash, after = body_end - 1)
writeLines(html, index)
cat("Static viewer written to", normalizePath(dest), "\n")

if ("--serve" %in% args) {
  httpuv::runStaticServer(dest, port = 8000)
}
