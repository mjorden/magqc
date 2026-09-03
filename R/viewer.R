#' Launch the interactive viewer
#'
#' A Shiny application for exploring a survey: load the simulated survey or
#' upload an XYZ flight file and base-station record, edit the acceptance
#' specification in the sidebar, re-run the checks, click a flight line on
#' the map to open its profile, browse the flag table, and download the
#' report or the flags as CSV.
#'
#' The same application is exported as a static shinylive site (webR) by
#' `tools/build-shinylive.R`, so it can be hosted on GitHub Pages with no
#' server; in that build the HTML report download is unavailable because
#' pandoc does not run in the browser.
#'
#' @param result Optional `magqc_result` to open with. When `NULL` the app
#'   starts by running the simulated survey.
#' @param ... Passed to [shiny::runApp()], e.g. `launch.browser = FALSE`.
#' @return Called for its side effect.
#' @examples
#' \dontrun{
#' run_viewer()
#' run_viewer(run_qc(sim_survey(seed = 3)))
#' }
#' @export
run_viewer <- function(result = NULL, ...) {
  rlang::check_installed(c("shiny", "bslib", "DT"), reason = "to run the viewer")
  if (!is.null(result) && !inherits(result, "magqc_result")) {
    stop("`result` must be a magqc_result from run_qc().", call. = FALSE)
  }
  app_dir <- system.file("shiny", package = "magqc")
  if (!nzchar(app_dir)) stop("Viewer app not found; is magqc installed?", call. = FALSE)
  shiny::shinyOptions(magqc_result = result)
  shiny::runApp(app_dir, ...)
}
