# Interactive viewer for the package.
#
# Runs two ways. Installed: `run_viewer()` launches this file with the
# package loaded. Static: tools/build-shinylive.R vendors the package's R/
# sources beside this file and exports a shinylive (webR) site, in which case
# the sources are sourced directly. The report download is only offered when
# pandoc is available, which it never is inside the browser.
#
# shinylive scans the app for package references and installs each one from
# the webR repository, and this package is not there - so the build script
# strips the package name from the vendored copy (the `pkg` line below is the
# only literal use) and the vendored branch never mentions it.
vendored <- file.exists("R/spec.R")
if (vendored) {
  library(dplyr)
  library(rlang)
  for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)
} else {
  pkg <- "magqc"
  do.call(library, list(pkg, character.only = TRUE))
}
library(shiny)
library(bslib)

can_render <- !vendored && requireNamespace("rmarkdown", quietly = TRUE) &&
  rmarkdown::pandoc_available("2.0")
defaults <- survey_spec()

spec_fields <- c(
  line_spacing = "Traverse spacing (m)", line_spacing_tol = "Spacing tolerance (fraction)",
  line_azimuth = "Traverse azimuth (deg)", sample_rate = "Sample rate (Hz)",
  ground_speed = "Ground speed (m/s)", max_deviation = "Max departure from line (m)",
  nominal_clearance = "Nominal clearance (m)", clearance_tol = "Clearance tolerance (m)",
  fourth_diff_tol = "4th-difference limit (nT)", spike_nsigma = "Spike threshold (sigma)",
  max_diurnal_dev = "Diurnal limit (nT)", diurnal_window = "Diurnal chord (s)",
  max_crossover_rms = "Crossover RMS limit (nT)", max_crossover_abs = "Single crossover limit (nT)",
  max_heading_error = "Heading error limit (nT)", max_gap = "Max data gap (s)")

spec_inputs <- lapply(names(spec_fields), function(f) {
  numericInput(f, spec_fields[[f]], value = defaults[[f]], width = "100%")
})

ui <- page_sidebar(
  title = "magqc — airborne magnetic survey QA/QC",
  theme = bs_theme(version = 5, base_font = "system-ui, -apple-system, 'Segoe UI', sans-serif",
                   primary = "#2a78d6"),
  sidebar = sidebar(
    width = 340,
    accordion(
      open = "Data",
      accordion_panel(
        "Data",
        radioButtons("source", NULL, c("Simulated survey" = "sim", "Upload XYZ" = "upload"),
                     inline = TRUE),
        conditionalPanel(
          "input.source == 'sim'",
          layout_columns(
            numericInput("seed", "Seed", 42, width = "100%"),
            numericInput("n_lines", "Traverse lines", 14, min = 12, step = 1, width = "100%")),
          checkboxInput("defects", "Inject defects", TRUE)),
        conditionalPanel(
          "input.source == 'upload'",
          fileInput("xyz", "XYZ flight file", accept = c(".xyz", ".txt", ".dat", ".csv")),
          layout_columns(
            textInput("col_x", "x column", "X"), textInput("col_y", "y column", "Y")),
          layout_columns(
            textInput("col_time", "time column", "TIME"), textInput("col_mag", "field column", "MAG")),
          textInput("col_alt", "radar altimeter column (blank if none)", "RALT"),
          dateInput("date", "Flight date (UTC; time column is seconds since midnight)"),
          layout_columns(
            numericInput("lon0", "Origin lon (optional)", NA, width = "100%"),
            numericInput("lat0", "Origin lat (optional)", NA, width = "100%")),
          fileInput("base", "Base station CSV (time, mag_base)", accept = c(".csv", ".txt")))),
      accordion_panel("Specification", spec_inputs)),
    actionButton("run", "Run QC", class = "btn-primary w-100"),
    if (can_render) downloadButton("report", "Download HTML report", class = "w-100") else
      helpText("HTML report download needs pandoc — run the viewer locally for it."),
    downloadButton("flags_csv", "Download flags (CSV)", class = "w-100")
  ),
  layout_columns(
    fill = FALSE,
    value_box("Line-km flown", textOutput("vb_km"), theme = "primary"),
    value_box("Line-km flagged", textOutput("vb_flagged"), theme = "secondary"),
    value_box("Accepted", textOutput("vb_pct"), theme = "secondary"),
    value_box("Checks passed", textOutput("vb_checks"), theme = "secondary")),
  navset_card_tab(
    nav_panel("Scorecard", DT::DTOutput("scorecard")),
    nav_panel(
      "Map & line profile",
      layout_columns(
        col_widths = c(5, 7),
        card(card_header("Flight lines — click a line to open its profile"),
             uiOutput("map_ui"), full_screen = TRUE),
        card(card_header(selectInput("line", NULL, choices = character(0), width = "220px")),
             plotly::plotlyOutput("profile", height = "560px"), full_screen = TRUE))),
    nav_panel("Flags", DT::DTOutput("flags")),
    nav_panel(
      "Diagnostics",
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("Fourth-difference noise"), plotly::plotlyOutput("d4", height = "320px")),
        card(card_header("Terrain clearance"), plotly::plotlyOutput("clearance", height = "320px")),
        card(card_header("Base station"), plotly::plotlyOutput("diurnal", height = "320px")),
        card(card_header("Crossover misfit"), plotly::plotlyOutput("crossovers", height = "320px")))),
    nav_panel("Lines", DT::DTOutput("lines")))
)

server <- function(input, output, session) {
  preset <- getShinyOption("magqc_result", NULL)

  spec <- reactive({
    args <- lapply(names(spec_fields), function(f) input[[f]])
    names(args) <- names(spec_fields)
    args <- Filter(function(v) !is.null(v) && !is.na(v), args)
    do.call(survey_spec, args)
  })

  load_upload <- function() {
    req(input$xyz)
    col_map <- c(x = input$col_x, y = input$col_y, time = input$col_time, mag_raw = input$col_mag)
    if (nzchar(trimws(input$col_alt))) col_map <- c(col_map, radar_alt = trimws(input$col_alt))
    origin <- if (is.finite(input$lon0) && is.finite(input$lat0)) c(lon = input$lon0, lat = input$lat0)
    survey <- read_xyz(input$xyz$datapath, col_map = col_map,
                       time_origin = as.POSIXct(input$date, tz = "UTC"), origin = origin)
    base <- NULL
    if (!is.null(input$base)) {
      base <- utils::read.csv(input$base$datapath, stringsAsFactors = FALSE)
      if (!all(c("time", "mag_base") %in% names(base))) {
        stop("Base station CSV needs `time` and `mag_base` columns.")
      }
      base$time <- as.POSIXct(base$time, tz = "UTC")
    }
    list(survey = survey, base = base)
  }

  res <- eventReactive(input$run, ignoreNULL = FALSE, {
    if (!is.null(preset) && input$run == 0) return(preset)
    withProgress(message = "Running QC checks…", value = 0.3, {
      tryCatch({
        if (input$source == "sim") {
          sim <- sim_survey(seed = input$seed, n_lines = max(12, input$n_lines),
                            defects = isTRUE(input$defects))
          run_qc(sim, spec = spec())
        } else {
          up <- load_upload()
          run_qc(up$survey, spec = spec(), base = up$base)
        }
      }, error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = 10)
        NULL
      })
    })
  })

  fmt <- function(x, d = 1) formatC(x, format = "f", digits = d, big.mark = ",")
  output$vb_km <- renderText({ req(res()); fmt(res()$stats$line_km) })
  output$vb_flagged <- renderText({ req(res()); fmt(res()$stats$flagged_km, 2) })
  output$vb_pct <- renderText({ req(res()); paste0(fmt(res()$stats$accepted_pct), "%") })
  output$vb_checks <- renderText({
    req(res()); sprintf("%d / %d", res()$stats$checks_passed, res()$stats$checks_total)
  })

  output$scorecard <- DT::renderDT({
    req(res())
    sc <- res()$scorecard
    tbl <- data.frame(
      check = sc$label, status = toupper(sc$status), statistic = sc$metric_label,
      value = ifelse(is.na(sc$metric), NA, signif(sc$metric, 3)), limit = signif(sc$threshold, 3),
      flags = sc$n_flags, `samples affected` = sc$n_samples,
      note = ifelse(is.na(sc$note), "", sc$note), check.names = FALSE)
    DT::datatable(tbl, rownames = FALSE, options = list(dom = "t", pageLength = 20)) |>
      DT::formatStyle("status", fontWeight = "bold",
                      color = DT::styleEqual(c("PASS", "FAIL", "N/A"), c("#006300", "#d03b3b", "#898781")))
  })

  output$map_ui <- renderUI({
    req(res())
    if (isTRUE(res()$stats$has_lonlat)) leaflet::leafletOutput("map", height = "560px") else
      plotly::plotlyOutput("map_xy", height = "560px")
  })
  output$map <- leaflet::renderLeaflet({ req(res(), isTRUE(res()$stats$has_lonlat)); qc_map(res()) })
  output$map_xy <- plotly::renderPlotly({ req(res(), !isTRUE(res()$stats$has_lonlat)); qc_map(res()) })

  observeEvent(res(), {
    lines <- res()$lines$line
    updateSelectInput(session, "line", choices = lines,
                      selected = if (isTRUE(input$line %in% lines)) input$line else lines[1])
  })
  observeEvent(input$map_shape_click, {
    id <- input$map_shape_click$id
    if (!is.null(id) && id %in% res()$lines$line) updateSelectInput(session, "line", selected = id)
  })
  output$profile <- plotly::renderPlotly({
    req(res(), input$line, input$line %in% res()$lines$line)
    plot_line_profile(res(), input$line)
  })

  output$flags <- DT::renderDT({
    req(res())
    DT::datatable(flag_table(res()), rownames = FALSE, filter = "top", class = "compact stripe",
                  options = list(pageLength = 15, dom = "tip", autoWidth = FALSE))
  })
  output$lines <- DT::renderDT({
    req(res())
    ln <- res()$lines
    DT::datatable(data.frame(
      line = ln$line, type = ln$line_type, samples = ln$n_samples,
      start = format(ln$start, "%H:%M:%S"), end = format(ln$end, "%H:%M:%S"),
      km = round(ln$length_km, 2), flagged_km = round(ln$flagged_km, 3),
      flags = ln$n_flags, status = ln$status),
      rownames = FALSE, class = "compact stripe", options = list(pageLength = 25, dom = "tip"))
  })

  output$d4 <- plotly::renderPlotly({ req(res()); plot_fourth_difference(res()) })
  output$clearance <- plotly::renderPlotly({
    req(res()); validate(need(!all(is.na(res()$survey$radar_alt)), "No radar altimeter column."))
    plot_clearance(res())
  })
  output$diurnal <- plotly::renderPlotly({
    req(res()); p <- plot_diurnal(res())
    validate(need(!is.null(p), "No base-station record supplied.")); p
  })
  output$crossovers <- plotly::renderPlotly({
    req(res()); p <- plot_crossovers(res())
    validate(need(!is.null(p), "No traverse/tie crossovers found.")); p
  })

  output$flags_csv <- downloadHandler(
    filename = function() "magqc-flags.csv",
    content = function(file) utils::write.csv(flag_table(res()), file, row.names = FALSE))
  if (can_render) {
    output$report <- downloadHandler(
      filename = function() "magqc-report.html",
      content = function(file) {
        withProgress(message = "Rendering report…", value = 0.5, qc_report(res(), file))
      })
  }
}

shinyApp(ui, server)
