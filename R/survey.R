#' Build a validated survey table
#'
#' Coerces a data frame of fiducial-level readings into the shape every check
#' expects and fails loudly if a required column is missing. Rows are sorted
#' by line and then time, and a row index `.i` is attached so that flagged
#' intervals can be mapped back to the exact samples they came from.
#'
#' Required columns: `line`, `time` (POSIXct in any time zone - the instants
#' are kept and re-stamped to render in UTC, which is how every report and
#' flag time is labelled), `x`, `y` (metres, any local or projected system),
#' and one of `mag_raw` or `mag` (nT). Optional:
#' `line_type` (`"traverse"` or `"tie"`; inferred from the line name when
#' absent - names beginning with `T` are ties), `fid`, `radar_alt` (m AGL),
#' `gps_alt` (m), `lon`, `lat`.
#'
#' @param df A data frame.
#' @param origin Optional `c(lon, lat)` used to derive `lon`/`lat` from
#'   `x`/`y` when they are not supplied. Without it the map in the report is
#'   drawn in local metres.
#' @return A tibble with class `magqc_survey`.
#' @export
as_survey <- function(df, origin = NULL) {
  df <- tibble::as_tibble(df)
  req <- c("line", "time", "x", "y")
  miss <- setdiff(req, names(df))
  if (length(miss)) {
    stop("Survey is missing required column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  }
  if (!("mag_raw" %in% names(df)) && !("mag" %in% names(df))) {
    stop("Survey needs a `mag_raw` or `mag` column (total field, nT).", call. = FALSE)
  }
  if (!inherits(df$time, "POSIXct")) {
    stop("`time` must be POSIXct.", call. = FALSE)
  }
  # POSIXct instants are absolute; only their rendering depends on the
  # tzone attribute. Everything downstream labels times "UTC", so render
  # them that way regardless of the zone the caller built them in (#8).
  attr(df$time, "tzone") <- "UTC"
  if (!"line_type" %in% names(df)) {
    df$line_type <- ifelse(grepl("^T", df$line, ignore.case = TRUE), "tie", "traverse")
  }
  bad_type <- setdiff(unique(df$line_type), c("traverse", "tie"))
  if (length(bad_type)) {
    stop("`line_type` must be 'traverse' or 'tie'; found: ", paste(bad_type, collapse = ", "), call. = FALSE)
  }
  df$line <- as.character(df$line)
  if (!"fid" %in% names(df)) df$fid <- seq_len(nrow(df))
  if (!"mag" %in% names(df)) df$mag <- df$mag_raw
  if (!"mag_raw" %in% names(df)) df$mag_raw <- df$mag
  if (!"radar_alt" %in% names(df)) df$radar_alt <- NA_real_
  if (!all(c("lon", "lat") %in% names(df))) {
    if (!is.null(origin)) {
      ll <- enu_to_lonlat(df$x, df$y, origin[["lon"]], origin[["lat"]])
      df$lon <- ll$lon; df$lat <- ll$lat
    } else {
      df$lon <- NA_real_; df$lat <- NA_real_
    }
  }
  df <- dplyr::arrange(df, .data$line_type == "tie", .data$line, .data$time)
  df$.i <- seq_len(nrow(df))
  class(df) <- c("magqc_survey", class(df))
  df
}

#' Read a Geosoft-style XYZ flight file
#'
#' Parses the ASCII XYZ interchange format used by most airborne contractors:
#' comment lines begin with `/`, each flight line is introduced by a
#' `Line <name>` or `Tie <name>` header, and the samples that follow are
#' whitespace-delimited. Column names are taken from the last comment line
#' before the first data row unless `columns` is given.
#'
#' @param path Path to the file.
#' @param columns Character vector of column names, in file order. Overrides
#'   any header comment.
#' @param col_map Named character vector mapping the survey columns this
#'   package needs (`x`, `y`, `mag_raw`, `radar_alt`, `time`, ...) to the names
#'   used in the file, e.g. `c(x = "X", y = "Y", mag_raw = "MAGRAW")`.
#' @param time_origin POSIXct date the file's time-of-day column is measured
#'   from. XYZ files usually store seconds since midnight; if the mapped
#'   `time` column is already POSIXct this is ignored.
#' @param ... Passed to [as_survey()].
#' @return A `magqc_survey` tibble.
#' @export
read_xyz <- function(path, columns = NULL, col_map = NULL,
                     time_origin = as.POSIXct("1970-01-01", tz = "UTC"), ...) {
  raw <- readLines(path, warn = FALSE)
  raw <- raw[nzchar(trimws(raw))]
  is_comment <- startsWith(trimws(raw), "/")
  is_header  <- grepl("^\\s*(line|tie)\\s+\\S+", raw, ignore.case = TRUE)
  first_data <- which(!is_comment & !is_header)[1]
  if (is.na(first_data)) stop("No data rows found in ", path, call. = FALSE)

  if (is.null(columns)) {
    hdr <- raw[is_comment & seq_along(raw) < first_data]
    if (!length(hdr)) stop("No header comment found; pass `columns`.", call. = FALSE)
    columns <- strsplit(trimws(sub("^/+", "", hdr[length(hdr)])), "\\s+")[[1]]
  }

  cur_line <- NA_character_; cur_type <- NA_character_
  line_of <- character(length(raw)); type_of <- character(length(raw))
  for (k in seq_along(raw)) {
    if (is_header[k]) {
      m <- regmatches(raw[k], regexec("^\\s*(line|tie)\\s+(\\S+)", raw[k], ignore.case = TRUE))[[1]]
      cur_type <- if (tolower(m[2]) == "tie") "tie" else "traverse"
      cur_line <- m[3]
    }
    line_of[k] <- cur_line; type_of[k] <- cur_type
  }
  keep <- !is_comment & !is_header
  body <- raw[keep]
  dat <- utils::read.table(text = body, header = FALSE, col.names = columns,
                           na.strings = c("*", "NA", "-99999", "1e+32"),
                           check.names = FALSE)
  dat$line <- line_of[keep]
  dat$line_type <- type_of[keep]

  if (!is.null(col_map)) {
    for (nm in names(col_map)) {
      src <- col_map[[nm]]
      if (!src %in% names(dat)) stop("col_map: column `", src, "` not in file.", call. = FALSE)
      dat[[nm]] <- dat[[src]]
    }
  }
  if ("time" %in% names(dat) && !inherits(dat$time, "POSIXct")) {
    dat$time <- time_origin + as.numeric(dat$time)
  }
  as_survey(dat, ...)
}
