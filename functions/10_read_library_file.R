# ==============================================================================
# Read and standardize sgRNA library file
# ==============================================================================

read_library_file <- function(library_path,
                              default_sublib = "ICS_1",
                              default_count = 1) {
  # ---------------------------------------------------------------------------
  # Dependencies
  # ---------------------------------------------------------------------------
  
  if (!requireNamespace("tools", quietly = TRUE)) {
    stop("The 'tools' package is required.")
  }
  
  if (!file.exists(library_path)) {
    stop("Library file does not exist: ", library_path)
  }
  
  # ---------------------------------------------------------------------------
  # Helpers
  # ---------------------------------------------------------------------------
  
  normalize_colnames <- function(x) {
    x <- trimws(x)
    x <- gsub("\\s+", "_", x)
    x <- gsub("\\.+", "_", x)
    x <- tolower(x)
    x
  }
  
  find_column <- function(df, canonical_name, alternatives, required = TRUE) {
    possible_names <- c(canonical_name, alternatives)
    possible_names <- normalize_colnames(possible_names)
    
    matches <- intersect(possible_names, colnames(df))
    
    if (length(matches) == 0) {
      if (required) {
        stop(
          "Library is missing required column '", canonical_name, "'. ",
          "Accepted names are: ",
          paste(possible_names, collapse = ", ")
        )
      } else {
        return(NULL)
      }
    }
    
    if (length(matches) > 1) {
      warning(
        "Multiple possible columns found for '", canonical_name, "': ",
        paste(matches, collapse = ", "),
        ". Using '", matches[1], "'.",
        call. = FALSE
      )
    }
    
    matches[1]
  }
  
  is_missing_string <- function(x) {
    is.na(x) | trimws(as.character(x)) == ""
  }
  
  make_unique_with_warning <- function(x) {
    x <- as.character(x)
    
    if (!anyDuplicated(x)) {
      return(x)
    }
    
    warning(
      "Column 'sgrna_id' contains duplicated values. ",
      "Making them unique by adding '_#' suffixes.",
      call. = FALSE
    )
    
    out <- x
    seen <- new.env(parent = emptyenv())
    
    for (i in seq_along(x)) {
      value <- x[i]
      
      if (!exists(value, envir = seen, inherits = FALSE)) {
        assign(value, 1L, envir = seen)
        out[i] <- value
      } else {
        n <- get(value, envir = seen, inherits = FALSE) + 1L
        assign(value, n, envir = seen)
        out[i] <- paste0(value, "_", n)
      }
    }
    
    out
  }
  
  coerce_integerish <- function(x, column_name, allow_na = TRUE) {
    original <- x
    
    x <- as.character(x)
    x <- trimws(x)
    
    missing <- is.na(original) | x == ""
    
    if (!allow_na && any(missing)) {
      stop("Column '", column_name, "' contains NA or empty values.")
    }
    
    out <- rep(NA_integer_, length(x))
    
    non_missing <- !missing
    
    valid <- grepl("^[0-9]+$", x[non_missing])
    
    if (any(!valid)) {
      bad_values <- unique(x[non_missing][!valid])
      stop(
        "Column '", column_name, "' must contain integers or strings with numbers only. ",
        "Invalid value(s): ",
        paste(head(bad_values, 10), collapse = ", ")
      )
    }
    
    out[non_missing] <- as.integer(x[non_missing])
    
    out
  }
  
  validate_no_na <- function(x, column_name) {
    if (any(is_missing_string(x))) {
      stop("Column '", column_name, "' contains NA or empty values, but none are allowed.")
    }
  }
  
  validate_nucleotide_sequence <- function(x, column_name) {
    validate_no_na(x, column_name)
    
    x <- toupper(trimws(as.character(x)))
    
    invalid <- !grepl("^[ACGTN]+$", x)
    
    if (any(invalid)) {
      bad_values <- unique(x[invalid])
      stop(
        "Column '", column_name, "' must contain nucleotide sequences using A/C/G/T/N. ",
        "Invalid value(s): ",
        paste(head(bad_values, 10), collapse = ", ")
      )
    }
    
    x
  }
  
  validate_control_type <- function(x) {
    validate_no_na(x, "control/type")
    
    allowed_types <- c(
      "targeting",
      "non_targeting_control",
      "targeting_control"
    )
    
    x <- trimws(as.character(x))
    
    invalid <- !(x %in% allowed_types)
    
    if (any(invalid)) {
      bad_values <- unique(x[invalid])
      stop(
        "Column 'control/type' contains invalid entries. ",
        "Allowed values are: ",
        paste(allowed_types, collapse = ", "),
        ". Invalid value(s): ",
        paste(head(bad_values, 10), collapse = ", ")
      )
    }
    
    x
  }
  
  # ---------------------------------------------------------------------------
  # Read file
  # ---------------------------------------------------------------------------
  
  ext <- tolower(tools::file_ext(library_path))
  
  if (ext == "tsv") {
    df <- utils::read.delim(
      library_path,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else if (ext == "csv") {
    df <- utils::read.csv(
      library_path,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop(
        "Reading Excel files requires the 'readxl' package. ",
        "Install it with: install.packages('readxl')"
      )
    }
    
    df <- as.data.frame(
      readxl::read_excel(library_path),
      stringsAsFactors = FALSE
    )
  } else if (ext == "rds") {
    df <- readRDS(library_path)
    
    if (!is.data.frame(df)) {
      stop("The .rds library file must contain a dataframe.")
    }
    
    df <- as.data.frame(df, stringsAsFactors = FALSE)
  } else {
    stop(
      "Unsupported library file format: .", ext, "\n",
      "Supported formats are: .tsv, .csv, .xlsx, .xls, .rds"
    )
  }
  
  if (nrow(df) == 0) {
    stop("Library file contains zero rows.")
  }
  
  colnames(df) <- normalize_colnames(colnames(df))
  
  # ---------------------------------------------------------------------------
  # Find required and optional columns
  # ---------------------------------------------------------------------------
  
  col_sgrna_id <- find_column(
    df,
    canonical_name = "sgrna_id",
    alternatives = c("sgrna"),
    required = TRUE
  )
  
  col_seq <- find_column(
    df,
    canonical_name = "sgrna_seq",
    alternatives = c("seq", "sgrna_sequence", "sequence"),
    required = TRUE
  )
  
  col_entrez <- find_column(
    df,
    canonical_name = "entrez_id",
    alternatives = c("entrez"),
    required = TRUE
  )
  
  col_symbol <- find_column(
    df,
    canonical_name = "gene_symbol",
    alternatives = c("gene", "symbol"),
    required = TRUE
  )
  
  col_control <- find_column(
    df,
    canonical_name = "control",
    alternatives = c("type"),
    required = TRUE
  )
  
  col_count <- find_column(
    df,
    canonical_name = "count",
    alternatives = character(),
    required = FALSE
  )
  
  col_sublib <- find_column(
    df,
    canonical_name = "sublib",
    alternatives = c("sublibrary"),
    required = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # Validate and standardize columns
  # ---------------------------------------------------------------------------
  
  sgrna_id <- as.character(df[[col_sgrna_id]])
  validate_no_na(sgrna_id, "sgrna_id/sgrna")
  sgrna_id <- trimws(sgrna_id)
  sgrna_id <- make_unique_with_warning(sgrna_id)
  
  seq <- validate_nucleotide_sequence(df[[col_seq]], "sgrna_seq/seq/sgrna_sequence/sequence")
  
  entrez <- coerce_integerish(
    df[[col_entrez]],
    column_name = "entrez_id/entrez",
    allow_na = TRUE
  )
  
  symbol <- as.character(df[[col_symbol]])
  symbol[is_missing_string(symbol)] <- NA_character_
  symbol <- trimws(symbol)
  
  type <- validate_control_type(df[[col_control]])
  
  if (is.null(col_count)) {
    count <- rep(as.integer(default_count), nrow(df))
  } else {
    count <- coerce_integerish(
      df[[col_count]],
      column_name = "count",
      allow_na = FALSE
    )
  }
  
  if (is.null(col_sublib)) {
    sublib <- rep(default_sublib, nrow(df))
  } else {
    sublib <- as.character(df[[col_sublib]])
    sublib[is_missing_string(sublib)] <- NA_character_
    sublib <- trimws(sublib)
    
    # Fill missing sublib with default_sublib.
    sublib[is.na(sublib)] <- default_sublib
  }
  
  # ---------------------------------------------------------------------------
  # Build output dataframe
  # ---------------------------------------------------------------------------
  
  out <- data.frame(
    sgrna_id = sgrna_id,
    seq = seq,
    symbol = symbol,
    entrez = entrez,
    sublib = sublib,
    Gene = symbol,
    count = count,
    type = type,
    stringsAsFactors = FALSE
  )
  
  # Enforce exact column order
  out <- out[, c("sgrna_id", "seq", "symbol", "entrez", "sublib", "Gene", "count", "type")]
  
  out
}