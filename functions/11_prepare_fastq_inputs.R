prepare_fastq_inputs <- function(
    fastq_dir,
    fastq_name_table_file_path = NULL,
    output_symlink_dir,
    manifest_output_path = NULL,
    strict_file_match = TRUE,
    overwrite_symlinks = TRUE
) {
  # -----------------------------
  # Helper functions
  # -----------------------------
  
  stop_with_message <- function(...) {
    stop(paste0(...), call. = FALSE)
  }
  
  normalize_optional_path <- function(x) {
    if (is.null(x)) return(NULL)
    if (length(x) == 0) return(NULL)
    if (is.na(x)) return(NULL)
    
    x <- trimws(as.character(x))
    
    if (x %in% c("", "NULL", "null", "NA", "None", "none")) {
      return(NULL)
    }
    
    x
  }
  
  get_fastq_extension <- function(x) {
    # Supported:
    # .txt, .txt.gz, .txt.gzip
    # .fq, .fq.gz, .fq.gzip
    # .fastq, .fastq.gz, .fastq.gzip
    
    m <- regexpr(
      "\\.(txt|fq|fastq)(\\.gz|\\.gzip)?$",
      basename(x),
      ignore.case = TRUE
    )
    
    ifelse(m == -1, NA_character_, regmatches(basename(x), m))
  }
  
  strip_fastq_extension <- function(x) {
    sub(
      "\\.(txt|fq|fastq)(\\.gz|\\.gzip)?$",
      "",
      basename(x),
      ignore.case = TRUE
    )
  }
  
  validate_bin <- function(x) {
    !is.na(x) & x %in% c("I", "L", "U")
  }
  
  validate_sublibrary <- function(x) {
    !is.na(x) & grepl("^L[0-9]+$", x)
  }
  
  validate_sample <- function(x) {
    !is.na(x) & nzchar(x) & !grepl("_", x)
  }
  
  validate_pipeline_name <- function(x) {
    grepl("^(I|L|U)_(L[0-9]+)_[^_]+$", x)
  }
  
  parse_pipeline_name <- function(pipeline_name) {
    parts <- strsplit(pipeline_name, "_", fixed = TRUE)
    
    data.frame(
      bin = vapply(parts, `[`, character(1), 1),
      sublibrary = vapply(parts, `[`, character(1), 2),
      sample = vapply(parts, `[`, character(1), 3),
      stringsAsFactors = FALSE
    )
  }
  
  make_error_list <- function(x, bullet = "  - ") {
    paste0(bullet, x, collapse = "\n")
  }
  
  # -----------------------------
  # Normalize inputs
  # -----------------------------
  
  fastq_name_table_file_path <- normalize_optional_path(
    fastq_name_table_file_path
  )
  
  if (missing(fastq_dir) || is.null(fastq_dir) || !nzchar(fastq_dir)) {
    stop_with_message("`fastq_dir` must be provided.")
  }
  
  if (!dir.exists(fastq_dir)) {
    stop_with_message(
      "Input FASTQ directory does not exist:\n  ",
      fastq_dir
    )
  }
  
  if (missing(output_symlink_dir) || is.null(output_symlink_dir) || !nzchar(output_symlink_dir)) {
    stop_with_message("`output_symlink_dir` must be provided.")
  }
  
  # -----------------------------
  # Detect input files
  # -----------------------------
  
  all_files <- list.files(
    fastq_dir,
    full.names = TRUE,
    recursive = FALSE,
    include.dirs = FALSE
  )
  
  all_files <- all_files[file.info(all_files)$isdir == FALSE]
  
  supported_ext <- get_fastq_extension(all_files)
  input_files <- all_files[!is.na(supported_ext)]
  
  if (length(input_files) == 0) {
    stop_with_message(
      "No supported input files were found in:\n  ",
      fastq_dir,
      "\n\nSupported extensions are:\n",
      "  - .txt\n",
      "  - .txt.gz\n",
      "  - .txt.gzip\n",
      "  - .fq\n",
      "  - .fq.gz\n",
      "  - .fq.gzip\n",
      "  - .fastq\n",
      "  - .fastq.gz\n",
      "  - .fastq.gzip"
    )
  }
  
  input_file_basename <- basename(input_files)
  input_file_ext <- get_fastq_extension(input_files)
  input_file_stem <- strip_fastq_extension(input_files)
  
  if (anyDuplicated(input_file_basename)) {
    duplicated_files <- unique(input_file_basename[duplicated(input_file_basename)])
    stop_with_message(
      "Duplicate input file names were found. File basenames must be unique:\n",
      make_error_list(duplicated_files)
    )
  }
  
  # -----------------------------
  # Case 1: xlsx parser provided
  # -----------------------------
  
  if (!is.null(fastq_name_table_file_path)) {
    if (!file.exists(fastq_name_table_file_path)) {
      stop_with_message(
        "The provided input FASTQ name parser file does not exist:\n  ",
        fastq_name_table_file_path
      )
    }
    
    if (!grepl("\\.xlsx$", fastq_name_table_file_path, ignore.case = TRUE)) {
      stop_with_message(
        "The input FASTQ name parser must be an .xlsx file:\n  ",
        fastq_name_table_file_path
      )
    }
    
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop_with_message(
        "Package `readxl` is required to read the input FASTQ name parser .xlsx file.\n",
        "Install it with:\n  install.packages(\"readxl\")"
      )
    }
    
    available_sheets <- readxl::excel_sheets(fastq_name_table_file_path)
    
    if (!"fastq_names" %in% available_sheets) {
      stop_with_message(
        "The input FASTQ name parser xlsx must contain a sheet named exactly `fastq_names`.\n\n",
        "Available sheets are:\n",
        make_error_list(available_sheets)
      )
    }
    
    fastq_names_sheet_index <- which(available_sheets == "fastq_names")
    
    if (fastq_names_sheet_index != 2) {
      stop_with_message(
        "The sheet `fastq_names` must be the second sheet in the xlsx file.\n\n",
        "Current sheet order is:\n",
        paste0(seq_along(available_sheets), ". ", available_sheets, collapse = "\n")
      )
    }
    
    parser <- readxl::read_excel(
      fastq_name_table_file_path,
      sheet = "fastq_names"
    )
    
    parser <- as.data.frame(parser, stringsAsFactors = FALSE)
    
    required_cols <- c("original_file", "bin", "sublibrary", "sample")
    missing_cols <- setdiff(required_cols, colnames(parser))
    
    if (length(missing_cols) > 0) {
      stop_with_message(
        "The parser xlsx is missing required columns:\n",
        make_error_list(missing_cols),
        "\n\nRequired columns are:\n",
        make_error_list(required_cols)
      )
    }
    
    parser <- parser[, required_cols, drop = FALSE]
    
    parser$original_file <- trimws(as.character(parser$original_file))
    parser$bin <- trimws(as.character(parser$bin))
    parser$sublibrary <- trimws(as.character(parser$sublibrary))
    parser$sample <- trimws(as.character(parser$sample))
    
    if (any(is.na(parser$original_file) | parser$original_file == "")) {
      stop_with_message("The parser xlsx contains empty values in `original_file`.")
    }
    
    if (anyDuplicated(parser$original_file)) {
      duplicated_original <- unique(parser$original_file[duplicated(parser$original_file)])
      stop_with_message(
        "The parser xlsx contains duplicate `original_file` entries:\n",
        make_error_list(duplicated_original)
      )
    }
    
    invalid_bin <- parser$original_file[!validate_bin(parser$bin)]
    invalid_sublibrary <- parser$original_file[!validate_sublibrary(parser$sublibrary)]
    invalid_sample <- parser$original_file[!validate_sample(parser$sample)]
    
    if (length(invalid_bin) > 0) {
      stop_with_message(
        "Invalid `bin` values in parser xlsx.\n\n",
        "Allowed bin values are: I, L, U\n\n",
        "Affected files:\n",
        make_error_list(invalid_bin)
      )
    }
    
    if (length(invalid_sublibrary) > 0) {
      stop_with_message(
        "Invalid `sublibrary` values in parser xlsx.\n\n",
        "Allowed sublibrary pattern is: L followed by digits, for example L1, L01, L002.\n\n",
        "Affected files:\n",
        make_error_list(invalid_sublibrary)
      )
    }
    
    if (length(invalid_sample) > 0) {
      stop_with_message(
        "Invalid `sample` values in parser xlsx.\n\n",
        "Sample names must be non-empty and must not contain underscores.\n\n",
        "Affected files:\n",
        make_error_list(invalid_sample)
      )
    }
    
    parser$pipeline_name <- paste(
      parser$bin,
      parser$sublibrary,
      parser$sample,
      sep = "_"
    )
    
    if (anyDuplicated(parser$pipeline_name)) {
      duplicated_pipeline_names <- unique(parser$pipeline_name[duplicated(parser$pipeline_name)])
      stop_with_message(
        "The parser xlsx generates duplicate pipeline names:\n",
        make_error_list(duplicated_pipeline_names),
        "\n\nEach combination of bin, sublibrary, and sample must be unique."
      )
    }
    
    # Match xlsx entries to files in the input directory.
    # We require original_file to match the basename of the files.
    files_in_dir <- input_file_basename
    files_in_xlsx <- parser$original_file
    
    missing_from_dir <- setdiff(files_in_xlsx, files_in_dir)
    missing_from_xlsx <- setdiff(files_in_dir, files_in_xlsx)
    
    if (length(missing_from_dir) > 0) {
      stop_with_message(
        "Some files listed in the parser xlsx were not found in the input directory:\n",
        make_error_list(missing_from_dir)
      )
    }
    
    if (strict_file_match && length(missing_from_xlsx) > 0) {
      stop_with_message(
        "Some supported input files in the directory are missing from the parser xlsx:\n",
        make_error_list(missing_from_xlsx),
        "\n\nBecause `strict_file_match = TRUE`, every supported input file must be listed in the parser xlsx."
      )
    }
    
    matched_idx <- match(parser$original_file, input_file_basename)
    
    manifest <- data.frame(
      original_file = input_files[matched_idx],
      original_file_basename = input_file_basename[matched_idx],
      original_extension = input_file_ext[matched_idx],
      bin = parser$bin,
      sublibrary = parser$sublibrary,
      sample = parser$sample,
      pipeline_name = parser$pipeline_name,
      stringsAsFactors = FALSE
    )
    
  } else {
    # -----------------------------
    # Case 2: no xlsx parser provided
    # Validate existing filenames directly
    # -----------------------------
    
    valid_name <- validate_pipeline_name(input_file_stem)
    
    if (!all(valid_name)) {
      invalid_files <- input_file_basename[!valid_name]
      
      stop_with_message(
        "Some input files do not match the required naming pattern.\n\n",
        "Required pattern:\n",
        "  BIN_SUBLIBRARY_SAMPLE.fastq.gz\n\n",
        "Where:\n",
        "  BIN        = I, L, or U\n",
        "  SUBLIBRARY = L followed by digits, e.g. L1, L01, L002\n",
        "  SAMPLE     = any string without underscores\n\n",
        "Examples:\n",
        "  I_L1_control.fastq.gz\n",
        "  L_L1_treated.fastq.gz\n",
        "  U_L1_treated.fastq.gz\n\n",
        "Invalid files:\n",
        make_error_list(invalid_files)
      )
    }
    
    parsed <- parse_pipeline_name(input_file_stem)
    
    manifest <- data.frame(
      original_file = input_files,
      original_file_basename = input_file_basename,
      original_extension = input_file_ext,
      bin = parsed$bin,
      sublibrary = parsed$sublibrary,
      sample = parsed$sample,
      pipeline_name = input_file_stem,
      stringsAsFactors = FALSE
    )
  }
  
  # -----------------------------
  # Create symlink names
  # -----------------------------
  
  manifest$symlink_file_basename <- paste0(
    manifest$pipeline_name,
    manifest$original_extension
  )
  
  if (anyDuplicated(manifest$symlink_file_basename)) {
    duplicated_symlinks <- unique(
      manifest$symlink_file_basename[duplicated(manifest$symlink_file_basename)]
    )
    
    stop_with_message(
      "Duplicate symlink names would be created:\n",
      make_error_list(duplicated_symlinks)
    )
  }
  
  manifest$symlink_file <- file.path(
    output_symlink_dir,
    manifest$symlink_file_basename
  )
  
  # -----------------------------
  # Prepare symlink directory
  # -----------------------------
  
  if (!dir.exists(output_symlink_dir)) {
    dir.create(output_symlink_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  if (!dir.exists(output_symlink_dir)) {
    stop_with_message(
      "Could not create output symlink directory:\n  ",
      output_symlink_dir
    )
  }
  
  existing_symlinks <- manifest$symlink_file[file.exists(manifest$symlink_file)]
  
  if (length(existing_symlinks) > 0) {
    if (!overwrite_symlinks) {
      stop_with_message(
        "Some target symlink files already exist:\n",
        make_error_list(existing_symlinks),
        "\n\nSet `overwrite_symlinks = TRUE` to replace them."
      )
    } else {
      unlink(existing_symlinks)
    }
  }
  
  # -----------------------------
  # Create symlinks
  # -----------------------------
  
  symlink_ok <- file.symlink(
    from = normalizePath(manifest$original_file, mustWork = TRUE),
    to = manifest$symlink_file
  )
  
  if (!all(symlink_ok)) {
    failed <- manifest$symlink_file[!symlink_ok]
    
    stop_with_message(
      "Failed to create one or more symlinks:\n",
      make_error_list(failed),
      "\n\nThis can happen on some systems if symlink creation is not permitted."
    )
  }
  
  # -----------------------------
  # Write manifest
  # -----------------------------
  
  if (!is.null(manifest_output_path)) {
    manifest_dir <- dirname(manifest_output_path)
    
    if (!dir.exists(manifest_dir)) {
      dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    write.table(
      manifest,
      file = manifest_output_path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  
  # -----------------------------
  # Return manifest invisibly/visibly
  # -----------------------------
  
  message(
    "Prepared ",
    nrow(manifest),
    " input files.\n",
    "Standardized symlinks written to:\n  ",
    output_symlink_dir
  )
  
  return(manifest)
}