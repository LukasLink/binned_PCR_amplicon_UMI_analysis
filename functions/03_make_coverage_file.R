#===============================================================================
# make coverage file
#===============================================================================
parse_coverage_file <- function(path) {
  lines <- readLines(path, warn = FALSE)
  
  mapped_output_folder <- get("mapped_output_folder", envir = .GlobalEnv)
  rds_output_folder <- get("rds_output_folder", envir = .GlobalEnv)
  
  # find each sample block start
  start_idx <- grep("^Checking wrong allignments for:", lines)
  if (length(start_idx) == 0) {
    stop("No sample blocks found. Expected lines starting with: 'Checking wrong allignments for:'")
  }
  
  # define block ends (right before next start, or end of file)
  end_idx <- c(start_idx[-1] - 1, length(lines))
  
  # small safe parsers
  parse_reads <- function(line) {
    if (is.na(line) || length(line) == 0) return(NA_real_)
    as.numeric(gsub(",", "", sub(".*?:\\s*([0-9,]+)\\s*\\(.*", "\\1", line)))
  }
  
  parse_perc <- function(line) {
    if (is.na(line) || length(line) == 0) return(NA_real_)
    as.numeric(sub(".*\\(([-0-9.]+)%\\).*", "\\1", line))
  }
  
  parse_coverage <- function(line) {
    if (is.na(line) || length(line) == 0) return(NA_real_)
    as.numeric(sub("^sgRNA coverage:\\s*([-0-9.]+)%.*", "\\1", line))
  }
  
  out <- lapply(seq_along(start_idx), function(k) {
    block <- lines[start_idx[k]:end_idx[k]]
    
    sample_name <- sub("^Checking wrong allignments for:\\s*", "", block[1])
    
    correct_line <- block[grep("^\\s*Correct", block)][1]
    wrong_line   <- block[grep("^\\s*Wrong", block)][1]
    cov_line     <- block[grep("^sgRNA coverage:", block)][1]
    
    correct_reads <- parse_reads(correct_line)
    correct_perc  <- parse_perc(correct_line)
    
    wrong_reads <- parse_reads(wrong_line)
    wrong_perc  <- parse_perc(wrong_line)
    
    coverage <- parse_coverage(cov_line)
    
    data.frame(
      sample_name   = sample_name,
      correct_reads = correct_reads,
      correct_perc  = paste0(sprintf("%.2f", correct_perc), "%"),
      wrong_reads   = wrong_reads,
      wrong_perc    = paste0(sprintf("%.2f", wrong_perc), "%"),
      coverage      = paste0(sprintf("%.2f", coverage), "%"),
      stringsAsFactors = FALSE
    ) %>% mutate(
      sample_name = sub(" ", "", sample_name)
    )
  })
  
  do.call(rbind, out)
}
add_star_log_stats <- function(df, mapped_output_folder) {
  
  mapped_output_folder <- get("mapped_output_folder", envir = .GlobalEnv)
  rds_output_folder <- get("rds_output_folder", envir = .GlobalEnv)  
  # Map STAR log "labels" -> dataframe column names
  metrics_map <- c(
    "Number of input reads"                         = "input_reads",
    "Average input read length"                     = "avg_input_read_length",
    "Uniquely mapped reads number"                  = "uniquely_mapped_reads",
    "Uniquely mapped reads %"                       = "uniquely_mapped_perc",
    "Average mapped length"                         = "avg_mapped_length",
    "Mismatch rate per base, %"                     = "mismatch_rate_per_base_perc",
    "Deletion rate per base"                        = "deletion_rate_per_base_perc",
    "Deletion average length"                       = "deletion_avg_length",
    "Insertion rate per base, %"                    = "insertion_rate_per_base_perc",
    "Insertion average length"                      = "insertion_avg_length",
    "Number of reads mapped to multiple loci"       = "multi_loci_reads",
    "% of reads mapped to multiple loci"            = "multi_loci_reads_perc",
    "Number of reads mapped to too many loci"       = "too_many_loci_reads",
    "% of reads mapped to too many loci"            = "too_many_loci_reads_perc",
    "Number of reads unmapped: too many mismatches" = "unmapped_too_many_mismatches",
    "% of reads unmapped: too many mismatches"      = "unmapped_too_many_mismatches_perc",
    "Number of reads unmapped: too short"           = "unmapped_too_short",
    "% of reads unmapped: too short"                = "unmapped_too_short_perc",
    "Number of reads unmapped: other"               = "unmapped_other",
    "% of reads unmapped: other"                    = "unmapped_other_perc"
  )
  
  # Helper: parse a single STAR Log.final.out file -> one-row data.frame
  parse_star_log_file <- function(sample_name) {
    log_path <- get_file_path(mapped_output_folder,
                              paste0(sample_name, "_Log.final.out"))
    
    # If file missing, return all NAs
    if (!file.exists(log_path)) {
      warning("STAR log not found for sample '", sample_name, "': ", log_path)
      na_vals <- setNames(rep(NA_real_, length(metrics_map)), metrics_map)
      return(as.data.frame(as.list(na_vals)))
    }
    
    lines <- readLines(log_path, warn = FALSE)
    
    # Extract numeric value from the line after the "|" and remove "%"
    extract_value <- function(key) {
      idx <- which(grepl(key, lines, fixed = TRUE))[1]
      if (length(idx) == 0 || is.na(idx)) return(NA_real_)
      
      parts <- strsplit(lines[idx], "\\|")[[1]]
      if (length(parts) < 2) return(NA_real_)
      
      val_str <- trimws(parts[2])
      as.numeric(gsub("%", "", val_str))
    }
    
    vals <- vapply(names(metrics_map), extract_value, numeric(1))
    names(vals) <- base::unname(metrics_map)
    
    as.data.frame(as.list(vals), stringsAsFactors = FALSE)
  }
  
  # Parse all samples and bind into one df
  stats_df <- do.call(
    rbind,
    lapply(df$sample_name, parse_star_log_file)
  )
  
  # Append columns to the original df
  cbind(df, stats_df)
}
