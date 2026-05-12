#!/usr/bin/env Rscript

# ==============================================================================
# Generate bcwithqc config.json
# ==============================================================================

suppressPackageStartupMessages({
  library(jsonlite)
})


# ------------------------------------------------------------------------------
# Main function
# ------------------------------------------------------------------------------

generate_bcwithqc_config <- function(
    merged_sgRNA_df,
    output_json_path,
    unknown_read_orientation = FALSE,
    single_end_reads = TRUE,
    stagger_list = c(
      "T",
      "AT",
      "GAT",
      "CGAT",
      "TCGAT",
      "ATCGAT"
    ),
    first_constant_sequence = "CTTGTGGAAAGGACGAAACACCG",
    stagger_min_length = 0,
    stagger_max_length = 4,
    second_constant_sequence = "GTTTAAGAGCTATGCTGGAAA",
    stagger_maxerrors = 0,
    first_constant_maxerrors = 1,
    sgrna_maxerrors = 1,
    second_constant_maxerrors = 1,
    keep_nonbarcode = TRUE
) {
  # ---------------------------------------------------------------------------
  # Basic checks
  # ---------------------------------------------------------------------------
  
  required_cols <- c("sgrna_id", "seq")
  
  missing_cols <- base::setdiff(required_cols, colnames(merged_sgRNA_df))
  
  if (length(missing_cols) > 0) {
    stop(
      "merged_sgRNA_df is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  if (any(is.na(merged_sgRNA_df$sgrna_id))) {
    stop("merged_sgRNA_df$sgrna_id contains NA values.")
  }
  
  if (any(is.na(merged_sgRNA_df$seq))) {
    stop("merged_sgRNA_df$seq contains NA values.")
  }
  
  if (anyDuplicated(merged_sgRNA_df$sgrna_id)) {
    warning(
      "merged_sgRNA_df$sgrna_id contains duplicated names. ",
      "bcwithqc may require unique barcode names.",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Construct barcode list from merged_sgRNA_df
  # ---------------------------------------------------------------------------
  
  sgrna_names <- as.character(merged_sgRNA_df$sgrna_id)
  sgrna_sequences <- as.character(merged_sgRNA_df$seq)
  
  # ---------------------------------------------------------------------------
  # Build config object
  # ---------------------------------------------------------------------------
  
  config <- list(
    unknown_read_orientation = unknown_read_orientation,
    single_end_reads = single_end_reads,
    
    barcode_struct_r1 = list(
      keep_nonbarcode = keep_nonbarcode,
      
      blocks = list(
        list(
          maxerrors = stagger_maxerrors,
          blockfunction = "discard",
          sequence = as.list(stagger_list),
          blocktype = "constantRegion",
          blockname = "stagger",
          keepblock = FALSE
        ),
        
        list(
          maxerrors = first_constant_maxerrors,
          blockfunction = "barcode",
          sequence = as.list(first_constant_sequence),
          blocktype = "constantRegion",
          blockname = "con1",
          keepblock = TRUE
        ),
        
        list(
          maxerrors = sgrna_maxerrors,
          blockfunction = "barcode",
          name = as.list(sgrna_names),
          sequence = as.list(sgrna_sequences),
          blocktype = "barcodeList",
          blockname = "sgrna",
          keepblock = TRUE
        ),
        
        list(
          maxerrors = second_constant_maxerrors,
          blockfunction = "barcode",
          sequence = as.list(second_constant_sequence),
          blocktype = "constantRegion",
          blockname = "con2",
          keepblock = TRUE
        )
      )
    )
  )
  
  # ---------------------------------------------------------------------------
  # Write JSON
  # ---------------------------------------------------------------------------
  
  output_dir <- dirname(output_json_path)
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  jsonlite::write_json(
    config,
    path = output_json_path,
    pretty = TRUE,
    auto_unbox = TRUE,
    null = "null"
  )
  
  message("Wrote bcwithqc config to: ", normalizePath(output_json_path))
  
  invisible(config)
}

# ==============================================================================
# Example usage
# ==============================================================================
# If running interactively or sourcing from an Rmd:
#
config <- generate_bcwithqc_config(
  merged_sgRNA_df = merged_sgRNA_df,
  output_json_path = file.path("/g/steinmetz/link/Amplicon_barcode_analysis/bcwithqc_test", "bcwithqc_config_stagger.json")
)
# 
# If loading merged_sgRNA_df from RDS:
# 
# merged_sgRNA_df <- readRDS("path/to/merged_sgRNA_df.rds")
# 
# generate_bcwithqc_config(
#   merged_sgRNA_df = merged_sgRNA_df,
#   output_json_path = "config.json"
# )