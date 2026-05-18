#===============================================================================
# get_count_df_long
#===============================================================================
read_file_to_df <- function(file_name,
                            folder_path = dedup_output_folder,
                            suffix_to_rm = "_dedup_idxstats.txt",
                            threshold_df = NULL,
                            check_alignments = TRUE
) {
  
  file_path <- file.path(folder_path, file_name)
  df <- read.table(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  name <- sub(suffix_to_rm, "", file_name)
  sub_lib <- stringr::str_match(name, "^[A-Za-z]+_(L\\d+)_[^_]+")[, 2]
  exp <- stringr::str_replace(name, "^([^_]+_)", "")
  
  if (is.null(threshold_df)) {
    threshold <- 0
  } else {
    threshold <- threshold_df$threshold[threshold_df$replicate == name]
  }
  
  df <- df %>% 
    dplyr::rename(sgRNA = V1, count = V3, length = V2) %>%
    dplyr::select(-V4) %>%
    dplyr::mutate(
      count = dplyr::if_else(count <= threshold, 0, count),
      exp = exp
    ) %>%
    dplyr::filter(sgRNA != "*")
  
  # ---------------------------------------------------------------------------
  # Prepare library annotation table
  # ---------------------------------------------------------------------------
  # New expected behavior:
  # merged_sgRNA_df should contain a `type` column with values:
  #   - non_targeting_control
  #   - targeting_control
  #   - targeting
  #
  # Backward-compatible fallback:
  # If `type` is missing, infer controls from sgRNA names using the old grep logic.
  # ---------------------------------------------------------------------------
  
  if ("type" %in% colnames(merged_sgRNA_df)) {
    
    check_df <- merged_sgRNA_df %>%
      dplyr::select(sgrna_id, sublib, type)
    
  } else {
    
    logger::log_warn(
      paste0(
        "`merged_sgRNA_df` does not contain a `type` column. ",
        "Falling back to grep-based control annotation using sgRNA names. ",
        "Please update the library file to include a `type` column with values: ",
        "non_targeting_control, targeting_control, targeting."
      )
    )
    
    check_df <- merged_sgRNA_df %>%
      dplyr::select(sgrna_id, sublib) %>%
      dplyr::mutate(
        type = dplyr::case_when(
          grepl("^CONTROL_C_NONTARG_", sgrna_id) ~ "non_targeting_control",
          grepl("^CONTROL_C_", sgrna_id) ~ "targeting_control",
          TRUE ~ "targeting"
        )
      )
  }
  
  # Optional but recommended: validate allowed type values
  allowed_types <- c(
    "non_targeting_control",
    "targeting_control",
    "targeting"
  )
  
  invalid_types <- unique(check_df$type[!check_df$type %in% allowed_types])
  
  if (length(invalid_types) > 0) {
    stop(
      "`merged_sgRNA_df$type` contains invalid values:\n",
      paste0("  - ", invalid_types, collapse = "\n"),
      "\n\nAllowed values are:\n",
      paste0("  - ", allowed_types, collapse = "\n"),
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Join counts to library annotation
  # ---------------------------------------------------------------------------
  
  joined_df <- df %>%
    dplyr::left_join(check_df, by = c("sgRNA" = "sgrna_id"))
  
  # Check for non-unique sgRNA in the result
  if (any(duplicated(joined_df$sgRNA))) {
    
    logger::log_warn(
      "Non-unique `sgRNA` after join; redoing join with 1:1 pairing by order."
    )
    
    df_idx <- df %>%
      dplyr::group_by(sgRNA) %>%
      dplyr::mutate(.pair_id = dplyr::row_number()) %>%
      dplyr::ungroup()
    
    check_idx <- check_df %>%
      dplyr::group_by(sgrna_id) %>%
      dplyr::mutate(.pair_id = dplyr::row_number()) %>%
      dplyr::ungroup()
    
    joined_df <- df_idx %>%
      dplyr::left_join(
        check_idx,
        by = c("sgRNA" = "sgrna_id", ".pair_id")
      ) %>%
      dplyr::select(-.pair_id)
    
    df <- joined_df
    
  } else {
    
    df <- joined_df
  }
  
  # ---------------------------------------------------------------------------
  # Alignment sanity checks
  # ---------------------------------------------------------------------------
  
  if (check_alignments == TRUE) {
    
    logger::log_info("Checking wrong alignments for: {name}")
    
    # Filter non-control rows with count > 0
    non_control_df <- df %>%
      dplyr::filter(type == "targeting", count > 0)
    
    correct_sum <- non_control_df %>%
      dplyr::filter(sublib == sub_lib) %>%
      dplyr::summarise(total = sum(count), .groups = "drop") %>%
      dplyr::pull(total)
    
    wrong_sum <- non_control_df %>%
      dplyr::filter(sublib != sub_lib) %>%
      dplyr::summarise(total = sum(count), .groups = "drop") %>%
      dplyr::pull(total)
    
    total_sum <- correct_sum + wrong_sum
    
    if (total_sum == 0) {
      logger::log_warn(
        "In file: {name} no reads/UMIs were detected. This is bad >.<"
      )
    } else if (total_sum < 1000) {
      logger::log_warn(
        paste0(
          "In file: {name} less than 1000 reads/UMIs were detected. ",
          "If you were not sequencing with low depth this is bad!"
        )
      )
    }
    
    logger::log_info("UMI/read count alignment stats:")
    
    if (total_sum > 0) {
      logger::log_info(
        "  Correct aligned to sgRNAs in sublibrary:     {correct_sum} ({round(100 * correct_sum / total_sum, 2)}%)"
      )
      logger::log_info(
        "  Wrong aligned to sgRNAs not in sublibrary:   {wrong_sum} ({round(100 * wrong_sum / total_sum, 2)}%)"
      )
      
      if (!is.na(wrong_sum / total_sum) && (wrong_sum / total_sum) > 0.5) {
        logger::log_warn(
          paste0(
            "More than 50% of reads aligned to sgRNAs from the wrong sublibrary. ",
            "Check whether the sublibraries in the library file and ",
            "fastq_name_table_xlsx have the same numbering."
          )
        )
      }
      
    } else {
      logger::log_info("  Correct aligned to sgRNAs in sublibrary:     0")
      logger::log_info("  Wrong aligned to sgRNAs not in sublibrary:   0")
    }
  }
  
  # ---------------------------------------------------------------------------
  # Filter to current sublibrary
  # ---------------------------------------------------------------------------
  
  same_controls_in_all_sublibraries <- get(
    "same_controls_in_all_sublibraries",
    envir = .GlobalEnv
  )
  
  if (same_controls_in_all_sublibraries == TRUE) {
    
    df <- df %>%
      dplyr::filter(
        sublib == sub_lib |
          type %in% c("non_targeting_control", "targeting_control")
      ) %>%
      dplyr::select(-sublib)
    
  } else {
    
    df <- df %>%
      dplyr::filter(sublib == sub_lib) %>%
      dplyr::select(-sublib)
  }
  n_non_targeting_controls <- df %>%
    dplyr::filter(type == "non_targeting_control") %>%
    dplyr::pull(sgRNA) %>%
    unique() %>%
    length()
  
  if (n_non_targeting_controls < 3) {
    logger::log_warn(
      paste0(
        "Only ", n_non_targeting_controls,
        " non-targeting controls detected after filtering. ",
        "MAUDE likely requires at least 3 non-targeting controls and may break downstream."
      )
    )
  } else if (n_non_targeting_controls < 50) {
    logger::log_warn(
      paste0(
        "Only ", n_non_targeting_controls,
        " non-targeting controls detected after filtering. ",
        "This may be statistically risky for normalization/modeling."
      )
    )
  }
  if (check_alignments == TRUE) {
    logger::log_info(
      "sgRNA coverage: {round(sum(df$count > 0) / nrow(df) * 100, 2)}%"
    )
    logger::log_info("--------------------------------------------")
  }
  
  # ---------------------------------------------------------------------------
  # group_category now comes from library type
  # ---------------------------------------------------------------------------
  
  df <- df %>%
    dplyr::mutate(group_category = type) %>%
    dplyr::select(-type)
  
  return(df)
}

process_folder_files <- function(folder_path,
                                 threshold_df = NULL,
                                 skip_list = c(),
                                 check_alignments = TRUE){
  
  data_type <- get("data_type", envir = .GlobalEnv)
  
  # List all matching files in the folder
  if (data_type == "umis") {
    files <- list.files(
      folder_path,
      pattern = "^[ILU]_L\\d+_[^_]+_dedup_idxstats\\.txt$",
      full.names = TRUE
    )
    
  } else if (data_type == "reads") {
    files <- list.files(
      folder_path,
      pattern = "^[ILU]_L\\d+_[^_]+_Aligned\\.sortedByCoord\\.out_idxstats\\.txt$",
      full.names = TRUE
    )
    
  } else {
    stop("Invalid data_type. data_type must be 'reads' or 'umis'")
  }
  
  
  
  # Read each file and combine them
  combined_df <- bind_rows(lapply(files, function(file) {
    # Extract the filename
    file_name <- basename(file)
    if (data_type == "umis"){
      name <- sub("_dedup_idxstats.txt","",file_name)
      # Extract condition (X), sublib (L#), and sample (#) using regex
      matches <- str_match(file_name, "^([ILU])_L(\\d+)_([^_]+)_dedup_idxstats\\.txt$")
      suffix_to_remove <- "_dedup_idxstats.txt"
    }
    if (data_type == "reads"){
      name <- sub("_Aligned.sortedByCoord.out_idxstats.txt","",file_name)
      # Extract condition (X), sublib (L#), and sample (#) using regex
      matches <- str_match(file_name, "^([ILU])_L(\\d+)_([^_]+)_Aligned\\.sortedByCoord\\.out_idxstats\\.txt$")
      suffix_to_remove <- "_Aligned.sortedByCoord.out_idxstats.txt"
    }
    
    # Skip entries from skip list
    skip_this_file <- any(vapply(skip_list, function(x) {
      if (startsWith(x, "re:")) {
        grepl(sub("^re:", "", x), name, perl = TRUE)
      } else {
        identical(x, name)
      }
    }, logical(1)))
    
    if (skip_this_file){
      print(paste(name, "skipped due to skip list"))
      return(NULL)
    }
    
    if (is.na(matches[1])) {
      stop(paste("Filename does not match expected pattern:", file_name))
    }
    
    condition <- case_when(
      matches[2] == "I" ~ "input",
      matches[2] == "L" ~ "lower",
      matches[2] == "U" ~ "upper"
    )
    
    sublib <- paste0("sublib_", matches[3])
    sample <- paste0("sample_", matches[4])
    
    # Read the file using your existing function
    df <- read_file_to_df(file_name = file_name,
                          folder_path = folder_path,
                          threshold_df = threshold_df,
                          check_alignments = check_alignments,
                          suffix_to_rm = suffix_to_remove)
    
    # Add the new columns
    df <- df %>%
      mutate(condition = condition, sublib = sublib, sample = sample) %>% 
      select(-length)
    
    return(df)
  }))
  
  return(combined_df)
}
read_bcwithqc_data <- function(dir_path, sgRNA_df) {
  # Step 1: Read barcodes.tsv.gz and get seq values
  
  barcodes_path <- file.path(dir_path, "barcodes.tsv.gz")
  barcode_seqs <- read_tsv(barcodes_path, col_names = FALSE, show_col_types = FALSE) %>%
    pull(1) %>%
    basename()
  
  # Step 2: Read matrix.mtx.gz (only diagonal values matter)
  matrix_path <- file.path(dir_path, "matrix.mtx.gz")
  mat <- readMM(matrix_path)
  
  if (!is(mat, "dgTMatrix")) {
    mat <- as(mat, "dgTMatrix")
  }
  
  # Extract diagonal (should match barcodes length)
  diag_values <- Matrix::diag(mat)
  
  # Step 3: Construct dataframe with seq and count
  count_df <- tibble(
    seq = barcode_seqs,
    count = diag_values
  )
  
  # Step 4: Join with sgRNA_df by seq
  joined_df <- sgRNA_df %>%
    select(seq, sgrna_id) %>% 
    inner_join(count_df, by = "seq") %>% 
    select(sgrna_id,count) %>% 
    mutate(count = as.integer(count))
  
  return(joined_df)
}

process_bcwithqc_data <- function(parent_dir,
                              merged_sgRNA_df = merged_sgRNA_df,
                              data_type = "reads",
                              skip_list = c(),
                              check_alignments = TRUE) {
  # Get list of matching subdirectories
  subdirs <- list.dirs(parent_dir, recursive = FALSE, full.names = TRUE)
  subdirs <- subdirs[grepl("_output_all$", basename(subdirs))]
  
  merged_sgRNA_df$sublib <- gsub("ICS", "sublib", merged_sgRNA_df$sublib)
  
  # Initialize list to collect data
  results <- list()
  
  for (dir_path in subdirs) {
    folder_name <- basename(dir_path)
    
    match <- str_match(folder_name, "^([A-Z])([A-Z])([0-9])([0-9])_output_all$")
    
    if (any(is.na(match))) {
      warning(paste("Skipping folder with unexpected name:", folder_name))
      next
    }
    name <- paste0(match[2],"_",match[3],match[4],"_",match[5])
    
    # Skip if in skip list.
    is_skipped <- any(vapply(skip_list, function(x) {
      if (startsWith(x, "re:")) {
        grepl(sub("^re:", "", x), name, perl = TRUE)
      } else {
        identical(x, name)
      }
    }, logical(1)))
    if (skip_this_file){next}
    
    condition <- case_when(
      match[2] == "I" ~ "input",
      match[2] == "L" ~ "lower",
      match[2] == "U" ~ "upper",
      TRUE ~ NA_character_
    )
    
    sublib <- paste0("sublib_", match[4])
    sample <- paste0("sample_", match[5])
    exp <- paste0(match[3], match[4], "_", match[5])
    
    
    # Path to raw_umis_bc_matrix
    if (data_type == "umis"){
      matrix_dir <- file.path(dir_path, "raw_umis_bc_matrix")
    }
    if (data_type == "reads"){
      matrix_dir <- file.path(dir_path, "raw_reads_bc_matrix")
    }
    
    # Safely read data
    df <- tryCatch({
      df <- read_bcwithqc_data(matrix_dir, merged_sgRNA_df)
      df <- complete_sgrna_df(df,
                              merged_sgRNA_df,
                              sublib,
                              name = name,
                              check_alignments = check_alignments)
      
      df %>%
        mutate(
          condition = condition,
          sublib = sublib,
          sample = sample,
          exp = exp
        )
    }, error = function(e) {
      cat("Skipping due to error in", folder_name, "->", e$message, "\n")
      return(NULL)
    })
    
    if (!is.null(df)) {
      # Rename and add group_category
      df <- df %>%
        rename(sgRNA = sgrna_id) %>%
        mutate(group_category = case_when(
          grepl("^CONTROL_C_NONTARG_", sgRNA) ~ "non_targeting_control",
          grepl("^CONTROL_C_", sgRNA) ~ "targeting_control",
          TRUE ~ "targeting"
        ))
      
      results[[length(results) + 1]] <- df
    } else {
      print(paste("WARNING:",matrix_dir,"yielded df -> NULL"))
    }
  }
  
  # Combine all data frames
  final_df <- bind_rows(results)
  
  return(final_df)
}
normalize_count_df_long <- function(count_df_long, norm_method = "control_median"){
  allowed_norm_methods <- c("control_median")
  if (!(norm_method %in% allowed_norm_methods)){
    cat("ERROR: ", norm_method, "is not an allowed normalization method.\n")
    cat("Implemented normalization methods: ", allowed_norm_methods,"\n")
    stop()
  }
  
  if (norm_method == "control_median"){
    
    # 1. Determine normalization factors based on control categories
    norm_fac <- count_df_long %>%
      filter(group_category %in% c("targeting_control", "non_targeting_control", "kept_control")) %>%
      group_by(condition, sublib, sample) %>%
      summarise(norm_factor = median(count), .groups = "drop")
    
    # 2. Global median across all counts (entire dataset)
    med_count <- median(count_df_long$count)
    
    count_df_long_continue <- count_df_long
    
    if (med_count == 0 | any(norm_fac$norm_factor == 0)) {
      cat("WARNING: at least one normalization factor is 0\n")
      cat("Median of all counts:", med_count, "\n")
      
      zero_nf <- norm_fac %>% 
        filter(norm_factor == 0)
      
      if (!(nrow(zero_nf) == 0)) {
        cat("Groups with norm_factor == 0:\n")
        print(norm_fac)# shows condition | sublib | sample | norm_factor
        cat("Adding global pseudocount (+1).\n")
        pseudocount_added <<- TRUE
        
        count_df_long_plus_one <- count_df_long %>% mutate(count = count + 1)
        norm_fac <- count_df_long_plus_one %>%
          filter(group_category %in% c("targeting_control", "non_targeting_control", "kept_control")) %>%
          group_by(condition, sublib, sample) %>%
          summarise(norm_factor = median(count), .groups = "drop")
        
        # 2. Global median across all counts (entire dataset) #
        med_count <- median(count_df_long_plus_one$count)
        
        count_df_long_continue <- count_df_long_plus_one
      }
      
      cat("--------------------------------------------\n")
    }
    
    # 3. Normalize counts using (count * med_count / norm_factor)
    return_df <- count_df_long_continue %>%
      inner_join(norm_fac, by = c("condition", "sublib", "sample")) %>%
      mutate(norm_count = (count * med_count) / norm_factor) %>% 
      mutate(count = round(norm_count,2)) %>%
      select(-c(norm_factor,norm_count))
  }
  return_df <- return_df %>% mutate(count = round(count))
  return(return_df)
}
