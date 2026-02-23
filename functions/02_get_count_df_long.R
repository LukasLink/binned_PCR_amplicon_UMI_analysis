#===============================================================================
# get_count_df_long
#===============================================================================
read_file_to_df <- function(file_name,
                            folder_path = dedup_output_folder,
                            suffix_to_rm = "_dedup_idxstats.txt",
                            threshold_df = NULL,
                            check_alignments = TRUE
){
  
  file_path <- file.path(folder_path,file_name)  
  df <- read.table(file_path, header=FALSE, sep="\t", stringsAsFactors=FALSE)
  name <- sub(suffix_to_rm,"",file_name)
  sub_lib <- str_match(name, "^[A-Za-z]+_L(\\d+)_[^_]+")[, 2]
  sub_lib <- paste0("ICS_",sub_lib)
  exp <- str_replace(name, "^([^_]+_)", "")
  if(is.null(threshold_df)){
    threshold <- 0
  } else {
    threshold <- threshold_df$threshold[threshold_df$replicate == name]
  }
  
  df <- df %>% 
    rename(sgRNA = V1, count = V3, length = V2) %>%
    select(-V4) %>%
    mutate(count = if_else(count <= threshold, 0, count),
           exp = exp) %>%
    filter(sgRNA != "*")
  
  check_df <- merged_sgRNA_df %>% select(sgrna_id, sublib)
  # First try: plain left_join
  joined_df <- df %>%
    left_join(check_df, by = c("sgRNA" = "sgrna_id"))
  
  # Check for non-unique sgRNA in the result
  if (any(duplicated(joined_df$sgRNA))) {
    warning("Non-unique `sgRNA` after join; redoing join with 1:1 pairing by order.\n")
    
    # Add within-id index on both sides
    df_idx <- df %>%
      group_by(sgRNA) %>%
      mutate(.pair_id = row_number()) %>%
      ungroup()
    
    check_idx <- check_df %>%
      group_by(sgrna_id) %>%
      mutate(.pair_id = row_number()) %>%
      ungroup()
    
    # Redo join using sgRNA/sgrna_id + the position
    joined_df <- df_idx %>%
      left_join(
        check_idx,
        by = c("sgRNA" = "sgrna_id", ".pair_id")
      ) %>%
      select(-.pair_id)      # clean up helper column
    df <- joined_df
  } else {
    df <- joined_df
  }
  
  if (check_alignments == TRUE) {
    cat("Checking wrong allignments for: ", name,"\n")
    # Filter non-control rows with count > 0
    non_control_df <- df %>%
      filter(!grepl("^CONTROL_", sgRNA), count > 0)
    
    # Count sum stats
    correct_sum <- non_control_df %>% filter(sublib == sub_lib) %>% summarise(total = sum(count)) %>% pull(total)
    wrong_sum <- non_control_df %>% filter(sublib != sub_lib) %>% summarise(total = sum(count)) %>% pull(total)
    total_sum <- correct_sum + wrong_sum
    
    cat("UMI/read count allignment Stats:\n")
    cat(sprintf("  Correct (aligned to sgRNA in sublibary):     %d (%.2f%%)\n", 
                correct_sum, 100 * correct_sum / total_sum))
    cat(sprintf("  Wrong   (aligned to sgRNA not in sublibary): %d (%.2f%%)\n\n", 
                wrong_sum, 100 * wrong_sum / total_sum))
  }
  
  same_controls_in_all_sublibraries <- get("same_controls_in_all_sublibraries", envir = .GlobalEnv)
  if (same_controls_in_all_sublibraries == TRUE){
    df <- df %>%
      filter(sublib == sub_lib | grepl("^CONTROL_", sgRNA)) %>%
      select(-sublib)
  } else {
    df <- df %>%
      filter(sublib == sub_lib) %>%
      select(-sublib)
  }
  
  if (check_alignments == TRUE) {
    cat(sprintf("sgRNA coverage: %.2f%%\n", sum(df$count > 0) / nrow(df) * 100))
    cat("--------------------------------------------\n")
  }
  
  df <- df %>%
    mutate(group_category = case_when(
      # grepl("^CONTROL_C_EGFP_", sgRNA) ~ "EGFP_control",
      grepl("^CONTROL_C_NONTARG_", sgRNA) ~ "non_targeting_control",
      grepl("^CONTROL_C_", sgRNA) ~ "targeting_control",
      TRUE ~ "targeting"
    ))
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
read_john_data <- function(dir_path, sgRNA_df) {
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
read_john_rf_data  <- function(file_path, sgRNA_df){
  count_df <- read.table(file_path, sep = "\t", header = FALSE, col.names = c("seq", "count"))
  
  joined_df <- sgRNA_df %>%
    select(seq, sgrna_id) %>% 
    inner_join(count_df, by = "seq") %>% 
    select(sgrna_id,count) %>% 
    mutate(count = as.integer(count))
  
  return(joined_df)
}

# Helper function for process_john_rf_data and process_john_data
complete_sgrna_df <- function(df,
                              merged_sgRNA_df,
                              sublib,
                              name = NULL,
                              check_alignments = TRUE) {
  
  expected_ids <- merged_sgRNA_df %>%
    filter(sublib == !!sublib | grepl("^CONTROL_", sgrna_id)) %>%
    pull(sgrna_id)
  
  missing_ids <- base::setdiff(expected_ids, df$sgrna_id)
  
  new_rows <- tibble(sgrna_id = missing_ids, count = 0)
  
  complete_df <- bind_rows(df, new_rows)
  
  if (check_alignments == TRUE) {
    if (is.null(name)) {name <- sublib}
    
    # Join with merged_sgRNA_df to get sublib info
    annotated_df <- complete_df %>%
      left_join(merged_sgRNA_df %>% select(sgrna_id, sublib), by = "sgrna_id") %>%
      filter(!grepl("^CONTROL_", sgrna_id), count > 0)
    
    correct_sum <- annotated_df %>% filter(sublib == !!sublib) %>% summarise(total = sum(count)) %>% pull(total)
    wrong_sum   <- annotated_df %>% filter(sublib != !!sublib) %>% summarise(total = sum(count)) %>% pull(total)
    total_sum   <- correct_sum + wrong_sum
    
    cat("Checking wrong alignments for: ", name, "\n")
    cat("UMI/read count alignment Stats:\n")
    cat(sprintf("  Correct (aligned to sgRNA in sublibrary):     %d (%.2f%%)\n", 
                correct_sum, 100 * correct_sum / total_sum))
    cat(sprintf("  Wrong   (aligned to sgRNA not in sublibrary): %d (%.2f%%)\n\n", 
                wrong_sum, 100 * wrong_sum / total_sum))
    cat("--------------------------------------------\n")
  }
  
  
  return(complete_df)
}

process_john_rf_data <- function(parent_dir,
                                 merged_sgRNA_df = merged_sgRNA_df,
                                 data_type = "umis",
                                 skip_list = c("L_L4_3","U_L4_3"),
                                 suffix_to_rm = "_results.tsv",
                                 check_alignments = TRUE) {
  
  # Initialize list to collect data
  results <- list()
  file_paths <- list.files(path = parent_dir, pattern = "_results\\.tsv$", full.names = TRUE)
  merged_sgRNA_df$sublib <- gsub("ICS", "sublib", merged_sgRNA_df$sublib)
  for (file in file_paths){
    file_name <- basename(file)
    name <- sub(suffix_to_rm,"",file_name)
    # Skip entries from skip list
    skip_this_file <- any(vapply(skip_list, function(x) {
      if (startsWith(x, "re:")) {
        grepl(sub("^re:", "", x), name, perl = TRUE)
      } else {
        identical(x, name)
      }
    }, logical(1)))
    
    if (skip_this_file){next}
    
    # Extract condition (X), sublib (L#), and sample (#) using regex
    pattern <- paste0("^([ILU])_L(\\d+)_([^_]+)", suffix_to_rm)
    matches <- str_match(file_name, pattern)
    
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
    exp <- paste0("L", matches[3], "_", matches[4])
    
    # Safely read data
    df <- tryCatch({
      df <- read_john_rf_data(file, merged_sgRNA_df)
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
      cat("Skipping due to error in", file_name, "->", e$message, "\n")
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
    }
  }
  
  # Combine all data frames
  final_df <- bind_rows(results)
  
  return(final_df)
}

process_john_data <- function(parent_dir,
                              merged_sgRNA_df = merged_sgRNA_df,
                              data_type = "umis",
                              skip_list = c("L_L4_3","U_L4_3"),
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
      df <- read_john_data(matrix_dir, merged_sgRNA_df)
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
