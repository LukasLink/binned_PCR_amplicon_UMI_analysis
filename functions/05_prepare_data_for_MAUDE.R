#Select wether to use targeting_control or non_targeting_control for normalization (or both I guess)
ctrl_selection <- c("targeting_control","non_targeting_control")
subsample_controls_func_old <- function(count_df_long, merged_sgRNA_df, percentage_kept_controls = 0.1){
  
  percentage_used_controls = 1 - percentage_kept_controls
  count_df_long_merged <- merge(count_df_long, merged_sgRNA_df[, c("sgrna_id", "entrez")], 
                                by.x = "sgRNA", by.y = "sgrna_id", 
                                all.x = TRUE)
  count_df_long_merged <- count_df_long_merged %>%
    filter(group_category != "targeting") %>%
    distinct(sgRNA) %>% 
    mutate(
      unique_names = case_when(
        grepl("CONTROL_C_[A-Za-z0-9]+_", sgRNA) & !grepl("CONTROL_C_NONTARG", sgRNA) ~ sub("CONTROL_C_([A-Za-z0-9]+)_.*", "\\1", sgRNA),
        grepl("CONTROL_C_NONTARG", sgRNA) ~ sgRNA,  # Keep NONTARG as it is
        TRUE ~ NA_character_  # Default case, if needed
      )
    )
  unique_entries <- unique(count_df_long_merged$unique_names)
  assignments <- sample(c(TRUE, FALSE), length(unique_entries), 
                        prob = c(percentage_kept_controls, percentage_used_controls), replace = TRUE)
  # Step 3: Create a data frame to map unique entries to the TRUE/FALSE assignments
  assignment_map <- data.frame(unique_names = unique_entries, keep = assignments)
  # Step 4: Merge this map back into the original data frame to expand the assignments
  count_df_long_merged <- count_df_long_merged %>%
    left_join(assignment_map, by = "unique_names")
  return_df <- count_df_long %>%
    left_join(count_df_long_merged %>% select(sgRNA, keep), by = "sgRNA")
  # Replace NAs in 'keep' with FALSE
  return_df <- return_df %>%
    mutate(keep = ifelse(is.na(keep), FALSE, keep)) %>%
    mutate(
      # Assign 'kept_control' to group_category when keep is TRUE
      group_category = ifelse(keep == TRUE, "kept_control", group_category)
    ) %>%
    select(-keep)  # Remove the 'keep' column
  
  return(return_df)
}
subsample_controls_func <- function(count_df_long, merged_sgRNA_df, percentage_kept_controls = 0.2){
  
  percentage_used_controls = 1 - percentage_kept_controls
  
  # Merge count_df_long with merged_sgRNA_df to include 'entrez'
  count_df_long_merged <- merge(count_df_long, merged_sgRNA_df[, c("sgrna_id", "entrez")], 
                                by.x = "sgRNA", by.y = "sgrna_id", 
                                all.x = TRUE)
  
  # Filter out the "targeting" group
  count_df_long_merged <- count_df_long_merged %>%
    filter(group_category != "targeting")
  
  # Select only unique sgRNA entries (no duplicates) for subsampling
  unique_entries <- unique(count_df_long_merged$sgRNA)
  
  # Randomly assign TRUE/FALSE based on percentage_kept_controls
  assignments <- sample(c(TRUE, FALSE), length(unique_entries), 
                        prob = c(percentage_kept_controls, percentage_used_controls), replace = TRUE)
  
  # Create a data frame that maps unique sgRNAs to their keep assignments
  assignment_map <- data.frame(sgRNA = unique_entries, keep = assignments)
  
  # Merge assignment_map with the original count_df_long
  return_df <- count_df_long %>%
    left_join(assignment_map, by = "sgRNA") %>%
    mutate(
      # Change group_category based on keep and conditions
      group_category = ifelse(keep == TRUE & group_category != "targeting", "kept_control", group_category)
    ) %>%
    select(-keep)  # Drop the 'keep' column after mutation
  
  return(return_df)
}

count_df_long_to_wide <- function(count_df_long,
                                  print = TRUE,
                                  drop_0s = TRUE,
                                  recover_input = FALSE,
                                  for_sum = FALSE){
  if (print == TRUE){
    print(paste("NAs in count_df_long:",sum(is.na(count_df_long))))
  }
  if (for_sum == TRUE){
    # recover_input <- FALSE
  }
  maude_counts_df <- count_df_long %>%
    mutate(isNontargeting = ifelse(group_category %in% ctrl_selection, T, F)) %>%  
    select(-c(group_category)) %>%
    pivot_wider(names_from = condition, values_from = count) %>%
    mutate(
      upper = ifelse(is.na(upper) & !is.na(lower), 0, upper), # Replace NA in upper with 0
      lower = ifelse(is.na(lower) & !is.na(upper), 0, lower)  # Replace NA in lower with 0
    ) %>% 
    as.data.frame()
  if (print == TRUE){
    print(paste("Dimension of wide_df:",dim(maude_counts_df)))
    print(paste("NA's in maude_counts_df:", sum(is.na(maude_counts_df))))
    print(paste("NA's in Input:", sum(is.na(maude_counts_df$input))))
    print(paste("NA's in Upper:", sum(is.na(maude_counts_df$upper))))
    print(paste("NA's in Lower:", sum(is.na(maude_counts_df$lower))))
  } 
  
  if (recover_input == TRUE){
    # Step 2: Impute 'input' values where both upper and lower exist
    # First, compute the average input for each sgRNA from non-NA values
    input_means <- maude_counts_df %>%
      group_by(sgRNA, sublib) %>%
      summarize(mean_input = mean(input, na.rm = TRUE), .groups = "drop")
    
    # Step 3: Join back the means and impute missing 'input' values if upper and lower exist
    maude_counts_df <- maude_counts_df %>%
      left_join(input_means, by = c("sgRNA","sublib")) %>%
      mutate(
        input = ifelse(is.na(input) & !is.na(upper) & !is.na(lower), mean_input, input)
      ) %>%
      select(-mean_input)  # clean up
  } else {
    if (for_sum == FALSE){
      maude_counts_df <- maude_counts_df %>% 
        drop_na() # Drop any row with NA that were introduced by pivot wider
    }
  }
  
  if (drop_0s == TRUE){
    maude_counts_df <- maude_counts_df %>%
      filter(!(input == 0 & upper == 0 & lower == 0))
  }
  
  if (for_sum == FALSE){
    maude_counts_df <- maude_counts_df %>% 
      mutate(input = input + 1,
             upper = upper + 1,
             lower = lower + 1)
  }
  
  if (print == TRUE){
    print(paste("Dimension of wide_df:",dim(maude_counts_df)))
    print(paste("NA's in maude_counts_df:", sum(is.na(maude_counts_df))))
    print(paste("NA's in Input:", sum(is.na(maude_counts_df$input))))
    print(paste("NA's in Upper:", sum(is.na(maude_counts_df$upper))))
    print(paste("NA's in Lower:", sum(is.na(maude_counts_df$lower))))
  }
  # MAUDE seems to not like floats, so we round to integers
  # TO DO: check how I can make maude use floats.
  maude_counts_df <- maude_counts_df %>% 
    mutate(input = round(input),
           lower = round(lower),
           upper = round(upper))
  return(maude_counts_df)
}
