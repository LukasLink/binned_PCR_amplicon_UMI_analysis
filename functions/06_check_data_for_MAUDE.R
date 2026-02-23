plot_maude_qc <- function(maude_counts_df, print = TRUE, add_0s = FALSE) {
  
  if (add_0s == TRUE){
    maude_counts_df <- maude_counts_df %>%
      mutate(
        upper = ifelse(is.na(upper), 0, upper),
        lower = ifelse(is.na(lower), 0, lower)
      )
  }
  # --- 1. Overall NA counts ---
  na_summary_overall <- maude_counts_df %>%
    pivot_longer(cols = c(input, lower, upper), names_to = "condition", values_to = "value") %>%
    mutate(is_na = is.na(value)) %>%
    group_by(condition) %>%
    summarise(na_count = sum(is_na), .groups = "drop")
  
  p1 <- ggplot(na_summary_overall, aes(x = condition, y = na_count, fill = condition)) +
    geom_col() +
    labs(title = "Missing Entries per Bin", y = "Number of NAs", x = "Condition") +
    theme_bw()
  
  # --- 2. NA counts by exp ---
  na_summary_by_exp <- maude_counts_df %>%
    pivot_longer(cols = c(input, lower, upper), names_to = "condition", values_to = "value") %>%
    mutate(is_na = is.na(value)) %>%
    group_by(exp, condition) %>%
    summarise(na_count = sum(is_na), .groups = "drop")
  
  p2 <- ggplot(na_summary_by_exp, aes(x = exp, y = na_count, fill = condition)) +
    geom_col(position = "dodge") +
    labs(title = "Missing Entries per Bin by Experiment", x = "Experiment", y = "Number of NAs") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # --- 3. sgRNA duplicate count distribution ---
  sgRNA_counts_per_exp <- maude_counts_df %>%
    group_by(exp, sgRNA) %>%
    summarise(dup_count = n(), .groups = "drop")
  
  p3 <- ggplot(sgRNA_counts_per_exp, aes(x = exp, y = dup_count)) +
    geom_boxplot(fill = "skyblue") +
    labs(
      title = "How Often Do Guides Appear Per Sample?",
      x = "Experiment",
      y = "Duplicated Guides"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (print == TRUE){
    print(p1)
    print(p2)
    print(p3)
  } else {
    return(list(
      overall_na_plot = p1,
      na_by_exp_plot = p2,
      sgRNA_duplicates_plot = p3
    ))
  }
  
  # Optionally return all plots
  
}

control_sanity_check <- function(maude_counts_df){
  # Create a list to store plots for each sublib
  plot_list <- list()
  
  # Loop through each unique sublib in the dataset
  for (sublib in unique(maude_counts_df$sublib)) {
    
    # Select data for the current sublib
    sublib_data <- maude_counts_df[maude_counts_df$sublib == sublib, ]
    
    # Select rows where sgRNA starts with "CONTROL_"
    control_data <- sublib_data[grepl("^CONTROL_", sublib_data$sgRNA), ]
    
    # Calculate log2 fold changes between upper/input and lower/input
    control_data$log2_upper_input <- log2(control_data$upper / control_data$input)
    control_data$log2_lower_input <- log2(control_data$lower / control_data$input)
    
    # Reshape the data for boxplot (long format with log2 fold change values)
    log2_values <- data.frame(
      sgRNA = rep(control_data$sgRNA, 2),
      log2_fold_change = c(control_data$log2_upper_input, control_data$log2_lower_input),
      comparison = rep(c("upper vs input", "lower vs input"), each = nrow(control_data))
    )
    
    # Create a boxplot for the current sublib
    p <- ggplot(log2_values, aes(x = comparison, y = log2_fold_change, fill = comparison)) +
      geom_boxplot(width = 0.5, alpha = 0.7) +
      scale_fill_manual(values = c("upper vs input" = "lightblue", "lower vs input" = "lightgreen")) +
      theme_bw() +
      labs(
        title = paste("Log2 Fold Change Between Upper/Lower and Input for", sublib),
        x = "",
        y = "Log2 Fold Change"
      ) +
      theme(
        axis.text.x = element_text(angle = 50, hjust = 1, face = "bold", size = 12),
        plot.title = element_text(size = 12),
        legend.position = "none"
      )
    
    # Store the plot for the current sublib
    plot_list[[sublib]] <- p
    
    # Print the plot for the current sublib
    print(p)
  }
  
  # Also, calculate for the whole dataset (without sublib filtering)
  # Select rows where sgRNA starts with "CONTROL_"
  control_data_all <- maude_counts_df[grepl("^CONTROL_", maude_counts_df$sgRNA), ]
  
  # Calculate log2 fold changes between upper/input and lower/input for the entire dataset
  control_data_all$log2_upper_input <- log2(control_data_all$upper / control_data_all$input)
  control_data_all$log2_lower_input <- log2(control_data_all$lower / control_data_all$input)
  
  # Reshape the data for boxplot (long format with log2 fold change values)
  log2_values_all <- data.frame(
    sgRNA = rep(control_data_all$sgRNA, 2),
    log2_fold_change = c(control_data_all$log2_upper_input, control_data_all$log2_lower_input),
    comparison = rep(c("upper vs input", "lower vs input"), each = nrow(control_data_all))
  )
  
  # Create a boxplot for the entire dataset
  p_all <- ggplot(log2_values_all, aes(x = comparison, y = log2_fold_change, fill = comparison)) +
    geom_boxplot(width = 0.5, alpha = 0.7) +
    scale_fill_manual(values = c("upper vs input" = "lightblue", "lower vs input" = "lightgreen")) +
    theme_bw() +
    labs(
      title = "Log2 Fold Change Between Upper/Lower and Input for All Data",
      x = "",
      y = "Log2 Fold Change"
    ) +
    theme(
      axis.text.x = element_text(angle = 50, hjust = 1, face = "bold", size = 12),
      plot.title = element_text(size = 12),
      legend.position = "none"
    )
  
  # Print the plot for the entire dataset
  print(p_all)
}