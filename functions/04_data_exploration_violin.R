# These functions are for creating violin plots for initial data inspection.
plot_violin_by_sublib_sample <- function(count_df_long, norm_method = NULL) {
  df <- count_df_long
  
  if (!is.null(norm_method)){
    df <- normalize_count_df_long(count_df_long, norm_method)
    norm_title <- paste0("(normalization: ",norm_method,")")
  } else {
    norm_title <- "(no normalization)"
  }
  count_type <- get("data_type", envir = .GlobalEnv)
  if (!(count_type %in% c("umis","reads"))){
    cat("ERROR ", count_type, " is not a viable data_type.")
    cat("data_type must be umis or reads")
    stop()
  }
  if (count_type == "umis"){
    title_prefix <- "UMI counts for"
  }
  if (count_type == "reads"){
    title_prefix <- "Read counts for"
  }
  # Use exact matching for group_type
  df <- df %>%
    mutate(group_type = ifelse(group_category == "targeting", "targeting", "non-targeting"))
  
  # Get all unique (sublib, sample) pairs
  unique_groups <- df %>%
    distinct(sublib, sample)
  
  plot_list <- list()
  
  for (i in 1:nrow(unique_groups)) {
    
    sublib_val <- unique_groups$sublib[i]
    sample_val <- unique_groups$sample[i]
    title <- paste0("L",sub("sublib_","",sublib_val),"_",sub("sample_","",sample_val))
    # Subset for that combination
    df_subset <- df %>%
      filter(sublib == sublib_val, sample == sample_val)
    
    # Filter 3 unique conditions from each group_type
    df_subset <- df_subset %>%
      rename(Type = group_type) %>% 
      group_by(Type, condition) %>%
      filter(n() > 0) %>%
      ungroup()
    
    targeting_conds <- df_subset %>%
      filter(Type == "targeting") %>%
      pull(condition) %>%
      unique() %>%
      head(3)
    
    nontargeting_conds <- df_subset %>%
      filter(Type == "non-targeting") %>%
      pull(condition) %>%
      unique() %>%
      head(3)
    
    df_filtered <- df_subset %>%
      filter((Type == "targeting" & condition %in% targeting_conds) |
               (Type == "non-targeting" & condition %in% nontargeting_conds))
    
    # Plot
    p <- ggplot(df_filtered, aes(x = condition, y = count, fill = Type)) +
      geom_violin(trim = FALSE, scale = "width") +
      labs(
        title = paste0(title_prefix," ",sublib_val,", ", sample_val," ",norm_title),
        x = "",
        y = "Count"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    plot_list[[paste(sublib_val, sample_val, sep = "_")]] <- p
  }
  
  return(plot_list)
}
get_grouped_summary_wide <- function(count_df_long, stat = c("median", "mean")) {
  stat <- match.arg(stat)  # ensure valid input
  
  # Define function to apply
  summary_fun <- switch(
    stat,
    median = function(x) median(x, na.rm = TRUE),
    mean = function(x) round(mean(x, na.rm = TRUE), 2)
  )
  
  # Assign group_type via exact matching
  count_df_long <- count_df_long %>%
    mutate(group_type = ifelse(group_category == "targeting", "targeting", "non-targeting"))
  
  unique_groups <- count_df_long %>%
    distinct(sublib, sample)
  
  summary_rows <- list()
  
  for (i in 1:nrow(unique_groups)) {
    sublib_i <- unique_groups$sublib[i]
    sample_i <- unique_groups$sample[i]
    
    df_sub <- count_df_long %>%
      filter(sublib == sublib_i, sample == sample_i)
    
    for (grp in c("targeting", "non-targeting")) {
      df_grp <- df_sub %>%
        filter(group_type == grp)
      
      # Take up to 3 sorted conditions
      conds <- sort(unique(df_grp$condition))[1:min(3, length(unique(df_grp$condition)))]
      
      # Calculate summary values (median or mean)
      summaries <- df_grp %>%
        filter(condition %in% conds) %>%
        group_by(condition) %>%
        summarise(val = summary_fun(count), .groups = "drop") %>%
        arrange(condition)
      
      values <- c(summaries$val, rep(NA, 3 - nrow(summaries)))
      names(values) <- c("condition_lower", "condition_middle", "condition_upper")[1:length(values)]
      
      row <- tibble(
        group_id = paste(sublib_i, sample_i, grp, sep = "_"),
        !!!as.list(values)
      )
      
      summary_rows[[length(summary_rows) + 1]] <- row
    }
  }
  
  summary_df <- bind_rows(summary_rows)
  return(summary_df)
}


plot_violin_by_group_category_split_by_sublib <- function(df,
                                                          include_targeting = TRUE,
                                                          norm_method = NULL,
                                                          y_limit = 60,
                                                          box_col = "white",
                                                          viol_col = "#E0E0E0") {
  
  
  
  
  if (!is.null(norm_method)){
    df <- normalize_count_df_long(df, norm_method)
    norm_title <- paste0("(normalization: ",norm_method,")")
  } else {
    norm_title <- "(no normalization)"
  }
  count_type <- get("data_type", envir = .GlobalEnv)
  if (!(count_type %in% c("umis","reads"))){
    cat("ERROR ", count_type, " is not a viable data_type.")
    cat("data_type must be umis or reads")
    stop()
  }
  if (count_type == "umis"){
    title_prefix <- "UMI counts for"
  }
  if (count_type == "reads"){
    title_prefix <- "Read counts for"
  }
  
  # Filter based on group_category
  filtered_df <- df %>%
    filter(
      if (include_targeting) group_category == "targeting"
      else group_category != "targeting"
    ) %>%
    mutate(group = interaction(condition, exp, drop = TRUE))
  
  # Compute count summary (output + for annotation)
  count_summary <- filtered_df %>%
    group_by(condition, exp) %>%
    summarise(
      total_count = sum(count),
      mean_count = round(mean(count), 1),
      sd_count = round(sd(count), 1),
      .groups = "drop"
    ) %>%
    mutate(group = interaction(condition, exp, drop = TRUE))
  
  # Create one plot per unique sublib
  plots <- list()
  
  for (sublib_name in unique(filtered_df$sublib)) {
    sub_df <- filtered_df %>% filter(sublib == sublib_name)
    # Levels used by this sub-plot’s x axis
    x_levels <- levels(interaction(sub_df$condition, sub_df$exp, drop = TRUE))
    
    # Build label_df with matching factor levels
    label_df <- count_summary %>%
      filter(interaction(condition, exp, drop = TRUE) %in%
               unique(interaction(sub_df$condition, sub_df$exp, drop = TRUE))) %>%
      mutate(
        group_fac = factor(interaction(condition, exp, drop = TRUE), levels = x_levels),
        x_center  = as.numeric(group_fac)   # centers now 1..k, aligned with this sub-plot
      )
    
    # Plot with violin, boxplot, and text
    p <- ggplot(sub_df, aes(x = interaction(condition, exp), y = count)) +
      geom_violin(fill = viol_col, color = NA, scale = "width", trim = FALSE) +
      geom_boxplot(width = 0.1, outlier.size = 0.2, fill = box_col) +
      geom_rect(
        data = label_df,
        aes(xmin = x_center - 0.3, xmax = x_center + 0.3,
            ymin = y_limit - y_limit/8 - 0.1, ymax = y_limit + 0.1),
        inherit.aes = FALSE,
        fill = "white", color = "black", linewidth = 0.15
      ) +
      geom_text(
        data = label_df,
        aes(x = group, y = y_limit - y_limit/26, label = total_count),
        size = 2.5, inherit.aes = FALSE
      ) +
      geom_text(
        data = label_df,
        aes(x = group, y = y_limit - y_limit/14, label = paste0("mean: ", mean_count)),
        size = 2.3, inherit.aes = FALSE
      ) +
      geom_text(
        data = label_df,
        aes(x = group, y = y_limit - y_limit/10, label = paste0("sd: ", sd_count)),
        size = 2.3, inherit.aes = FALSE
      ) +
      labs(
        title = paste(title_prefix,
                      ifelse(include_targeting, "targeting sgRNA", "non-targeting sgRNA"),
                      "–",
                      sublib_name, norm_title),
        x = "",
        y = "Count"
      ) +
      coord_cartesian(ylim = c(0, y_limit)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    plots[[sublib_name]] <- p
  }
  
  return(list(
    plots = plots,
    count_summary = count_summary
  ))
}
generate_all_violin_plots_and_summaries <- function(df,
                                                    norm_method=NULL,
                                                    targeting = TRUE,
                                                    non_targeting = FALSE,
                                                    summary_df = FALSE,
                                                    y_limit = 60) {
  # --- Print all 8 plots ---
  if (targeting == TRUE){
    # Run plotting function for targeting
    targeting_results <- plot_violin_by_group_category_split_by_sublib(df,
                                                                       include_targeting = TRUE,
                                                                       norm_method=norm_method,
                                                                       y_limit = y_limit,
                                                                       viol_col = "lightgreen")
    cat("Targeting plots:\n")
    for (sublib in names(targeting_results$plots)) {
      print(targeting_results$plots[[sublib]])
    }
    if (summary_df == TRUE){
      # --- Print summary tables ---
      cat("\nTargeting count summary:\n")
      print(targeting_results$count_summary)
    }
    
  }
  
  if (non_targeting == TRUE){
    non_targeting_results <- plot_violin_by_group_category_split_by_sublib(df,
                                                                           include_targeting = FALSE,
                                                                           norm_method=norm_method,
                                                                           y_limit = y_limit,
                                                                           viol_col = "lightblue")
    cat("\nNon-targeting plots:\n")
    for (sublib in names(non_targeting_results$plots)) {
      print(non_targeting_results$plots[[sublib]])
    }
    if (summary_df == TRUE){
      cat("\nNon-targeting count summary:\n")
      print(non_targeting_results$count_summary)
    }
  }
}