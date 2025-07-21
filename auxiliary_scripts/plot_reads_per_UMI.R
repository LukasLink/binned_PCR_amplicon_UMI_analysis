library(tidyverse)
library(conflicted)
library(Biostrings)
library(readxl)
library(conflicted)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::count)

# Change this to match the output folder
########
output_folder <- "/g/steinmetz/link/Amplicon_barcode_analysis/NB_EXP030/test"
########

make_clean_dir <- function(base_path, sub_path) {
  
  full_path <- paste0(base_path, "/", sub_path)
  
  full_path <- gsub("///", "/", full_path, fixed = TRUE)
  full_path <- gsub("//", "/", full_path, fixed = TRUE)
  
  if (!dir.exists(full_path)) {
    dir.create(full_path, recursive = TRUE)
  }
  return(full_path)
}
get_file_path <- function(folder_path, file_name){
  full_path <- paste0(folder_path, "/",file_name)
  full_path <- gsub("///", "/", full_path, fixed = TRUE)
  full_path <- gsub("//", "/", full_path, fixed = TRUE)
  return(full_path)
}

genome_output_folder <- make_clean_dir(output_folder, "/ref/NOPE/")
grouped_output_folder <- make_clean_dir(output_folder, "/grouped/")
rds_output_folder <- make_clean_dir(output_folder, "/rds/")


plot_reads_per_UMI = function(threshold_df_short = NULL, folder_path = grouped_output_folder){
  
  file_paths <- list.files(path = folder_path, pattern = "\\.tsv$", full.names = TRUE)
  file_names <- list.files(path = folder_path, pattern = "\\.tsv$", full.names = FALSE)
  file_names = gsub("\\.tsv$", "", file_names)
  graph_list = list()
  for (i in 1:length(file_paths)){
    # if (i > 1){
    #   next
    # }
    name = file_names[[i]]
    title = paste0("Knee Plot of Reads per UMI: ",name)
    df_1 = read.table(file_paths[[i]], sep="\t", header=TRUE)
    df_1 <- df_1 %>% select(contig,umi,umi_count,final_umi_count,unique_id)
    # Remove duplicate (unique_id, umi) pairs
    df_1 <- df_1 %>%
      distinct(unique_id, umi, .keep_all = TRUE) 
    df_2 <- df_1 %>% distinct(unique_id, .keep_all = TRUE)
    # Create a ranked index (1 = largest count)
    # Rank both datasets (highest value = rank 1)
    df_1 <- df_1 %>% mutate(rank = rank(-umi_count, ties.method = "first"))
    df_2 <- df_2 %>% mutate(rank = rank(-final_umi_count, ties.method = "first"))
    
    # Convert to a dataframe for plotting
    # Combine datasets into a single long format dataframe for ggplot
    df <- bind_rows(
      df_1 %>% rename(reads_per_umi = umi_count) %>% mutate(line = "pre-grouping"),
      df_2 %>% rename(reads_per_umi = final_umi_count) %>% mutate(line = "post-grouping")
    )
    
    p = ggplot(df, aes(x = rank, y = reads_per_umi, color = line)) +
      geom_line() + 
      scale_x_log10() +  # Logarithmic x-axis
      scale_y_log10() +
      labs(title = title,
           x = "UMI rank", 
           y = "Reads per UMI (log10)") +
      theme_bw()
    
    if (!is.null(threshold_df_short = NULL)){
      p <- p +
        geom_hline(aes(yintercept = threshold_df_short$threshold[threshold_df_short$replicate == name], 
                       linetype = "Threshold"), 
                   color = "red", size = 1) +
        scale_linetype_manual(name = "", values = c("Threshold" = "dashed")) # Add threshold to legend
    }
    
    graph_list[[i]] <- p
    
  }
  return(graph_list)
}

read_count_plots = plot_reads_per_UMI()
saveRDS(read_count_plots,
        get_file_path(rds_output_folder, "read_count_plots.rds"))
