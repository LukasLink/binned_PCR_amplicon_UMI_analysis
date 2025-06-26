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

# Change this to your custom thresholds, and file names
threshold_df <- data.frame(
  replicate = c("I_L1_1", "I_L1_2", "I_L1_3", "I_L2_1", "I_L2_2", "I_L3_1", "I_L3_2", "I_L4_1", "I_L4_2",
                "L_L1_1", "L_L1_2", "L_L1_3", "L_L2_1", "L_L2_2", "L_L3_1", "L_L3_2", "L_L3_3", "L_L4_1",
                "L_L4_2", "L_L4_3", "U_L1_1", "U_L1_2", "U_L1_3", "U_L2_1", "U_L2_2", "U_L3_1", "U_L3_2",
                "U_L3_3", "U_L4_1", "U_L4_2", "U_L4_3"),
  threshold = c(3000, 2000, 2000, 3000, 3000, 2000, 3000, 3000, 3000,
                3000, 1500, 3000, 3000, 2000, 30000, 2000, 2000, 3000,
                3000, 2000, 1000, 1000, 3000, 2000, 800, 1000, 2000,
                1000, 1000, 1000, 2000)
)

#######

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

filter_UMIs_with_few_reads = function(threshold_df,
                                      folder_path = grouped_output_folder,
                                      count_only = FALSE){
  
  print("Beginning Generation of Read Filter *_threshold_umis.txt files")
  file_paths <- list.files(path = folder_path, pattern = "\\.tsv$", full.names = TRUE)
  file_names <- list.files(path = folder_path, pattern = "\\.tsv$", full.names = FALSE)
  file_names = gsub("\\.tsv$", "", file_names)
  graph_list = list()
  excluded_reads_total <- 0
  for (i in 1:length(file_paths)){
    
    excluded_reads_sample <- 0
    df = read.table(file_paths[[i]], sep="\t", header=TRUE)
    df <- df %>% select(read_id,final_umi_count)
    # Extract the value from threshold_df$threshold
    threshold_value <- threshold_df %>%
      filter(replicate == file_names[[i]]) %>%
      pull(threshold)
    # Filter rows where final_umi_count is less than or equal to X
    df <- df %>%
      filter(final_umi_count <= threshold_value)
    
    excluded_reads_sample <- nrow(df)
    excluded_reads_total <- excluded_reads_total + excluded_reads_sample
    if (count_only == TRUE){
      save_path = sub(".tsv","_threshold_umis.txt",file_paths[[i]])
      write.table(df$read_id, save_path, row.names = FALSE, col.names = FALSE, quote = FALSE)
    }

    print(paste("Completed",i,"out of",length(file_paths)))
    print(paste("Excluded reads for sample:",excluded_reads_sample))
  }
  print(paste("TOTAL excluded reads:",excluded_reads_total))
  return(df)
}

UMI_read_filter_df = filter_UMIs_with_few_reads(threshold_df, folder_path = grouped_output_folder, count_only = TRUE)