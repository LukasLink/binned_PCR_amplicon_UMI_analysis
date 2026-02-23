# R/project_setup.R

project_setup <- function(project_root_dir,
                          user_opts,
                          envir = .GlobalEnv,
                          use_old_suffix_construction = FALSE) {
  stopifnot(is.character(project_root_dir), length(project_root_dir) == 1)
  stopifnot(is.list(user_opts))
  
  # ---- global knitr / options ----
  knitr::opts_chunk$set(echo = FALSE)
  options(bitmapType = "cairo")
  
  # ---- packages ----
  library(tidyverse)
  library(Matrix)
  library(conflicted)
  library(MAUDE)
  library(ggplot2)
  library(optparse)
  library(ggrepel)
  library(writexl)
  library(stringr)
  
  conflicts_prefer(dplyr::rename)
  conflicts_prefer(dplyr::filter)
  conflicts_prefer(dplyr::select)
  conflicts_prefer(dplyr::slice)
  
  #===============================================================================
  option_list <- list(
    make_option(c("--first_time"), type = "logical", default = user_opts$first_time,
                help = "First run flag, will load results from previous run with same settings if set to FALSE (default: TRUE)", metavar = "LOGICAL"),
    make_option(c("--output_folder"), type = "character", default = user_opts$output_folder,
                help = "Output folder",
                metavar = "PATH"),
    make_option(c("--john_folder"), type = "character", default = user_opts$john_folder,
                help = "John folder",
                metavar = "PATH"),
    make_option(c("--john_rf_folder"), type = "character", default = user_opts$john_rf_folder,
                help = "John RF folder",
                metavar = "PATH"),
    make_option(c("--skip_list"), type = "character", default = paste(user_opts$skip_list, collapse = ","),
                help = "If any files are to be skipped list their names here (comma-separated) (default: '')",
                metavar = "LIST"),
    make_option(c("--skip_list_sublib"), type = "character", default = paste(user_opts$skip_list_sublib, collapse = ","),
                help = "If any sublibraries are to be skipped list their names here (comma-separated) (default: '')",
                metavar = "LIST"),
    make_option(c("--skip_list_sample"), type = "character", default = paste(user_opts$skip_list_sample, collapse = ","),
                help = "If any samples are to be skipped list their names here (comma-separated) (default: '')",
                metavar = "LIST"),
    make_option(c("--include_controls_list"), type = "character", default = paste(user_opts$include_controls_list, collapse = ","),
                help = "If any control guides should be treated like they are targeting list them here. \nNames without a '_' will be applied to all versions with an '_' so 'EGFP' will affect both 'EGFP_1' and 'EGFP_2'",
                metavar = "LIST"),
    make_option(c("--pipeline"), type = "character", default = user_opts$pipeline,
                help = "Pipeline name (default: 'lukas')", metavar = "CHARACTER"),
    make_option(c("--data_type"), type = "character", default = user_opts$data_type,
                help = "Data type (default: 'umis')", metavar = "CHARACTER"),
    make_option(c("--method"), type = "character", default = user_opts$method,
                help = "Method (default: '')", metavar = "CHARACTER"),
    make_option(c("--norm_method"), type = "character", default = user_opts$norm_method,
                help = "Normalization method (default: '')", metavar = "CHARACTER"),
    make_option(c("--recover_input"), type = "logical", default = user_opts$recover_input,
                help = "Recover missing input (default: TRUE)", metavar = "LOGICAL"),
    make_option(c("--subsample_controls"), type = "logical", default = user_opts$subsample_controls,
                help = "Subsample control guides so some appear in the results (default: TRUE)", metavar = "LOGICAL"),
    make_option(c("--use_custom_bins"), type = "logical", default = user_opts$use_custom_bins,
                help = "Overrides upper_lower_percentage by using a custom binning system. (default: FALSE)", metavar = "LOGICAL"),
    make_option(c("--same_controls_in_all_sublibraries"), type = "logical", default = user_opts$same_controls_in_all_sublibraries,
                help = "Set to FALSE if the individual sublibraries/replicates have different control sgRNAs (default: TRUE)", metavar = "LOGICAL"),
    make_option(c("--extra_suffix"), type = "character", default = user_opts$extra_suffix,
                help = "Suffix for additional options (default: '')", metavar = "CHARACTER"),
    make_option(c("--simplified_rf_threshold"), type = "integer", default = user_opts$simplified_rf_threshold,
                help = "Threshold for simplified filtering (default: 1000)", metavar = "INTEGER"),
    make_option(c("--upper_lower_percentage"), type = "double", default = user_opts$upper_lower_percentage,
                help = "Fraction of the lower&upper bin (default: 0.10)", metavar = "DOUBLE"),
    make_option(c("--drop_0s"), type = "logical", default = user_opts$drop_0s,
                help = "Drop rows where input, upper, and lower are all 0 (default: TRUE)", metavar = "LOGICAL"),
    make_option(c("--combine_for_guide_stats"), type = "character", default = user_opts$combine_for_guide_stats,
                help = "can be '', 'sample', or 'sublib'.", metavar = "CHARACTER"),
    make_option(c("--combine_for_gene_stats"), type = "character", default = user_opts$combine_for_gene_stats,
                help = "can be 'none', 'all', 'sample', or 'sublib'", metavar = "CHARACTER"),
    
    make_option(c("--simplified_cf_threshold"), type = "integer", default = user_opts$simplified_cf_threshold,
                help = "threshold for simplified filtering", metavar = "INTEGER"),
    
    make_option(c("--strict_mode"), type = "logical", default = user_opts$strict_mode,
                help = "bolean, removes all guides which have a count of 0 in any bin.", metavar = "LOGICAL"),
    
    make_option(c("--min_guides_per_gene"), type = "integer", default = user_opts$min_guides_per_gene,
                help = "Minimum number of guides required to detect a gene, per replicate not total", metavar = "INTEGER"),
    
    make_option(c("--auto_combine_replicates"), type = "logical", default = user_opts$auto_combine_replicates,
                help = "automatically combines replicas", metavar = "LOGICAL")
  )
  
  #=============================================================================
  # Parse the command-line arguments
  #=============================================================================
  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)
  
  #=============================================================================
  # Override the defaults with CLI values
  #=============================================================================
  first_time                        <- opt$first_time
  output_folder                     <- opt$output_folder
  john_folder                       <- opt$john_folder
  john_rf_folder                    <- opt$john_rf_folder
  skip_list                         <- strsplit(opt$skip_list, ",")[[1]]
  skip_list_sublib                  <- strsplit(opt$skip_list_sublib, ",")[[1]]
  skip_list_sample                  <- strsplit(opt$skip_list_sample, ",")[[1]]
  include_controls_list             <- strsplit(opt$include_controls_list, ",")[[1]]
  
  pipeline                          <- opt$pipeline
  data_type                         <- opt$data_type
  method                            <- opt$method
  norm_method                       <- opt$norm_method
  combine_for_guide_stats           <- opt$combine_for_guide_stats
  combine_for_gene_stats            <- opt$combine_for_gene_stats
  
  recover_input                     <- opt$recover_input
  subsample_controls                <- opt$subsample_controls
  use_custom_bins                   <- opt$use_custom_bins
  same_controls_in_all_sublibraries <- opt$same_controls_in_all_sublibraries
  
  extra_suffix                      <- opt$extra_suffix
  simplified_rf_threshold           <- opt$simplified_rf_threshold
  simplified_cf_threshold           <- opt$simplified_cf_threshold
  upper_lower_percentage            <- opt$upper_lower_percentage
  
  drop_0s                           <- opt$drop_0s
  strict_mode                       <- opt$strict_mode
  min_guides_per_gene               <- opt$min_guides_per_gene
  auto_combine_replicates           <- opt$auto_combine_replicates
  
  #=============================================================================
  # Print the options (for verification)
  #=============================================================================
  cat("first_time:                        ", first_time, "\n")
  cat("output_folder:                     ", output_folder, "\n")
  cat("john_folder:                       ", john_folder, "\n")
  cat("john_rf_folder:                    ", john_rf_folder, "\n")
  
  cat("skip_list:                         ", paste(skip_list, collapse = ", "), "\n")
  cat("skip_list_sublib:                  ", paste(skip_list_sublib, collapse = ", "), "\n")
  cat("skip_list_sample:                  ", paste(skip_list_sample, collapse = ", "), "\n")
  cat("include_controls_list:             ", paste(include_controls_list, collapse = ", "), "\n")
  
  cat("pipeline:                          ", pipeline, "\n")
  cat("data_type:                         ", data_type, "\n")
  cat("method:                            ", method, "\n")
  cat("norm_method:                       ", norm_method, "\n")
  cat("combine_for_guide_stats:           ", combine_for_guide_stats, "\n")
  cat("combine_for_gene_stats:            ", combine_for_gene_stats, "\n")
  
  cat("recover_input:                     ", recover_input, "\n")
  cat("subsample_controls:                ", subsample_controls, "\n")
  cat("use_custom_bins:                   ", use_custom_bins, "\n")
  cat("same_controls_in_all_sublibraries: ", same_controls_in_all_sublibraries, "\n")
  
  cat("extra_suffix:                      ", extra_suffix, "\n")
  cat("simplified_rf_threshold:           ", simplified_rf_threshold, "\n")
  cat("simplified_cf_threshold:           ", simplified_cf_threshold, "\n")
  cat("upper_lower_percentage:            ", upper_lower_percentage, "\n")
  
  cat("drop_0s:                           ", drop_0s, "\n")
  cat("strict_mode:                       ", strict_mode, "\n")
  cat("min_guides_per_gene:               ", min_guides_per_gene, "\n")
  cat("auto_combine_replicates:           ", auto_combine_replicates, "\n")
  
  #===============================================================================
  # Construct special suffixes
  #===============================================================================
  
  if (recover_input) {
    recover_input_suffix <- "RI"
  } else {
    recover_input_suffix <- ""
  }
  if (subsample_controls) {
    subsample_controls_suffix <- "ss_cntrl"
  } else {
    subsample_controls_suffix <- ""
  }
  if (use_custom_bins) {
    custom_bins_suffix <- "custom_bins"
  } else {
    custom_bins_suffix <- ""
  }
  if (drop_0s) {
    drop_0s_suffix <- "D0"
  } else {
    drop_0s_suffix <- ""
  }
  if (strict_mode) {
    strict_mode_suffix <- "strict"
  } else {
    strict_mode_suffix <- ""
  }
  if (auto_combine_replicates){
    auto_combine_replicates_suffix <- "acr"
  } else {
    auto_combine_replicates_suffix <- ""
  }
  if (min_guides_per_gene > 0){
    min_guides_per_gene_suffix <- paste0("min_guides_",min_guides_per_gene)
  } else {
    min_guides_per_gene_suffix <- ""
  }
  if (combine_for_guide_stats == ""){
    combine_for_guide_stats_suffix <- ""
  } else {
    combine_for_guide_stats_suffix <- paste0("comb_", combine_for_guide_stats)
  }
  if (combine_for_gene_stats == ""){
    combine_for_gene_stats_suffix <- ""
  } else {
    combine_for_gene_stats_suffix <- paste0("comb_", combine_for_gene_stats)
  }
  #=============================================================================
  # Construct the skip list
  #=============================================================================
  
  skip_list_and_suffix <- create_skip_list_and_suffix(skip_list,
                                                      skip_list_sublib,
                                                      skip_list_sample)
  skip_list <- skip_list_and_suffix[[1]]
  skip_suffix <- skip_list_and_suffix[[2]]
  
  #=============================================================================
  # Construct the file suffix
  #=============================================================================
  if (use_old_suffix_construction){
    fs_parts <- c(pipeline,
                  data_type,
                  method,
                  norm_method,
                  recover_input_suffix,
                  drop_0s_suffix,
                  strict_mode_suffix,
                  min_guides_per_gene_suffix,
                  combine_for_guide_stats_suffix,
                  auto_combine_replicates_suffix,
                  skip_suffix,
                  extra_suffix)
  } else {
    fs_parts <- c(pipeline,
                  data_type,
                  method,
                  norm_method,
                  recover_input_suffix,
                  subsample_controls_suffix,
                  drop_0s_suffix,
                  strict_mode_suffix,
                  custom_bins_suffix,
                  min_guides_per_gene_suffix,
                  combine_for_guide_stats_suffix,
                  combine_for_gene_stats_suffix,
                  auto_combine_replicates_suffix,
                  skip_suffix,
                  extra_suffix)
  }

  
  # keep only non-empty parts
  fs_parts <- fs_parts[fs_parts != ""]
  
  file_suffix <- paste0("_", paste(fs_parts, collapse = "_"), ".rds")

  file_info_suffix <- paste(fs_parts, collapse = "_")
  
  # Print the final file suffix
  cat("File Suffix: ", file_suffix, "\n")
  
  #=============================================================================
  # Construct File Paths
  #=============================================================================
  
  
  
  data_dir <- get_file_path(project_root_dir,"data")
  genome_output_folder <- make_clean_dir(output_folder, "/ref/")
  dedup_output_folder <- make_clean_dir(output_folder, "/dedup/")
  mapped_output_folder <- make_clean_dir(output_folder, "/mapped/")
  rds_output_folder <- make_clean_dir(output_folder, "/rds/")
  results_output_folder <- make_clean_dir(output_folder, "/results/")
  
  merged_sgRNA_df <- readRDS(get_file_path(rds_output_folder,
                                           "merged_sgRNA_df.rds"))
  
  # return everything you want the Rmd to use
  return_list <- list(
    # =========================
    # CLI-overridable options
    # =========================
    first_time                        = first_time,
    output_folder                     = output_folder,
    john_folder                       = john_folder,
    john_rf_folder                    = john_rf_folder,
    skip_list                         = skip_list,
    skip_list_sublib                  = skip_list_sublib,
    skip_list_sample                  = skip_list_sample,
    include_controls_list             = include_controls_list,
    
    pipeline                          = pipeline,
    data_type                         = data_type,
    method                            = method,
    norm_method                       = norm_method,
    combine_for_guide_stats           = combine_for_guide_stats,
    combine_for_gene_stats            = combine_for_gene_stats,
    
    recover_input                     = recover_input,
    subsample_controls                = subsample_controls,
    use_custom_bins                   = use_custom_bins,
    same_controls_in_all_sublibraries = same_controls_in_all_sublibraries,
    
    extra_suffix                      = extra_suffix,
    simplified_rf_threshold           = simplified_rf_threshold,
    simplified_cf_threshold           = simplified_cf_threshold,
    upper_lower_percentage            = upper_lower_percentage,
    
    drop_0s                           = drop_0s,
    strict_mode                       = strict_mode,
    min_guides_per_gene               = min_guides_per_gene,
    auto_combine_replicates           = auto_combine_replicates,
    
    # =========================
    # Derived values and paths
    # =========================
    skip_suffix       = skip_suffix,
    file_suffix       = file_suffix,
    file_info_suffix  = file_info_suffix,
    
    data_dir              = data_dir,
    genome_output_folder  = genome_output_folder,
    dedup_output_folder   = dedup_output_folder,
    mapped_output_folder  = mapped_output_folder,
    rds_output_folder     = rds_output_folder,
    results_output_folder = results_output_folder,
    
    merged_sgRNA_df = merged_sgRNA_df
  )
  # Assign them as globals
  list2env(return_list, envir = envir)
  invisible(return_list)
}