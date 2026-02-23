#===============================================================================
# Imports and loads
#===============================================================================
library(fs)

#===============================================================================
# Setup block functions
#===============================================================================
# Extends the skip_list to skip entire sublibraries or sample numbers.
create_skip_list_and_suffix <- function(skip_list, skip_list_sublib, skip_list_sample){
  skip_suffix <- ""
  
  bins <- c("I","L","U")
  
  # helper to safely inject literal values into regex (important if A could contain regex chars)
  esc_re <- function(x) gsub("([][{}()+*^$|\\\\.?])", "\\\\\\1", x, perl = TRUE)
  
  if (length(skip_list_sublib) > 0) {
    # regex: [ILU]_A_[^_]+   (A from skip_list_sublib)
    pats <- as.vector(outer(bins, esc_re(skip_list_sublib),
                            FUN = function(b, a) paste0("^", b, "_", a, "_[^_]+$")))
    skip_list <- c(skip_list, paste0("re:", pats))
    skip_suffix <- paste0(skip_suffix, "no_sublib_", paste0(skip_list_sublib, collapse = "_"))
  }
  
  if (length(skip_list_sample) > 0) {
    if (length(skip_list_sublib) > 0) skip_suffix <- paste0(skip_suffix, "_")
    
    # regex: [ILU]_L\\d+_A   (A from skip_list_sample)
    # (anchored; allows L + digits, any width)
    pats <- as.vector(outer(bins, esc_re(skip_list_sample),
                            FUN = function(b, a) paste0("^", b, "_L\\d+_", a, "$")))
    skip_list <- c(skip_list, paste0("re:", pats))
    skip_suffix <- paste0(skip_suffix, "no_sample_", paste0(skip_list_sample, collapse = "_"))
  }
  
  skip_list <- unique(skip_list)
  return(list(skip_list, skip_suffix))
}

# Create a directory
make_clean_dir <- function(base_path, sub_path = NULL) {
  full_path <- if (is.null(sub_path) || sub_path == "") base_path else fs::path(base_path, sub_path)
  fs::dir_create(full_path, recurse = TRUE)
  fs::path_norm(full_path)
}

# construct a file path
get_file_path <- function(folder_path, file_name) {
  fs::path_norm(fs::path(folder_path, file_name))
}
# Old homebrew versions of the two

# make_clean_dir <- function(base_path, sub_path) {
#   
#   full_path <- paste0(base_path, "/", sub_path)
#   
#   full_path <- gsub("///", "/", full_path, fixed = TRUE)
#   full_path <- gsub("//", "/", full_path, fixed = TRUE)
#   
#   if (!dir.exists(full_path)) {
#     dir.create(full_path, recursive = TRUE)
#   }
#   return(full_path)
# }
# get_file_path <- function(folder_path, file_name){
#   full_path <- paste0(folder_path, "/",file_name)
#   full_path <- gsub("///", "/", full_path, fixed = TRUE)
#   full_path <- gsub("//", "/", full_path, fixed = TRUE)
#   return(full_path)
# }