add_info_to_gene_stats <- function(maude_guide_stats, maude_gene_stats){
  
  maude_guide_stats <- maude_guide_stats %>%
    left_join(
      merged_sgRNA_df %>%
        mutate(entrez = as.character(entrez)) %>% 
        select(sgrna_id, entrez, seq, symbol) %>%
        distinct(sgrna_id, .keep_all = TRUE),
      by = c("sgRNA" = "sgrna_id")
    ) %>%
    mutate(entrez = coalesce(as.character(entrez), sgRNA))
  
  maude_guide_stats <- maude_guide_stats %>%
    mutate(abs_meanZ = abs(mean)) %>%
    group_by(entrez) %>%
    slice_max(order_by = abs_meanZ, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  
  maude_gene_stats <- maude_gene_stats %>%
    left_join(maude_guide_stats %>%
                select(entrez, seq, sgRNA, symbol), by = "entrez")
  
  export_df <- maude_gene_stats %>% 
    select(c(symbol, entrez, numGuides,stoufferZ,meanZ,significanceZ,p.value, FDR, seq, sgRNA)) %>% 
    arrange(significanceZ)
  # Add Liangfus expression data
  HepG2_cpm_Xie <- readRDS(file.path(data_dir, "Xie_hepato_df_cpm_RNAseq.rds"))
  
  export_df <- export_df %>%
    left_join(
      HepG2_cpm_Xie %>%
        filter(!is.na(entrez_id)) %>%
        distinct(entrez_id, .keep_all = TRUE) %>% 
        select(entrez_id,
               HepG2_CPM_Xie = HepG2,
               primary_hepato_CPM_Xie = primary),
      by = c("entrez" = "entrez_id")
    )
  HepG2_tpm_Xie <- readRDS(file.path(data_dir, "Xie_hepato_df_tpm_RNAseq.rds"))
  
  export_df <- export_df %>%
    left_join(
      HepG2_tpm_Xie %>%
        filter(!is.na(entrez_id)) %>%
        distinct(entrez_id, .keep_all = TRUE) %>% 
        select(entrez_id,
               HepG2_TPM_Xie = HepG2,
               primary_hepato_TPM_Xie = primary),
      by = c("entrez" = "entrez_id")
    )
  
  # Add Nicos HepG2 RNAseq expression data
  HepG2_tpm_Battisti <- readRDS(file.path(data_dir, "Battisti_HepG2_dual_rep_RNAseq.rds"))
  
  export_df <- export_df %>%
    left_join(
      HepG2_tpm_Battisti %>%
        filter(!is.na(entrez)) %>%
        distinct(entrez, .keep_all = TRUE) %>% 
        select(entrez,
               HepG2_CPM_Battisti = CPM_WT,
               HepG2_14D9_CPM_Battisti = CPM_dual_rep_14D9,
               HepG2_14B11_CPM_Battisti = CPM_dual_rep_14B11),
      by = "entrez"
    )
  
  # Add information about essential genes (no cardiomyocytes so far)
  essential_df <- readRDS(file.path(data_dir, "dependency_df.rds")) %>%
    mutate(entrez = as.character(entrez_id)) %>% 
    select(-c(gene_symbol,entrez_id))
  
  export_df <- export_df %>%
    left_join(essential_df, by = "entrez")
  
  cardio_cpm_Xie <- readRDS(file.path(data_dir, "Xie_cardio_df_cpm_RNAseq.rds"))
  
  export_df <- export_df %>%
    left_join(
      cardio_cpm_Xie %>%
        filter(!is.na(entrez_id)) %>%
        distinct(entrez_id, .keep_all = TRUE) %>% 
        select(entrez_id,
               cardio_IPSC_CPM = IPS_d,
               cardio_primary_CPM = primary),
      by = c("entrez" = "entrez_id")
    )
  cardio_tpm_Xie <- readRDS(file.path(data_dir, "Xie_cardio_df_tpm_RNAseq.rds"))
  
  export_df <- export_df %>%
    left_join(
      cardio_tpm_Xie %>%
        filter(!is.na(entrez_id)) %>%
        distinct(entrez_id, .keep_all = TRUE) %>% 
        select(entrez_id,
               cardio_IPSC_TPM = IPS_d,
               cardio_primary_TPM = primary),
      by = c("entrez" = "entrez_id")
    )
  cardio_cpm_Xie_Liangfu_own_data <- readRDS(file.path(data_dir, "Xie_cardio_df_cpm_Liangfus_own_RNAseq.rds"))
  
  export_df <- export_df %>%
    left_join(
      cardio_cpm_Xie_Liangfu_own_data %>%
        filter(!is.na(entrez_id)) %>%
        distinct(entrez_id, .keep_all = TRUE) %>%
        mutate(both = rowMeans(across(c(Xie_CM_Rep_1, Xie_CM_Rep_2)), na.rm = TRUE)) %>% 
        select(entrez_id,
               cardio_IPSC_CPM_Xie = both),
      by = c("entrez" = "entrez_id")
    )
  cardio_tpm_Xie_Liangfus_own_data <- readRDS(file.path(data_dir, "Xie_cardio_df_tpm_Liangfus_own_RNAseq.rds"))
  
  export_df <- export_df %>%
    left_join(
      cardio_tpm_Xie_Liangfus_own_data %>%
        filter(!is.na(entrez_id)) %>%
        distinct(entrez_id, .keep_all = TRUE) %>%
        mutate(both = rowMeans(across(c(Xie_CM_Rep_1, Xie_CM_Rep_2)), na.rm = TRUE)) %>% 
        select(entrez_id,
               cardio_IPSC_TPM_Xie = both),
      by = c("entrez" = "entrez_id")
    )
  return(export_df)
}

add_info_wrapper <- function(suffix,
                             folder = get("rds_output_folder", envir = .GlobalEnv),
                             preffix_gene = "MAUDE_gene_stats",
                             preffix_guide = "MAUDE_guide_stats"){
  gene <- readRDS(file.path(folder,paste0(preffix_gene,"_",suffix,".rds")))
  guide <- readRDS(file.path(folder,paste0(preffix_guide,"_",suffix,".rds")))
  result <- add_info_to_gene_stats(guide, gene)
  
  return(result)
}
check_UMI_counts <- function(guide_df, reference_df) {
  # Step 1: Merge reference_df's entrez column to guide_df based on matching sgRNA and sgrna_id
  guide_df <- guide_df %>%
    left_join(reference_df %>%
                select(sgrna_id, entrez), 
              by = c("sgRNA" = "sgrna_id"))
  
  # Step 2: Group by unique combinations of sublib and sample
  result_df <- guide_df %>%
    group_by(sublib, sample) %>%
    summarise(
      input_sum = sum(input, na.rm = TRUE),
      lower_sum = sum(lower, na.rm = TRUE),
      upper_sum = sum(upper, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Return the result as a dataframe
  return(result_df)
}