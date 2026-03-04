
perform_pathway_analysis <- function(results_df, FDR_threshold = 0.05, go_enrichment = "BP", go_terms_qvalueCutoff = 0.05){
  
  if (!(go_enrichment %in% c("BP","MF","CC"))){
    stop("go_enrichment must be 'BP', 'MF', or 'CC'")
  }
  
  sig_genes <- results_df %>%
    filter(!is.na(FDR)) %>%
    filter(FDR < FDR_threshold) %>%
    pull(entrez) %>%
    as.character() %>%
    unique()
  
  universe_genes <- results_df %>%
    pull(entrez) %>% as.character() %>% unique()
  
  ego <- clusterProfiler::enrichGO(
    gene          = sig_genes,
    universe      = universe_genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = go_enrichment,
    pAdjustMethod = "BH",
    qvalueCutoff  = go_terms_qvalueCutoff,
    readable      = TRUE
  )
  
  geneList_Z <- results_df %>%
    filter(!is.na(significanceZ), !is.na(entrez)) %>%
    mutate(entrez = as.character(entrez)) %>%
    group_by(entrez) %>%
    mutate(score = significanceZ) %>%
    arrange(desc(score)) %>%
    { setNames(.$score, .$entrez) }
  
  gsea <- clusterProfiler::gseGO(
    geneList     = geneList_Z,
    OrgDb        = org.Hs.eg.db,
    keyType      = "ENTREZID",
    ont          = go_enrichment,
    pAdjustMethod= "BH",
    verbose      = FALSE
  )
  return(list(
    gsea = gsea,
    ego = ego
  ))
}