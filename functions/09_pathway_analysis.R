run_consensus_gsea <- function(file_info_suffix,
                               runs = 20,
                               consensus_threshold = 10,
                               go_enrichment = "BP") {
  
  results_list <- vector("list", runs)
  
  for (i in seq_len(runs)) {
    
    if (i == 1){
      results_df <- add_info_wrapper(file_info_suffix)
    } else {
      file_info_suffix_rep <- paste0(file_info_suffix,"_rep",i-1)
      results_df <- add_info_wrapper(file_info_suffix_rep)
    }
    
    geneList_Z <- results_df %>%
      filter(!is.na(significanceZ), !is.na(entrez)) %>%
      mutate(entrez = as.character(entrez)) %>%
      group_by(entrez) %>%
      mutate(score = significanceZ) %>%
      arrange(desc(score)) %>%
      { setNames(.$score, .$entrez) }
    
    gsea <- clusterProfiler::gseGO(
      geneList      = geneList_Z,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENTREZID",
      ont           = go_enrichment,
      pAdjustMethod = "BH",
      verbose       = FALSE,
      eps           = 0
    )
    
    if (!is.null(gsea) && nrow(as.data.frame(gsea)) > 0) {
      df <- as.data.frame(gsea)
      df$run <- i
      results_list[[i]] <- df
    }
    cat(paste0(
      "Completed run ",
      i,
      " of ",
      runs,
      " for consensus calling.\n"))
  }
  
  all_results <- dplyr::bind_rows(results_list)
  
  # Count occurrences
  counts <- all_results %>%
    dplyr::count(ID, name = "appearance_count")
  
  # Keep first instance per GO term
  representatives <- all_results %>%
    dplyr::group_by(ID) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  # Add counts
  all_terms <- representatives %>%
    dplyr::left_join(counts, by = "ID") %>%
    dplyr::arrange(desc(appearance_count))
  
  consensus_terms <- all_terms %>%
    dplyr::filter(appearance_count >= consensus_threshold)
  
  return(list(
    consensus_terms = consensus_terms,
    all_terms = all_terms
  ))
}

perform_pathway_analysis <- function(file_info_suffix,
                                     FDR_threshold = 0.05,
                                     go_enrichment = "BP",
                                     method = "gseGO",
                                     go_terms_qvalueCutoff = 0.05,
                                     consensus = TRUE,
                                     runs = 20, 
                                     consensus_threshold = 10,
                                     enrichGO_consensus_df = NULL){
  # This function identifies enriched pathways, either
  # using a hit list vs. the "universe" of all checked genes. 
  # This is the enrichGO method
  # OR it uses significanceZ to compare the entire screen and find enriched pathways. 
  # This is called gseGO. 
  # For consensus = FALSE; both methods simply load the results from file_info_suffix
  # For consensus = TRUE; enrichGO needs to be provided with a specific list of hits, 
  # while gseGO loads replicates named 'rep1' etc.
  
  results_df <- add_info_wrapper(file_info_suffix)
  
  if (!(go_enrichment %in% c("BP","MF","CC", "ALL"))){
    stop("go_enrichment must be 'BP', 'MF', 'ALL' or 'CC'")
  }
  if (!(method %in% c("enrichGO","gseGO"))){
    stop("method must be 'enrichGO' or 'gseGO'")
  }
  
  if (consensus){
    if (method == "gseGO"){
      
      consensus_results <- run_consensus_gsea(
        file_info_suffix,
        go_enrichment = go_enrichment,
        runs = runs,
        consensus_threshold = consensus_threshold
      )
      
      consensus_df <- consensus_results$consensus_terms
      all_df <- consensus_results$all_terms
      
      return(list(
        gsea_consensus = consensus_df,
        gsea_all = all_df
      ))
    }
    if (method == "enrichGO"){
      if (is.null(enrichGO_consensus_df)) {
        stop("Consensus calling for enrichGO requires the hit list 'enrichGO_consensus_df' to not be NULL") 
      }
      universe_genes <- results_df %>%
        pull(entrez) %>% as.character() %>% unique()
      
      sig_genes <- enrichGO_consensus_df %>%
        pull(entrez) %>%
        as.character() %>%
        unique()
      
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
      return(as.data.frame(ego))
      
    }
  } else {
    
    if (method == "enrichGO"){
      
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
      return(as.data.frame(ego))
    }
    if (method == "gseGO"){
      geneList_Z <- results_df %>%
        filter(!is.na(significanceZ), !is.na(entrez)) %>%
        mutate(entrez = as.character(entrez)) %>%
        group_by(entrez) %>%
        mutate(score = significanceZ) %>%
        arrange(desc(score)) %>%
        { setNames(.$score, .$entrez) }
      
      gsea <- clusterProfiler::gseGO(
        geneList = geneList_Z,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = go_enrichment,
        pAdjustMethod= "BH",
        verbose = FALSE,
        eps = 0 )
      return(as.data.frame(gsea))
    }

  }
}

add_STRING_id <- function(hits_df, use_local = FALSE){
  hits_df$entrez <- as.character(hits_df$entrez)
  
  string_db_path <- if (use_local) data_dir else ""

  # Create our string object
  string_db <- STRINGdb$new(
    version = "12.0",
    species = 9606,
    score_threshold = 400,
    input_directory = string_db_path
  )
  # get the string IDs for all our hits
  string_mapped <- string_db$map(
    hits_df,
    "entrez",
    removeUnmappedRows = TRUE
  ) 
  
  return(list(string_mapped = string_mapped,
              string_db = string_db))
}

get_gene_centrality_one_GO_term <- function(go_term, go_genes_df, string_db) {
  
  genes <- go_genes_df %>%
    filter(GO_term == go_term) %>%
    distinct(STRING_id, .keep_all = TRUE)
  
  # Default: no interactions
  genes$centrality <- NA
  genes$interaction_rank <- NA
  
  if (nrow(genes) >= 2) {
    
    interactions <- string_db$get_interactions(genes$STRING_id)
    
    interactions <- interactions %>%
      filter(from %in% genes$STRING_id & to %in% genes$STRING_id)
    
    if (nrow(interactions) > 0) {
      
      g <- igraph::graph_from_data_frame(
        interactions[, c("from","to")],
        directed = FALSE
      )
      
      g <- igraph::simplify(g)
      
      deg <- igraph::degree(g)
      
      centrality_df <- tibble::tibble(
        STRING_id = names(deg),
        centrality = as.numeric(deg)
      )
      
      genes <- genes %>%
        left_join(centrality_df, by = "STRING_id", suffix = c("", ".new")) %>%
        mutate(
          centrality = dplyr::coalesce(centrality.new, centrality)
        ) %>%
        select(-centrality.new)
      
      genes <- genes %>%
        mutate(
          interaction_rank = rank(-centrality, ties.method = "min")
        )
    }
  }
  return(genes)
}


get_gene_centrality <- function(file_info_suffix,
                           consensus_pathways,
                           use_local= FALSE) {
  
  full_hits <- add_info_wrapper(file_info_suffix)

  enriched_go_terms <- consensus_pathways$ID

  # add STRING_id to hits
  add_STRING_id_return_list <- add_STRING_id(full_hits, use_local = use_local)
  full_hits <- add_STRING_id_return_list$string_mapped
  string_db <- add_STRING_id_return_list$string_db
  
  # Filter for only those genes with enriched go terms
  go_genes_df <- full_hits %>%
    mutate(GO_term = GO_terms) %>% 
    separate_rows(GO_term, sep = ";") %>% 
    filter(GO_term %in% enriched_go_terms)
  
  if (nrow(go_genes_df) < 3) stop("Less than 3 enriched go terms in hits_df")
  
  core_genes <- lapply(
    enriched_go_terms,
    get_gene_centrality_one_GO_term,
    go_genes_df = go_genes_df,
    string_db = string_db
  )
  
  core_genes <- dplyr::bind_rows(core_genes)
  
  return(core_genes)
}

select_validation_genes <- function(core_genes,
                                    hits_df,
                                    consensus_pathways,
                                    top_n = 2) {

  hits_entrez <- unique(hits_df$entrez)
  
  # Keep only genes that are hits
  candidates <- core_genes %>%
    filter(entrez %in% hits_entrez)
  
  pathways <- unique(candidates$GO_term)
  
  selected <- list()
  used_genes <- c()
  
  for (pw in pathways) {
    
    pw_genes <- candidates %>%
      filter(GO_term == pw)
    
    n_specific <- sum(pw_genes$go_count <= 2, na.rm = TRUE)
    
    if (n_specific >= top_n) {
      pw_genes <- pw_genes %>%
        arrange(go_count, interaction_rank)
    } else {
      pw_genes <- pw_genes %>%
        arrange(interaction_rank, go_count)
    }
    
    pw_genes <- pw_genes %>%
      filter(!(entrez %in% used_genes))
    
    chosen <- pw_genes %>%
      slice_head(n = top_n)
    
    if (nrow(chosen) > 0) {
      
      selected[[pw]] <- chosen
      
      used_genes <- base::union(used_genes, chosen$entrez)
    }
  }
  
  selected_df <- bind_rows(selected)
  
  # Unique genes for validation
  validation_genes <- selected_df %>%
    distinct(entrez, .keep_all = TRUE)
  
  # Hits not assigned to pathways
  hits_no_pathway <- hits_df %>%
    filter(!(entrez %in% core_genes$entrez))
  
  hits_with_pathway <- hits_df %>%
    filter(entrez %in% core_genes$entrez)
  
  hits_not_selected <- hits_with_pathway %>%
    filter(!(entrez %in% selected_df$entrez))
  
  # Pathways with no selected genes
  pathways_no_hits <- setdiff(pathways, selected_df$GO_term)
  
  pathways_no_hits_df <- consensus_pathways %>%
    filter(
      ID %in% pathways_no_hits,
      !ID %in% unlist(stringr::str_extract_all(validation_genes$GO_terms, "GO:\\d+"))
    )
  
  return(list(
    validation_genes = validation_genes,
    hits_without_pathway = hits_no_pathway,
    hits_not_selected = hits_not_selected,
    pathways_without_hits = pathways_no_hits_df
  ))
}

filter_candidates_expression <- function(hits_df,
                                         expression_threshold = 0,
                                         method = "AND") {
  
  if (!(method %in% c("AND"))) {
    stop(paste(method, "is not a valid method"))
  }
  
  if (method == "AND") {
    
    remove_condition <- hits_df$HepG2_CPM_Battisti <= expression_threshold & 
      hits_df$HepG2_14B11_CPM_Battisti <= expression_threshold
    
    removed_df <- hits_df %>% 
      filter(remove_condition)
    
    kept_df <- hits_df %>% 
      filter(!remove_condition) %>% 
      mutate(
        warning = ifelse(
          (HepG2_CPM_Battisti < expression_threshold) != 
            (HepG2_14B11_CPM_Battisti < expression_threshold),
          "expression",
          NA
        )
      )
  }
  
  return(list(
    kept = kept_df,
    removed = removed_df
  ))
}

plot_significanceZ_histogram <- function(df,
                                         normalize = FALSE,
                                         fill_color = "steelblue",
                                         plot_title = NULL) {
  if (!"significanceZ" %in% colnames(df)) {
    stop("Input dataframe must contain a column named 'significanceZ'")
  }
  
  values <- df$significanceZ
  values <- values[!is.na(values)]
  
  if (length(values) == 0) {
    stop("No non-NA values found in 'significanceZ'")
  }
  
  # Custom breaks:
  # (-Inf, -10), then steps of 2 until 10, then (10, Inf)
  breaks <- c(-Inf, seq(-10, 10, by = 2), Inf)
  
  hist_obj <- hist(values, breaks = breaks, plot = FALSE, right = FALSE)
  
  plot_df <- data.frame(
    xmin = hist_obj$breaks[-length(hist_obj$breaks)],
    xmax = hist_obj$breaks[-1],
    xmid = seq_len(length(hist_obj$counts)),
    count = hist_obj$counts
  )
  
  # Labels for bins
  bin_labels <- c(
    "< -10",
    paste(seq(-10, 8, by = 2), "to", seq(-8, 10, by = 2)),
    ">= 10"
  )
  plot_df$bin_label <- bin_labels
  
  if (normalize) {
    plot_df$y <- plot_df$count / nrow(df) * 100
    plot_df$label <- paste0(round(plot_df$y, 1), "%")
    ylab_text <- "Percentage"
  } else {
    plot_df$y <- plot_df$count
    plot_df$label <- as.character(plot_df$count)
    ylab_text <- "Count"
  }
  
  if (is.null(plot_title)) {
    plot_title <- if (normalize) {
      "Histogram of significanceZ (normalized)"
    } else {
      "Histogram of significanceZ"
    }
  }
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = xmid, y = y)) +
    ggplot2::geom_col(fill = fill_color, color = "black") +
    ggplot2::geom_text(ggplot2::aes(label = label), vjust = -0.4, size = 4) +
    ggplot2::scale_x_continuous(
      breaks = plot_df$xmid,
      labels = plot_df$bin_label
    ) +
    ggplot2::labs(
      title = plot_title,
      x = "significanceZ bins",
      y = ylab_text
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::coord_cartesian(ylim = c(0, max(plot_df$y) * 1.15))
  
  return(p)
}

run_pathway_by_category <- function(df, file_info_suffix,
                                    category_col = "category",
                                    FDR_threshold = 0.05,
                                    go_enrichment = "ALL",
                                    method = "enrichGO",
                                    consensus = TRUE,
                                    go_terms_qvalueCutoff = 0.05) {
  
  # split dataframe by category
  split_dfs <- split(df, df[[category_col]])
  
  # run pathway analysis for each subset
  results <- lapply(names(split_dfs), function(cat) {
    message("Running pathway analysis for category: ", cat)
    
    perform_pathway_analysis(
      file_info_suffix = file_info_suffix,
      FDR_threshold = FDR_threshold,
      go_enrichment = go_enrichment,
      method = method,
      consensus = consensus,
      enrichGO_consensus_df = split_dfs[[cat]],
      go_terms_qvalueCutoff = go_terms_qvalueCutoff
    )
  })
  
  # add pathway analysis for all hits combined
  message("Running pathway analysis for category: All_Hits")
  results[["All_Hits"]] <- perform_pathway_analysis(
    file_info_suffix = file_info_suffix,
    FDR_threshold = FDR_threshold,
    go_enrichment = go_enrichment,
    method = method,
    consensus = consensus,
    enrichGO_consensus_df = df,
    go_terms_qvalueCutoff = go_terms_qvalueCutoff
  )
  
  names(results) <- names(split_dfs)
  return(results)
}

assign_broad_classes_and_go <- function(df,
                                         output_col = "broad_class",
                                         rules = NULL,
                                         priority = NULL) {
  
  # Add Go term information
  go_map <- suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys     = as.character(all_conjugates_overlap_df$entrez),
      keytype  = "ENTREZID",
      columns  = c("GO", "ONTOLOGY")
    )
  )
  
  go_collapsed <- go_map %>%
    mutate(ENTREZID = as.character(ENTREZID)) %>%
    group_by(ENTREZID) %>%
    summarise(
      GO_terms = paste(unique(na.omit(GO)), collapse = ";"),
      GO_ontology = paste(unique(na.omit(ONTOLOGY)), collapse = ";"),
      go_count = n_distinct(na.omit(GO)),   # annotation specificity
      .groups = "drop"
    )
  
  go_col = "GO_terms"
  
  df <- df %>%
    left_join(go_collapsed, by = c("entrez" = "ENTREZID"))
  
  # Default broad class rules
  if (is.null(rules)) {
    rules <- list(
      "Transcription factors" = c(
        "transcription", "dna-binding transcription factor", "transcription regulator"
      ),
      "RNA processing / splicing" = c(
        "rna splicing", "mrna processing", "splice", "rna binding", "ribonucleoprotein", "sirna binding", "sirna processing"
      ),
      "Chromatin regulation" = c(
        "chromatin", "histone", "nucleosome", "epigenetic", "chromosome organization"
      ),
      "Receptors" = c(
        "receptor activity", "transmembrane signaling receptor"
      ),
      "Signaling" = c(
        "signal transduction", "signaling", "kinase activity", "phosphatase activity"
      ),
      "Transporters / channels" = c(
        "transporter activity", "channel activity", "ion transport"
      ),
      "Vesicle trafficking" = c(
        "vesicle", "endocyt", "exocyt", "golgi", "lysosom", "vacuole", "recycling","endosom"
      ),
      "Cytoskeleton / adhesion" = c(
        "cytoskeleton", "microtubule", "actin", "cell adhesion", "extracellular matrix"
      ),
      "Metabolism" = c(
        "metabolic process", "biosynthetic process", "catabolic process", "mitochond"
      ),
      "Protein homeostasis" = c(
        "protein folding", "chaperone", "ubiquitin", "proteasome", "phagy"
      )
    )
  }
  
  # Default priority order
  if (is.null(priority)) {
    priority <- c(
      "Vesicle trafficking",
      "RNA processing / splicing",
      "Signaling",
      "Receptors",
      "Transporters / channels",
      "Cytoskeleton / adhesion",
      "Protein homeostasis",
      "Other / unclassified",
      "Transcription factors",
      "Chromatin regulation",
      "Metabolism"
    )
  }
  
  # Check consistency between rules and priority
  missing_from_priority <- setdiff(names(rules), priority)
  if (length(missing_from_priority) > 0) {
    stop(
      "These rule names are missing from 'priority': ",
      paste(missing_from_priority, collapse = ", ")
    )
  }
  
  # Helper: assign one class based on GO term descriptions
  assign_primary_class <- function(term_text, rules, priority) {
    if (is.na(term_text) || term_text == "") {
      return("Other / unclassified")
    }
    
    term_text_lower <- tolower(term_text)
    
    hits <- names(rules)[vapply(rules, function(patterns) {
      any(vapply(patterns, function(pat) {
        stringr::str_detect(term_text_lower, stringr::fixed(pat, ignore_case = TRUE))
      }, logical(1)))
    }, logical(1))]
    
    if (length(hits) == 0) {
      return("Other / unclassified")
    }
    
    priority[priority %in% hits][1]
  }
  
  # Extract all GO IDs from the GO_terms column
  split_go <- strsplit(as.character(df[[go_col]]), ";", fixed = TRUE)
  all_go_ids <- unique(trimws(unlist(split_go)))
  all_go_ids <- all_go_ids[!is.na(all_go_ids) & all_go_ids != ""]
  all_go_ids <- all_go_ids[grepl("^GO:\\d+$", all_go_ids)]
  
  # Build lookup from GO ID -> GO term name
  if (length(all_go_ids) > 0) {
    go_lookup <- suppressMessages(
      AnnotationDbi::select(
        GO.db,
        keys = all_go_ids,
        keytype = "GOID",
        columns = c("TERM")
      )
    )
    
    go_lookup <- go_lookup |>
      dplyr::filter(!is.na(GOID), !is.na(TERM)) |>
      dplyr::distinct(GOID, TERM)
    
    goid_to_term <- stats::setNames(go_lookup$TERM, go_lookup$GOID)
  } else {
    goid_to_term <- character(0)
  }
  
  # Expand GO IDs to GO term names for each row
  expanded_terms <- vapply(split_go, function(go_vec) {
    go_vec <- trimws(go_vec)
    go_vec <- go_vec[!is.na(go_vec) & go_vec != ""]
    go_vec <- unique(go_vec)
    
    terms <- base::unname(goid_to_term[go_vec])
    terms <- terms[!is.na(terms) & terms != ""]
    terms <- unique(terms)
    
    paste(terms, collapse = "; ")
  }, character(1))
  
  # Assign classes
  broad_class <- vapply(expanded_terms, assign_primary_class,
                        rules = rules, priority = priority,
                        character(1))
  
  # Return original df with added columns
  df$GO_terms_expanded <- expanded_terms
  df[[output_col]] <- broad_class
  
  return(df)
}


select_balanced_extremes <- function(sub_df, score_col, n_select) {
  
  if (n_select <= 0 || nrow(sub_df) == 0) {
    stop("n_select and nrow(sub_df) must both be larger than 0.")
  }
  
  if (n_select > nrow(sub_df)) {
    warning("Requested more genes than available in subset; returning all rows.")
    return(sub_df)
  }
  
  priority <- c(
    "Vesicle trafficking",
    "RNA processing / splicing",
    "Signaling",
    "Receptors",
    "Transporters / channels",
    "Cytoskeleton / adhesion",
    "Protein homeostasis",
    "Other / unclassified",
    "Transcription factors",
    "Chromatin regulation",
    "Metabolism"
  )
  
  # split equally into low and high tails, with slight preference for low/negative scores
  n_high <- floor(n_select / 2)
  n_low  <- n_select - n_high
  
  # target counts per broad_class
  class_counts <- sub_df %>%
    count(broad_class, name = "N_class") %>%
    mutate(
      prop = N_class / sum(N_class),
      n_total_target = round(prop * n_select),
      priority_rank = match(broad_class, priority),
      priority_rank = ifelse(is.na(priority_rank), length(priority) + 1, priority_rank)
    )
  
  # fix rounding so totals sum exactly to n_select
  diff_total <- n_select - sum(class_counts$n_total_target)
  
  if (diff_total != 0) {
    
    if (diff_total > 0) {
      # add to highest-priority classes first
      ord <- order(class_counts$priority_rank, decreasing = FALSE)
      
      for (i in seq_len(diff_total)) {
        idx <- ord[((i - 1) %% length(ord)) + 1]
        class_counts$n_total_target[idx] <- class_counts$n_total_target[idx] + 1
      }
      
    } else {
      # subtract from lowest-priority classes first, but only where count > 0
      ord <- order(class_counts$priority_rank, decreasing = TRUE)
      
      for (i in seq_len(abs(diff_total))) {
        possible <- ord[class_counts$n_total_target[ord] > 0]
        if (length(possible) == 0) {
          stop("Cannot reduce n_total_target further without going below zero.")
        }
        idx <- possible[1]
        class_counts$n_total_target[idx] <- class_counts$n_total_target[idx] - 1
      }
    }
  }
  
  # split per class into low/high
  class_counts <- class_counts %>%
    mutate(
      n_high_target = floor(n_total_target / 2),
      n_low_target = n_total_target - n_high_target
    )
  
  # adjust low counts so they sum to n_low
  diff_low <- n_low - sum(class_counts$n_low_target)
  
  if (diff_low != 0) {
    
    if (diff_low > 0) {
      # add low counts to highest-priority classes first
      # only possible if low can still increase by shifting from high
      ord <- order(class_counts$priority_rank, decreasing = FALSE)
      
      for (i in seq_len(diff_low)) {
        possible <- ord[class_counts$n_high_target[ord] > 0]
        if (length(possible) == 0) {
          stop("Cannot increase n_low_target further because no class has remaining n_high_target to shift from.")
        }
        idx <- possible[1]
        class_counts$n_low_target[idx] <- class_counts$n_low_target[idx] + 1
        class_counts$n_high_target[idx] <- class_counts$n_total_target[idx] - class_counts$n_low_target[idx]
      }
      
    } else {
      # subtract low counts from lowest-priority classes first
      ord <- order(class_counts$priority_rank, decreasing = TRUE)
      
      for (i in seq_len(abs(diff_low))) {
        possible <- ord[class_counts$n_low_target[ord] > 0]
        if (length(possible) == 0) {
          stop("Cannot reduce n_low_target further without going below zero.")
        }
        idx <- possible[1]
        class_counts$n_low_target[idx] <- class_counts$n_low_target[idx] - 1
        class_counts$n_high_target[idx] <- class_counts$n_total_target[idx] - class_counts$n_low_target[idx]
      }
    }
  }
  
  # low tail: smallest / most negative scores
  low_selected <- sub_df %>%
    arrange(.data[[score_col]]) %>%
    left_join(class_counts %>% select(broad_class, n_low_target), by = "broad_class") %>%
    group_by(broad_class) %>%
    group_modify(~ head(.x, .x$n_low_target[1])) %>%
    ungroup() %>%
    select(-n_low_target)
  
  # remove already selected rows before taking high tail
  remaining_df <- sub_df %>%
    anti_join(low_selected, by = "entrez")
  
  high_selected <- remaining_df %>%
    arrange(desc(.data[[score_col]])) %>%
    left_join(class_counts %>% select(broad_class, n_high_target), by = "broad_class") %>%
    group_by(broad_class) %>%
    group_modify(~ head(.x, .x$n_high_target[1])) %>%
    ungroup() %>%
    select(-n_high_target)
  
  selected <- bind_rows(low_selected, high_selected)
  
  # if rounding / class-size limits caused slight underfill, top up from remaining extremes
  if (nrow(selected) < n_select) {
    n_missing <- n_select - nrow(selected)
    
    filler <- sub_df %>%
      anti_join(selected, by = "entrez") %>%
      mutate(
        extreme_rank = pmin(
          rank(.data[[score_col]], ties.method = "first"),
          rank(-.data[[score_col]], ties.method = "first")
        )
      ) %>%
      arrange(extreme_rank) %>%
      slice_head(n = n_missing) %>%
      select(-extreme_rank)
    
    selected <- bind_rows(selected, filler)
  }
  
  return(selected)
}

simple_validation_selection <- function(df, n = 300){

  selected_df <- df %>%
    filter(!category %in% c("unique_PA", "unique_DCA", "unique_GalNAc"))
  
  n_left <- n - nrow(selected_df)
  
  if (n_left <= 0){
    stop("shared hits already larger than n")
  } else {
    cat(nrow(selected_df),"shared genes selected.\n")
  }
  
  pa_df <- df %>% filter(category == "unique_PA")
  dca_df <- df %>% filter(category == "unique_DCA")
  gal_df <- df %>% filter(category == "unique_GalNAc")
  
  n_each <- floor(n_left / 3)
  remainder <- n_left %% 3
  
  n_pa  <- n_each + ifelse(remainder >= 1, 1, 0)
  n_gal <- n_each + ifelse(remainder >= 2, 1, 0)
  n_dca <- n_each
  
  pa_sel <- select_balanced_extremes(pa_df, "significanceZ_PA", n_pa)
  dca_sel <- select_balanced_extremes(dca_df, "significanceZ_DCA", n_dca)
  gal_sel <- select_balanced_extremes(gal_df, "significanceZ_GalNAc", n_gal)
  
  final_df <- bind_rows(selected_df, pa_sel, dca_sel, gal_sel)
  
  cat(nrow(pa_sel), "PA genes selected\n")
  cat(nrow(dca_sel), "DCA genes selected\n")
  cat(nrow(gal_sel), "GalNAc genes selected\n")
  cat(nrow(final_df), "genes total\n")
  
  return(final_df)
}

calculate_category_based_Z_score <- function(df){
  df <- df %>%
    mutate(
      significanceZ = case_when(
        category == "shared_all_three" ~ rowMeans(
          pick(significanceZ_PA, significanceZ_DCA, significanceZ_GalNAc),
          na.rm = TRUE
        ),
        category == "shared_DCA_PA" ~ rowMeans(
          pick(significanceZ_DCA, significanceZ_PA),
          na.rm = TRUE
        ),
        category == "shared_DCA_GalNAc" ~ rowMeans(
          pick(significanceZ_DCA, significanceZ_GalNAc),
          na.rm = TRUE
        ),
        category == "shared_PA_GalNAc" ~ rowMeans(
          pick(significanceZ_PA, significanceZ_GalNAc),
          na.rm = TRUE
        ),
        category == "unique_DCA" ~ significanceZ_DCA,
        category == "unique_PA" ~ significanceZ_PA,
        category == "unique_GalNAc" ~ significanceZ_GalNAc,
        TRUE ~ NA_real_
      )
    )
  return(df)
}

# Validation List selection method 3: split equally into broad classes. 

select_validation_hits_balanced_classes <- function(
    df,
    n = 300,
    score_col = "significanceZ",
    category_col = "category",
    class_col = "broad_class",
    gene_id_col = "entrez",
    target_neg_fraction = 0.5,
    shared_categories = NULL,
    unique_categories = c("unique_PA", "unique_DCA", "unique_GalNAc"),
    priority = c(
      "Vesicle trafficking",
      "RNA processing / splicing",
      "Signaling",
      "Receptors",
      "Transporters / channels",
      "Cytoskeleton / adhesion",
      "Protein homeostasis",
      "Other / unclassified",
      "Transcription factors",
      "Chromatin regulation",
      "Metabolism"
    )
) {
  
  
  if (target_neg_fraction < 0 || target_neg_fraction > 1) {
    stop("target_neg_fraction must be between 0 and 1.")
  }
  target_pos_fraction <- 1 - target_neg_fraction
  
  required_cols <- c(score_col, category_col, class_col, gene_id_col)
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # if shared categories not explicitly given, define them as all non-unique categories
  if (is.null(shared_categories)) {
    shared_categories <- setdiff(unique(df[[category_col]]), unique_categories)
  }
  
  # helper: class priority rank
  get_priority_rank <- function(x) {
    out <- match(x, priority)
    out[is.na(out)] <- length(priority) + 1
    out
  }
  
  # helper: calculate exact total targets per class summing to n
  make_class_targets <- function(priority, n_total) {
    k <- length(priority)
    base_n <- floor(n_total / k)
    remainder <- n_total %% k
    
    tibble(
      !!class_col := priority,
      class_target = base_n + ifelse(seq_along(priority) <= remainder, 1, 0),
      priority_rank = seq_along(priority)
    )
  }
  
  # helper: select within one class using low/high extremes
  select_by_sign_and_extremeness <- function(sub_df, score_col, n_pos, n_neg, gene_id_col) {
    
    pos_df <- sub_df %>%
      filter(.data[[score_col]] > 0) %>%
      arrange(desc(.data[[score_col]]))
    
    neg_df <- sub_df %>%
      filter(.data[[score_col]] < 0) %>%
      arrange(.data[[score_col]])
    
    pos_sel <- pos_df %>% slice_head(n = min(n_pos, nrow(pos_df)))
    neg_sel <- neg_df %>% slice_head(n = min(n_neg, nrow(neg_df)))
    
    selected <- bind_rows(pos_sel, neg_sel) %>%
      distinct(.data[[gene_id_col]], .keep_all = TRUE)
    return(selected)
  }
  
  # helper: one candidate from a class and sign
  get_one_candidate <- function(pool_df, class_name, sign_needed, score_col) {
    sub <- pool_df %>%
      filter(.data[[class_col]] == class_name)
    
    if (sign_needed == "neg") {
      sub <- sub %>%
        filter(.data[[score_col]] < 0) %>%
        arrange(.data[[score_col]])
    } else if (sign_needed == "pos") {
      sub <- sub %>%
        filter(.data[[score_col]] > 0) %>%
        arrange(desc(.data[[score_col]]))
    } else {
      sub <- sub %>%
        mutate(
          extreme_rank = pmin(
            rank(.data[[score_col]], ties.method = "first"),
            rank(-.data[[score_col]], ties.method = "first")
          )
        ) %>%
        arrange(extreme_rank) %>%
        select(-extreme_rank)
    }
    
    slice_head(sub, n = 1)
  }
  
  # split shared vs unique
  shared_df <- df %>%
    filter(.data[[category_col]] %in% shared_categories)
  
  unique_df <- df %>%
    filter(.data[[category_col]] %in% unique_categories)
  
  if (nrow(shared_df) > n) {
    stop("Number of shared hits (", nrow(shared_df), ") is larger than n = ", n, ".")
  }
  
  message(nrow(shared_df), " shared genes selected first.")
  message(nrow(shared_df %>% filter(significanceZ > 0)),
          " positive and ",
          nrow(shared_df %>% filter(significanceZ < 0)),
          " negative genes selected")
  
  # class targets across ALL desired hits
  class_targets <- make_class_targets(priority = priority, n_total = n)
  
  # Figure out how many positive and negative hits are left to assign
  target_n_neg <- round(n * target_neg_fraction)
  target_n_pos <- n - target_n_neg
  
  shared_n_pos <- sum(shared_df[[score_col]] > 0, na.rm = TRUE)
  shared_n_neg <- sum(shared_df[[score_col]] < 0, na.rm = TRUE)
  
  remaining_pos_quota <- max(target_n_pos - shared_n_pos, 0)
  remaining_neg_quota <- max(target_n_neg - shared_n_neg, 0)
  
  message("Remaining quota after shared hits: ",
          remaining_pos_quota, " positive and ",
          remaining_neg_quota, " negative.")
  
  # how many shared hits already occupy each class
  shared_counts <- shared_df %>%
    count(.data[[class_col]], name = "n_shared") %>%
    rename(!!class_col := 1)
  
  class_targets <- class_targets %>%
    left_join(shared_counts, by = class_col) %>%
    mutate(
      n_shared = ifelse(is.na(n_shared), 0L, n_shared),
      n_to_fill = pmax(class_target - n_shared, 0L)
    )
  
  # initial fill per class from unique hits
  initial_fill_list <- list()
  
  for (i in seq_len(nrow(class_targets))) {
    class_name <- class_targets[[class_col]][i]
    n_need <- class_targets$n_to_fill[i]
    
    if (n_need <= 0) next
    
    sub_df <- unique_df %>%
      filter(.data[[class_col]] == class_name)
    
    if (nrow(sub_df) == 0) next
    
    # initial desired split for this class
    n_pos_class <- floor(n_need * target_pos_fraction)
    n_neg_class <- n_need - n_pos_class
    
    # enforce remaining global quotas
    n_pos_class <- min(n_pos_class, remaining_pos_quota)
    n_neg_class <- min(n_neg_class, remaining_neg_quota)
    
    # if class still has capacity after one sign is limited, only fill from signs that still have quota
    n_assigned_class <- n_pos_class + n_neg_class
    n_missing_class <- n_need - n_assigned_class
    
    if (n_missing_class > 0) {
      pos_avail <- sum(sub_df[[score_col]] > 0, na.rm = TRUE)
      neg_avail <- sum(sub_df[[score_col]] < 0, na.rm = TRUE)
      
      extra_neg_possible <- min(
        n_missing_class,
        remaining_neg_quota - n_neg_class,
        pos_avail * 0 + neg_avail - n_neg_class
      )
      if (extra_neg_possible > 0) {
        n_neg_class <- n_neg_class + extra_neg_possible
        n_missing_class <- n_need - (n_pos_class + n_neg_class)
      }
      
      extra_pos_possible <- min(
        n_missing_class,
        remaining_pos_quota - n_pos_class,
        pos_avail - n_pos_class
      )
      if (extra_pos_possible > 0) {
        n_pos_class <- n_pos_class + extra_pos_possible
      }
    }
    
    class_sel <- select_by_sign_and_extremeness(
      sub_df = sub_df,
      score_col = score_col,
      n_pos = n_pos_class,
      n_neg = n_neg_class,
      gene_id_col = gene_id_col
    )
    
    remaining_pos_quota <- remaining_pos_quota - sum(class_sel[[score_col]] > 0, na.rm = TRUE)
    remaining_neg_quota <- remaining_neg_quota - sum(class_sel[[score_col]] < 0, na.rm = TRUE)
    
    initial_fill_list[[length(initial_fill_list) + 1]] <- class_sel
  }
  
  initial_fill <- bind_rows(initial_fill_list)
  
  selected_df <- bind_rows(shared_df, initial_fill) %>%
    distinct(.data[[gene_id_col]], .keep_all = TRUE)
  
  # remaining unique pool after initial fill
  remaining_pool <- unique_df %>%
    anti_join(selected_df %>% select(all_of(gene_id_col)), by = gene_id_col)
  
  n_left <- n - nrow(selected_df)
  
  message(nrow(initial_fill), " genes selected in class-balancing step.")
  message(nrow(initial_fill %>% filter(significanceZ > 0)),
          " positive and ",
          nrow(initial_fill %>% filter(significanceZ < 0)),
          " negative genes selected")
  message(n_left,
          " slots remain after class targets were filled as much as possible.")
  
  # iterative fill based on current sign imbalance
  while (nrow(selected_df) < n && nrow(remaining_pool) > 0) {
    
    n_pos <- sum(selected_df[[score_col]] > 0, na.rm = TRUE)
    n_neg <- sum(selected_df[[score_col]] < 0, na.rm = TRUE)
    
    remaining_pos_quota <- target_n_pos - n_pos
    remaining_neg_quota <- target_n_neg - n_neg
    
    if (remaining_neg_quota > 0 && remaining_pos_quota <= 0) {
      sign_needed <- "neg"
    } else if (remaining_pos_quota > 0 && remaining_neg_quota <= 0) {
      sign_needed <- "pos"
    } else if (remaining_neg_quota > remaining_pos_quota) {
      sign_needed <- "neg"
    } else if (remaining_pos_quota > remaining_neg_quota) {
      sign_needed <- "pos"
    } else {
      sign_needed <- "either"
    }
    
    added_this_round <- FALSE
    
    for (class_name in priority) {
      if (nrow(selected_df) >= n) break
      
      candidate <- get_one_candidate(
        pool_df = remaining_pool,
        class_name = class_name,
        sign_needed = sign_needed,
        score_col = score_col
      )
      
      if (nrow(candidate) == 0) {
        next
      }
      
      if (nrow(candidate) > 0) {
        selected_df <- bind_rows(selected_df, candidate)
        remaining_pool <- remaining_pool %>%
          anti_join(candidate %>% select(all_of(gene_id_col)), by = gene_id_col)
        added_this_round <- TRUE
      }
    }
    
    # if no class contributed in this round, take most extreme remaining hits overall
    if (!added_this_round && nrow(remaining_pool) > 0) {
      filler <- remaining_pool %>%
        mutate(
          extreme_rank = pmin(
            rank(.data[[score_col]], ties.method = "first"),
            rank(-.data[[score_col]], ties.method = "first")
          ),
          priority_rank = get_priority_rank(.data[[class_col]])
        ) %>%
        arrange(priority_rank, extreme_rank) %>%
        slice_head(n = min(n - nrow(selected_df), nrow(.))) %>%
        select(-extreme_rank, -priority_rank)
      
      selected_df <- bind_rows(selected_df, filler)
      remaining_pool <- remaining_pool %>%
        anti_join(filler %>% select(all_of(gene_id_col)), by = gene_id_col)
    }
  }
  
  # summaries
  final_class_counts <- selected_df %>%
    count(.data[[class_col]], name = "n_selected") %>%
    rename(!!class_col := 1) %>%
    right_join(class_targets %>% select(all_of(class_col), class_target, n_shared), by = class_col) %>%
    mutate(n_selected = ifelse(is.na(n_selected), 0L, n_selected))
  
  n_pos_final <- sum(selected_df[[score_col]] > 0, na.rm = TRUE)
  n_neg_final <- sum(selected_df[[score_col]] < 0, na.rm = TRUE)
  
  message(nrow(selected_df), " genes total selected.")
  message(n_pos_final, " positive and ", n_neg_final, " negative genes in final set.")
  
  return(list(
    selected_df = selected_df,
    class_summary = final_class_counts,
    n_positive = n_pos_final,
    n_negative = n_neg_final,
    n_shared_selected = nrow(shared_df)
  ))
}
plot_class_distribution <- function(
    df_all,
    df_validation,
    sign_split = FALSE
) {

  class_levels <- rev(c(
    "Transcription factors",
    "RNA processing / splicing",
    "Chromatin regulation",
    "Receptors",
    "Signaling",
    "Transporters / channels",
    "Vesicle trafficking",
    "Cytoskeleton / adhesion",
    "Metabolism",
    "Protein homeostasis",
    "Other / unclassified"
  ))
  
  if (!sign_split) {
    
    plot_df <- bind_rows(
      df_all %>%
        count(broad_class, name = "n") %>%
        mutate(
          group = "All",
          percentage = n / sum(n) * 100
        ),
      df_validation %>%
        count(broad_class, name = "n") %>%
        mutate(
          group = "Validation",
          percentage = n / sum(n) * 100
        )
    ) %>%
      select(broad_class, group, percentage) %>%
      mutate(
        broad_class = factor(broad_class, levels = class_levels),
        group = factor(group, levels = c("All", "Validation"))
      )
    
    p <- ggplot(plot_df, aes(x = percentage, y = broad_class, fill = group)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      geom_text(
        aes(label = paste0(round(percentage, 1), "%")),
        position = position_dodge(width = 0.8),
        hjust = -0.1,
        size = 3.5
      ) +
      scale_fill_manual(values = c("All" = "forestgreen", "Validation" = "steelblue")) +
      scale_x_continuous(
        limits = c(0, max(plot_df$percentage) * 1.15),
        labels = label_number(suffix = "%")
      ) +
      labs(
        x = "Percentage of genes",
        y = NULL,
        fill = NULL
      ) +
      theme_classic()
    
  } else {
    
    plot_df <- bind_rows(
      df_all %>%
        mutate(sign_group = ifelse(significanceZ > 0, "Pos", "Neg")) %>%
        count(broad_class, sign_group, name = "n") %>%
        mutate(
          percentage = n / nrow(df_all) * 100,
          group = paste0("All_", sign_group)
        ),
      
      df_validation %>%
        mutate(sign_group = ifelse(significanceZ > 0, "Pos", "Neg")) %>%
        count(broad_class, sign_group, name = "n") %>%
        mutate(
          percentage = n / nrow(df_validation) * 100,
          group = paste0("Validation_", sign_group)
        )
    ) %>%
      select(broad_class, group, percentage) %>%
      mutate(
        broad_class = factor(broad_class, levels = class_levels),
        group = factor(
          group,
          levels = c("All_Pos", "All_Neg", "Validation_Pos", "Validation_Neg")
        )
      )
    
    p <- ggplot(plot_df, aes(x = percentage, y = broad_class, fill = group)) +
      geom_col(position = position_dodge(width = 0.85), width = 0.7) +
      geom_text(
        aes(label = paste0(round(percentage, 1), "%")),
        position = position_dodge(width = 0.85),
        hjust = -0.1,
        size = 2
      ) +
      scale_fill_manual(values = c(
        "All_Pos" = "forestgreen",
        "All_Neg" = "darkgreen",
        "Validation_Pos" = "steelblue",
        "Validation_Neg" = "navy"
      )) +
      scale_x_continuous(
        limits = c(0, max(plot_df$percentage) * 1.15),
        labels = label_number(suffix = "%")
      ) +
      labs(
        x = "Percentage of genes",
        y = NULL,
        fill = NULL
      ) +
      theme_classic()
  }
  
  return(p)
}
