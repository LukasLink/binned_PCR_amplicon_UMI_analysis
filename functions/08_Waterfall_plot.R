plot_significance_by_rank <- function(Hits_df,
                                      mark_cntrl = TRUE,
                                      mark_special = NULL,
                                      mark_N_top_hits = 0,
                                      box_padding = 0.8,
                                      no_text = FALSE,
                                      signif_lines = FALSE,
                                      mark_all_signif_level = NULL) {
  # Create rank column
  Hits_df <- Hits_df %>%
    mutate(rank = rank(significanceZ, ties.method = "first"))
  
  # Define the Limits for my two label boxes
  y_min   <- min(Hits_df$significanceZ, na.rm = TRUE)
  y_max   <- max(Hits_df$significanceZ, na.rm = TRUE)
  y_range <- y_max - y_min
  
  # Box geometry (relative to overall range)
  box_height <- 0.06 * y_range   # 8% of range high
  box_gap    <- 0.02 * y_range   # 3% of range between boxes
  
  # Centers of the two boxes (on y-axis)
  y_center1 <- y_min           # lower box
  y_center2 <- y_center1 + box_height + box_gap    # upper box
  
  # x range (can reuse for both boxes)
  x_min_box <- max(Hits_df$rank) - 3000
  x_max_box <- max(Hits_df$rank) + 500
  
  # Start base plot
  p <- ggplot(Hits_df, aes(x = rank, y = significanceZ)) +
    geom_point(color = "lightgrey", size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
    theme_bw() +
    labs(
      x = "Gene rank",
      y = "Significance Z-score",
      title = "Significance Z-score by gene rank"
    )
  
  # Mark control (NA in symbol column)
  if (mark_cntrl && "symbol" %in% names(Hits_df)) {
    p <- p + geom_point(
      data = Hits_df %>% filter(is.na(symbol)),
      aes(x = rank, y = significanceZ),
      color = "blue", size = 1.5
    ) +
      annotate("rect",
               xmin = x_min_box,
               xmax = x_max_box,
               ymin = y_center2 - box_height / 2,
               ymax = y_center2 + box_height / 2,
               fill = "white", color = "black", size = 0.3) +
      annotate("point",
               x = x_min_box + 200,
               y = y_center2,
               color = "blue", size = 2) +
      annotate("text",
               x = x_min_box + 400,
               y = y_center2,
               label = "Control sgRNA", hjust = 0, size = 3)
  }
  
  if (!is.null(mark_all_signif_level)) {
    signif_df <- Hits_df %>%
      filter(FDR <= mark_all_signif_level)
    signif_label_text <- paste0("Hits: FDR <= ",mark_all_signif_level*100,"%")
    
    if (nrow(signif_df) > 0) {
      # Highlight significant points
      p <- p +
        geom_point(
          data = signif_df,
          aes(x = rank, y = significanceZ),
          color = "red",
          size = 1.5
        ) +
        annotate("rect",
                 xmin = x_min_box,
                 xmax = x_max_box,
                 ymin = y_center1 - box_height / 2,
                 ymax = y_center1 + box_height / 2,
                 fill = "white", color = "black", size = 0.3) +
        annotate("point",
                 x = x_min_box + 200,
                 y = y_center1,
                 color = "red", size = 2) +
        annotate("text",
                 x = x_min_box + 400,
                 y = y_center1,
                 label = signif_label_text, hjust = 0, size = 3)
    }
    
    # Add lines if requested
    if (isTRUE(signif_lines)) {
      
      ## Lowest positive significanceZ (closest to zero, > 0)
      pos_df <- signif_df %>% filter(significanceZ > 0)
      if (nrow(pos_df) > 0) {
        y_pos <- min(pos_df$significanceZ, na.rm = TRUE)
        p <- p +
          geom_hline(
            yintercept = y_pos,
            linetype = "dashed",
            color = "red",
            size = 0.3
          )
      }
      
      ## Negative significanceZ with smallest absolute value (closest to zero)
      neg_df <- signif_df %>% filter(significanceZ < 0)
      if (nrow(neg_df) > 0) {
        # largest negative value, e.g. -1 is "closer to zero" than -3
        y_neg <- max(neg_df$significanceZ, na.rm = TRUE)
        p <- p +
          geom_hline(
            yintercept = y_neg,
            linetype = "dashed",
            color = "red",
            size = 0.3
          )
      }
    }
  }
  
  
  # Mark top hits (both ends), excluding special symbols
  if (mark_N_top_hits >= 1 && "symbol" %in% names(Hits_df)) {
    top_hits_df <- Hits_df %>%
      filter(!(symbol %in% mark_special)) %>%  # exclude special symbol
      arrange(significanceZ) %>%
      slice(c(1:mark_N_top_hits, (n() - mark_N_top_hits + 1):n()))
    
    top_hits_pos <- top_hits_df %>% filter(significanceZ >= 0)
    top_hits_neg <- top_hits_df %>% filter(significanceZ < 0)
    
    if (no_text == TRUE){
      p <- p +
        geom_point(data = top_hits_df, aes(x = rank, y = significanceZ),
                   color = "red", size = 1.5)
    } else {
      p <- p +
        geom_point(data = top_hits_df, aes(x = rank, y = significanceZ),
                   color = "red", size = 1.5) +
        geom_text_repel(
          data = top_hits_df,
          aes(x = rank, y = significanceZ, label = symbol),
          color = "red",           # or "orange" for special symbol
          size = 3,
          box.padding = box_padding,       # space around the label box
          point.padding = 0.6,     # space around the point
          force = 10,               # strength of repulsion between labels
          force_pull = 0.5,          # how strongly labels are pulled back to points
          max.overlaps = Inf,      # allow all overlaps to be considered
          max.time = 5,           # more time for optimization
          min.segment.length = 0,
          segment.color = "red",# line color from label to point
          segment.size = 0.3       # line thickness
        )
    }
  }
  # Mark special symbols (can be multiple)
  if (!is.null(mark_special) && length(mark_special) > 0 && "symbol" %in% names(Hits_df)) {
    special_df <- Hits_df %>% filter(symbol %in% mark_special)
    
    if (nrow(special_df) > 0) {
      p <- p +
        geom_point(
          data = special_df,
          aes(x = rank, y = significanceZ),
          shape = 17,
          color = "orange",
          size = 4
        ) +
        geom_text_repel(
          data = special_df,
          aes(x = rank, y = significanceZ, label = symbol),
          color = "orange",
          size = 4,
          box.padding = 1,
          nudge_x = ifelse(special_df$significanceZ > 0, -500, 500),
          nudge_y = ifelse(special_df$significanceZ > 0, y_range * (0.05), y_range * (-0.05)),
          force = 4
        )
    }
  }
  return(p)
}