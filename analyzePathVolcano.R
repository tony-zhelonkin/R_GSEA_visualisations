##################
### Custom Volcano plot refactored (V6): re-layered
##################

library(ggplot2)
library(dplyr)
library(stringr)
library(ggrepel)


analyze_pathway_volcano_custom <- function(pathway_name, gsea_results, de_results, 
                                           p_cutoff = 0.05, fc_cutoff = 2.0, 
                                           label_method = "default", max_overlaps = 100) {
  # Step 1: Extract the relevant genes for the specified pathway
  message("Extracting genes for the pathway: ", pathway_name)
  
  pathway_genes <- as.data.frame(gsea_results) %>%
    filter(Description == pathway_name) %>%
    pull(core_enrichment) %>% 
    str_split("/", simplify = TRUE) %>% 
    as.vector()
  
  # Step 2: Prepare DE results
  message("Preparing DE results...")
  
  de_results <- de_results %>%
    mutate(
      in_pathway = rownames(.) %in% pathway_genes,
      significant = adj.P.Val < p_cutoff & abs(logFC) > fc_cutoff,
      highlight = in_pathway & significant,
      color = case_when(
        highlight ~ "highlight",
        in_pathway ~ "pathway",
        TRUE ~ "grey"
      )
    )
  
  # Step 3: Determine which genes to label based on label_method
  label_data <- switch(
    label_method,
    "default" = de_results %>% filter(highlight),
    "fc" = de_results %>% filter(in_pathway & abs(logFC) > fc_cutoff),
    "p" = de_results %>% filter(in_pathway & adj.P.Val < p_cutoff),
    "all" = de_results %>% filter(in_pathway & (abs(logFC) > fc_cutoff | adj.P.Val < p_cutoff)),
    stop("Invalid label_method. Please use 'default', 'fc', 'p', or 'all'.")
  )
  
  # Step 4: Generate Volcano Plot
  message("Generating Volcano Plot...")
  
  colors <- c("grey" = "#999999", "pathway" = "#56B4E9", "highlight" = "#E69F00")
  
  # Plot the background (non-pathway genes)
  p <- ggplot(de_results, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(data = de_results %>% filter(!in_pathway), 
               aes(color = color, alpha = 0.6), size = 1) +
    scale_color_manual(values = colors) +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
    custom_minimal_theme_with_grid() +
    labs(
      title = paste("Volcano Plot -", pathway_name),
      subtitle = paste("p-value cutoff:", p_cutoff, "| FC cutoff:", fc_cutoff),
      x = "Log2 Fold Change",
      y = "-Log10(Adj. P-value)"
    ) +
    scale_x_continuous(breaks = seq(floor(min(de_results$logFC)), ceiling(max(de_results$logFC)), by = 2)) +
    theme(legend.position = "none")
  
  # Plot the pathway genes with higher alpha and size
  p <- p + geom_point(data = de_results %>% filter(in_pathway), 
                      aes(color = color, alpha = 0.7), size = 1.5) 
  
  # Adding labels for selected genes
  if (nrow(label_data) > 0) {
    p <- p + geom_label_repel(
      data = label_data,
      aes(label = rownames(label_data)),
      box.padding = 0.5,
      point.padding = 0.5,
      force = 5,
      segment.color = 'grey50',
      max.overlaps = max_overlaps
    )
  }
  
  return(p)
}
