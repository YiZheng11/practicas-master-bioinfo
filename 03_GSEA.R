# --------------------
# Import necessary packages
# --------------------
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(igraph)
require(DOSE)
library(tidyr)
library(dplyr)
library(igraph)
library(ggraph)
library(grid)

# Set organism database
organism <- "org.Hs.eg.db"
library(organism, character.only = TRUE)


# --------------------
# Input variables
# --------------------
org <- "human"
mitocarta <- "before_mitocarta"
tissue <- "heart_lv"

DESEQ_DIR <- paste0("./data/05_deseq/", org, "/", tissue, "/", mitocarta)
OUT_DIR <- paste0("./figuras/", org, "/", tissue, "/", mitocarta, "/gsea/")

genotype <- "male"

if (org == "mouse") {
  message("MOUSE")
  if (genotype == "SKA111") {
    message("SKA111 comparisons")
    # SKA111
    comparisons <- list(
      c("G1", "G2"),           # MALE OLD VS FEMALE OLD
      c("G5", "G9"),           # MALE SHAM VS FEMALE SHAM
      c("G6", "G5"),           # MALE ORQ VS MALE SHAM
      c("G10", "G9"),          # FEMALE OVX (1 MONTH) VS FEMALE SHAM
      c("G14", "G13")          # FEMALE OVX (7 MONTH) VS FEMALE SHAM
    )
  } else if (genotype == "SKA111_EXTRA") {
    message("SKA111 EXTRA comparisons")
    # SKA111
    comparisons <- list(
      c("G13", "G2"),          # FEMALE SHAM (7 MONTH) VS FEMALE OLD
      c("G14", "G2"),          # FEMALE OVX (7 MONTH) VS FEMALE OLD
      c("G5", "G1"),           # MALE SHAM VS MALE OLD
      c("G6", "G1"),           # MALE ORQ VS MALE OLD
      c("G14", "G5"),          # FEMALE OVX (7 MONTH) VS MALE SHAM
      c("G6", "G9")            # MALE ORQ VS FEMALE SHAM
    )
  } else if (genotype == "SKA113") {
    message("SKA113 comparisons")
    # SKA113
    comparisons <- list(
      c("G3", "G4"),           # MALE OLD VS FEMALE OLD
      c("G7", "G11"),          # MALE SHAM VS FEMALE SHAM
      c("G8", "G7"),           # MALE ORQ VS MALE SHAM
      c("G12", "G11")          # FEMALE OVX VS FEMALE SHAM
    )
  } else if (genotype == "SKA111_vs_SKA113") {
    message("SKA111 VS SKA113 comparisons")
    # SKA111 VS SKA113
    comparisons <- list(
      c("G10", "G12"),         # FEMALE OVX (1 MONTH)
      c("G5", "G7"),           # MALE SHAM
      c("G6", "G8"),           # MALE ORQ
      c("G9", "G11")           # FEMALE SHAM
    )
  }
} else if (org == "human") {
  message("HUMAN")
  if (genotype == "male_vs_female") {
    message("Men vs women comparisons")
    comparisons <- list(
      c("M_20_39", "F_20_39"),
      c("M_40_49", "F_40_49"),
      c("M_50_59", "F_50_59"),
      c("M_60_79", "F_60_79")
    )
  } else if (genotype == "female") {
    message("Within women comparisons")
    comparisons <- list(
      c("F_40_49", "F_20_39"),
      c("F_50_59", "F_20_39"),
      c("F_60_79", "F_20_39")
    )
  } else if (genotype  == "male") {
    message("Within men comparisons")
    comparisons <- list(
      c("M_40_49", "M_20_39"),
      c("M_50_59", "M_20_39"),
      c("M_60_79", "M_20_39")
    )
  }
}


# --------------------
# Loop over comparisons
# --------------------
for(comp in comparisons) {
  group1 <- comp[1]
  group2 <- comp[2]

  # Load DESeq results for the comparison
  df_path <- paste0(DESEQ_DIR, "/", genotype, "/deseq_all_", group1, "_vs_", group2, ".tsv")
  
  if(!file.exists(df_path)) {
    message("File not found: ", df_path)
    next
  }
  
  df <- read.csv(df_path, header=TRUE, sep='\t')
  df <- df[!is.na(df$stat), ] # remove rows with missing statistics

  # Prepare ranked gene list for GSEA
  gene_list <- df %>%
    arrange(desc(stat))    # sort according to test statistic
  
  gene_vector <- gene_list$stat
  names(gene_vector) <- gene_list$symbol # name vector with gene symbols
  
  
  # Run GSEA
  set.seed(21514)
  
  ego <- gseGO(geneList = gene_vector, 
               OrgDb = "org.Hs.eg.db", 
               keyType = "SYMBOL", 
               ont= "ALL",             # Use all GO ontologies (BP, MF, CC)
               verbose = TRUE,
               seed = TRUE)
  
  # Proceed only if we got at least one enriched term
  if (!is.null(ego) && nrow(ego@result) > 0) {
    # Calculate term-term similarity
    ego_termsim <- pairwise_termsim(ego, method = "JC", showCategory=nrow(ego@result))
    
    sim_matrix <- as.data.frame(ego_termsim@termsim)
    
    # Identify higly similar GO term pairs (> 0.8)
    high_sim <- which(sim_matrix > 0.8 & upper.tri(sim_matrix), arr.ind = TRUE)
    
    pairs_high_GO <- data.frame(
      GO1 = rownames(sim_matrix)[high_sim[,1]],
      GO2 = colnames(sim_matrix)[high_sim[,2]],
      sim = sim_matrix[high_sim]
    )
    
    # Group similar GO terms
    groups <- list()
    for (i in seq_len(nrow(pairs_high_GO))) {
      pair <- pairs_high_GO[i, ]
      group <- unique(c(pair$GO1, pair$GO2))
      groups[[i]] <- group
    }
    
    # Merge overlapping groups into non-redundant clusters
    groups <- Reduce(function(x, y) {
      merged <- FALSE
      for (i in seq_along(x)) {
        if (length(intersect(x[[i]], y)) > 0) {
          x[[i]] <- unique(c(x[[i]], y))
          merged <- TRUE
          break
        }
      }
      if (!merged) x <- c(x, list(y))
      x
    }, 
    groups, 
    list()
    )
    
    # Choose representative term per group (lowest padj)
    get_representative <- function(group, result_df) {
      sub <- result_df[result_df$Description %in% group, ]
      sub$Description[which.min(sub$p.adjust)]
    }
    
    reps <- sapply(groups, get_representative, ego@result)
    
    # Collect all core-enriched genes for each group
    get_genes_group <- function(group, ego_result) {
      genes <- ego_result$core_enrichment[ego_result$Description %in% group]
      genes_sep <- unlist(strsplit(as.character(genes), "/"))
      unique_genes <- unique(as.character(genes_sep))
      paste(unique_genes, collapse = ";")  # collapse unique genes into a single string separated by ;
    }

    genes_per_group <- sapply(groups, get_genes_group, ego_result = ego_termsim@result)

    # Build and save grouping table (representative, grouped terms, genes)
    agrupaciones <- data.frame(
      representative = reps,
      grouped_terms = sapply(groups, function(g) paste(g, collapse = ";")),
      genes = genes_per_group
    )
    
    GSEA_EXTRA_DIR <- paste0("./extra/", org, "/", tissue, "/", mitocarta, "/", genotype, "/gsea")
    
    dir.create(GSEA_EXTRA_DIR, showWarnings = FALSE, recursive = TRUE)
    
    write.table(
      agrupaciones, 
      file.path(GSEA_EXTRA_DIR, paste0("GSEA_groups_", group1, "_vs_", group2, ".tsv")), 
      sep = "\t", 
      quote = FALSE, 
      row.names = FALSE
      )

    # Build final non-redundant term set
    all_grouped_terms <- unique(unlist(groups))
    non_merged_terms <- setdiff(ego@result$Description, all_grouped_terms)
    
    final_terms_GO <- c(reps, non_merged_terms)
    
    # Subset result table to final (non‑redundant) terms
    final_res <- ego_termsim@result[ego_termsim@result$Description %in% final_terms_GO, ]
    
    # Ordered by adjusted p-value (most significant first)
    final_res <- final_res[order(final_res$p.adjust, -abs(final_res$NES)), ]
    
    # Save final non-redundant terms
    write.table(
      final_res,
      file.path(GSEA_EXTRA_DIR, paste0("GSEA_final_terms_", group1, "_vs_", group2, ".tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    # New ordered term vector
    final_terms_GO_ordered <- final_res$Description
    
    # Create output directories for plots
    dir.create(file.path(OUT_DIR, "dotplot", genotype), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(OUT_DIR, "emapplot", genotype), showWarnings = FALSE, recursive = TRUE)
    
    # Dotplot of top terms
    ego_termsim@result$sign <- ifelse(ego_termsim@result$NES > 0, "Activated", "Suppressed") # NES sign used to split activated vs suppressed
    
    dot_GO <- dotplot(ego_termsim, showCategory=final_terms_GO_ordered[1:10], split = "sign") +
      facet_grid(. ~ sign) +
      theme(
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 22),      # Texto de la leyenda
        legend.title = element_text(size = 22),     # Título de la leyenda
        strip.text = element_text(size = 22)        # Texto de facetas (Activated/Suppressed)
      )
    
    dotplot_path <- file.path(OUT_DIR, "dotplot", genotype, 
                               paste0("dotplot_", group1, "_vs_", group2, ".png"))
    
    ggsave(dotplot_path,
           plot = dot_GO, width = 20, height = 15, dpi = 120)
    
    # Enrichment map for top GO terms
    # emap_GO <- emapplot(ego_termsim, 
    #                     showCategory = final_terms_GO_ordered[1:10],
    #                     layout.params = list(layout = 'kk'), # Kamada–Kawai layout
    #                     edge.params = list(min = 0.5),       # minimum similarity for edges
    #                     cex.params = list(category_label = 1.3, category_node = 3)
    #                     ) +
    #   scale_edge_width(range = c(10,20)) +
    #   theme_graph() +
    #   theme(
    #     legend.text = element_text(size = 22),      # Texto de la leyenda
    #     legend.title = element_text(size = 22),     # Título de la leyenda
    #     )
    
    # Opcional: filtrar a los términos que quieras mostrar
    terms_to_show <- final_terms_GO_ordered[1:10]
    sim_mat_sub <- sim_matrix[terms_to_show, terms_to_show]
    
    # 1.2. Generar lista de edges con similitud > umbral (ej. 0.5)
    thresh <- 0.1
    
    edges_idx <- which(sim_mat_sub > thresh & upper.tri(sim_mat_sub), arr.ind = TRUE)
    edges_df <- data.frame(
      from       = rownames(sim_mat_sub)[edges_idx[, "row"]],
      to         = colnames(sim_mat_sub)[edges_idx[, "col"]],
      similarity = sim_mat_sub[edges_idx]
    )
    
    edges_df$show_label <- edges_df$similarity >= 0.5

    # 1.3. Crear grafo igraph
    g <- graph_from_data_frame(edges_df, directed = FALSE)
    
    total_genes_universe <- length(gene_vector)  # o el universo que definas
    res_df <- ego_termsim@result
    res_df$GeneCount <- vapply(
      res_df$core_enrichment,
      function(x) {
        genes <- unlist(strsplit(as.character(x), "/"))
        length(unique(genes))
      },
      numeric(1)
    )
    
    res_sub <- res_df[match(V(g)$name, res_df$Description), ]
    
    V(g)$GeneCount <- res_sub$GeneCount
    V(g)$p.adjust <- res_sub$p.adjust
    
    set.seed(123)
    
    g_plot <- ggraph(g, layout = "kk") +
      # Aristas: grosor ~ similarity, label separado de la línea
      geom_edge_link(
        aes(width = similarity,
            label  = ifelse(show_label, round(similarity, 2), "")),
        colour     = "grey60",
        alpha       = 0.8,
        angle_calc  = "along",                    # texto alineado con la arista
        label_dodge = unit(0, "mm"),            # separación del label respecto a la arista
        label_size = 4,
        show.legend = TRUE
      ) +
      # Nodes: tamaño por |NES|, color por p.adjust
      geom_node_point(
        aes(size = GeneCount, 
            colour = p.adjust)
        ) +
      # geom_node_text(
      #   aes(label = name),
      #   repel = TRUE,
      #   size = 5
      #   ) +
      geom_node_label(
        aes(label = name),
        repel         = TRUE,
        size          = 4,
        label.padding = grid::unit(2, "mm"),  # espacio interno del texto
        label.size    = NA,                    # sin borde
        fill          = NA,                   # sin fondo (solo texto)
        point.padding = grid::unit(0, "mm")   # separación respecto al nodo
      ) +
      scale_edge_width(
        name  = "Similarity",
        range = c(1, 4)
        ) +
      # Escalar tamaño de nodos (hazlos más grandes en general)
      scale_size_continuous(
        name = "Number of genes", 
        range = c(4, 12)
        ) +  # antes algo como c(1,6)
      scale_colour_gradient(
        name = "p.adjust", 
        low = "red",
        high = "blue",
        trans = "reverse" # opcional: rojo = más significativo
        ) +
      theme_void(base_size = 16) +
      theme(
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 12)
      )
    
    emapplot_path <- file.path(OUT_DIR, "emapplot", genotype, 
                               paste0("emapplot_", group1, "_vs_", group2, ".png"))
    
    ggsave(emapplot_path,
           plot = g_plot, width = 8, height = 8, dpi = 120)
  }
}