# --------------------
# Import packages
# --------------------
library(dplyr)
library(tidyverse)
library(DESeq2)
library(PCAtools)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(ggplot2)
library(ggrepel)
library(scales)
library(ggVennDiagram)

# --------------------
# Input variables
# --------------------
mitocarta <- "after_mitocarta"
HTSEQ_DIR <- paste0("./data/04_htseq/mouse/liver/", mitocarta) # Input folder
OUT_DIR <- paste0("./figuras/mouse/liver/", mitocarta)

genotype <- "SKA111"

DESEQ_OUT_DIR <- paste0("./data/05_deseq/mouse/liver/", mitocarta, "/", genotype)

# SKA111
comparisons <- list(
  c("G1", "G2"),           # MALE OLD VS FEMALE OLD
  c("G5", "G9"),           # MALE SHAM VS FEMALE SHAM
  c("G6", "G5"),           # MALE ORQ VS MALE SHAM
  c("G10", "G9"),          # FEMALE OVX (1 MONTH) VS FEMALE SHAM
  c("G14", "G13")          # FEMALE OVX (7 MONTH) VS FEMALE SHAM
)

# # SKA111
# comparisons <- list(
#   c("G13", "G2"),          # FEMALE SHAM (7 MONTH) VS FEMALE OLD
#   c("G14", "G2"),          # FEMALE OVX (7 MONTH) VS FEMALE OLD
#   c("G5", "G1"),           # MALE SHAM VS MALE OLD
#   c("G6", "G1"),           # MALE ORQ VS MALE OLD
#   c("G14", "G5"),          # FEMALE OVX (7 MONTH) VS MALE SHAM
#   c("G6", "G9")            # MALE ORQ VS FEMALE SHAM
# )

# # SKA113
# comparisons <- list(
#   c("G3", "G4"),           # MALE OLD VS FEMALE OLD
#   c("G7", "G11"),          # MALE SHAM VS FEMALE SHAM
#   c("G8", "G7"),           # MALE ORQ VS MALE SHAM
#   c("G12", "G11")          # FEMALE OVX VS FEMALE SHAM
# )

# # SKA111 VS SKA113
# comparisons <- list(
#   c("G10", "G12"),         # FEMALE OVX (1 MONTH)
#   c("G5", "G7"),           # MALE SHAM
#   c("G6", "G8"),           # MALE ORQ
#   c("G9", "G11")           # FEMALE SHAM
# )

min_count <- 10                     # Minimum number of reads
smallestGroupSize <- 3

# --------------------
# Creating DESeq Dataset
# --------------------
sampleFiles <- grep("SC_L", list.files(HTSEQ_DIR), value=TRUE)

condition <- sub(".*(G[0-9]+).*", "\\1", sampleFiles)

sampleTable <- data.frame(sampleName=sub("(^.*)\\.tsv", "\\1", sampleFiles),
                          fileName=sampleFiles,
                          condition=condition)

dds <- DESeqDataSetFromHTSeqCount(sampleTable,
                                  directory=file.path(HTSEQ_DIR),
                                  design = ~ condition)

table(dds$condition)


# --------------------
# Filtering low count genes
# --------------------
keep <- rowSums(counts(dds) >= min_count) >= smallestGroupSize
dds <- dds[keep, ]

dds <- DESeq(dds)

deg_lists <- list()

for (comp in comparisons) {
  group1 <- comp[1]                 # Sample group 1 (level to compare)
  group2 <- comp[2]                 # Sample group 2 (control)

  # --------------------
  # Exploratory PCA
  # --------------------
  dds$condition <- relevel(dds$condition, ref = group2)
  dds <- nbinomWaldTest(dds)
  # resultsNames(dds)
  
  vsd <- varianceStabilizingTransformation(dds)
  vsd_subset <- vsd[, vsd$condition %in% c(group1, group2)]
  
  expr_mat <- assay(vsd_subset)
  metadata <- as.data.frame(colData(vsd_subset))
  
  pca_expl <- pca(expr_mat, metadata = metadata, removeVar = 0)
  
  dir.create(file.path(OUT_DIR, "pca", genotype), showWarnings = FALSE, recursive = TRUE)
  png(file.path(OUT_DIR, "pca", genotype, paste0("screeplot_", group1, "_vs_", group2, ".png")),
      width = 1000, height = 800, res = 120)
  
  # Calculate individual and cumulative variance explained
  var_ind <- pca_expl$variance / sum(pca_expl$variance) * 100
  cum_var <- cumsum(pca_expl$variance) / sum(pca_expl$variance) * 100
  pc <- seq_along(cum_var)
  
  # Draw blue bars and save bar centers
  bar_centers <- barplot(var_ind,
                         names.arg = pc,
                         col = "dodgerblue",
                         xlab = "Principal component",
                         ylab = "Explained variance (%)",
                         ylim = c(0, max(cum_var) * 1.1),
                         main = paste0("SCREE plot. Genotype ", genotype, ": ", group1, " vs ", group2)
                         )
  
  # Add the red cumulative variance line and points aligned with the bars
  lines(bar_centers, cum_var, type = "b", lwd = 4, col = "red", pch = 19)
  
  # Add a legend
  legend("right", legend = c("Individual variance", "Cumulative variance"),
         fill = c("dodgerblue", NA), border = NA, lty = c(NA, 1), pch = c(NA, 19), col = c("dodgerblue", "red"))
  
  dev.off()
  
  
  # Extraer datos de PCA y metadatos
  pca_data <- data.frame(
    PC1 = pca_expl$rotated[, "PC1"],
    PC2 = pca_expl$rotated[, "PC2"],
    condition = pca_expl$metadata$condition,
    sample = rownames(pca_expl$rotated)
  )
  
  # Calcular centroides por grupo (por ejemplo, ageBracket)
  centroids <- aggregate(cbind(PC1, PC2) ~ condition, data = pca_data, FUN = mean)
  
  # Distancia entre centroides (por exploración)
  dist_matrix <- dist(centroids[, c("PC1","PC2")])
  
  # Armar el biplot estilo ggplot2
  biplot_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    stat_ellipse(level = 0.95, size = 1.5) +
    geom_point(
      data = centroids,
      aes(x = PC1, y = PC2),
      shape = 4,          # Puedes elegir otro shape si lo prefieres (por ejemplo: 16 para círculo)
      size = 3,           # Tamaño más pequeño (ajústalo a gusto)
      stroke = 2,         # Grosor de línea si usas un shape con borde
      alpha = 0.6         # 0.4 significa bastante transparente; puedes ajustar entre 0 (invisible) y 1 (opaco)
    ) +
    labs(
      title = paste0("Biplot. Genotype ", genotype, ": ", group1, " vs ", group2),
      subtitle = paste0("PC1 (", round(pca_expl$variance["PC1"], 2), "%) vs PC2 (", round(pca_expl$variance["PC2"], 2), "%)"),
      x = paste0("PC1 (", round(pca_expl$variance["PC1"], 2), "% explained variance)"),
      y = paste0("PC2 (", round(pca_expl$variance["PC2"], 2), "% explained variance)"),
      color = "Condition"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 22, face = "bold"),
      plot.subtitle = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 25),
      axis.text = element_text(size = 22),
      legend.title = element_text(size = 22),
      legend.text = element_text(size = 22),
      legend.position = "right"
    )
  
  ggsave(
    filename = file.path(OUT_DIR, "pca", genotype, paste0("biplot_", group1, "_vs_", group2, ".png")),
    plot = biplot_plot,
    width = 10, height = 8, units = "in", dpi = 120
  )
  
  
  pairsplot_obj <- pairsplot(
    pca_expl,
    components = getComponents(pca_expl, 1:5),
    colby = "condition",
    pointSize = 1.5,
    axisLabSize = 12,
    title = paste0("Pairsplot. Genotype ", genotype, ": ", group1, " vs ", group2)
  )
  ggsave(
    filename = file.path(OUT_DIR, "pca", genotype, paste0("pairsplot_", group1, "_vs_", group2, ".png")),
    plot = pairsplot_obj,
    width = 15, height = 15, units = "in", dpi = 120
  )
  
  
  # --------------------
  # Differential Expression Analysis
  # --------------------
  # plotDispEsts(dds)
  
  normalized_counts_DESeq <- counts(dds, normalized=TRUE)
  write.table(normalized_counts_DESeq,
              file=paste0("./data/05_deseq/mouse/liver/", mitocarta, "/", "deseq2_normalizedCounts.tsv"),
              sep = '\t')
  
  contrast <- c("condition", group1, group2)
  
  res <- results(dds,
                 contrast = contrast,
                 alpha = 0.05)
  
  res_shrunken <- lfcShrink(dds,
                            coef = paste0("condition_", group1, "_vs_", group2),
                            type = "apeglm",
                            res = res)
  
  
  # --------------------
  # Formatting and saving deseq results
  # --------------------
  res_df <- as.data.frame(res)
  
  res_df$symbol <- rownames(res_df)
  
  gene_info <- AnnotationDbi::select(org.Mm.eg.db,
                                     keys = res_df$symbol,
                                     keytype = "SYMBOL",
                                     columns = c("SYMBOL", "GENENAME"))
  
  res_df <- res_df %>%
    left_join(gene_info, by = c("symbol" = "SYMBOL")) %>%
    mutate(FoldChange = 2^(log2FoldChange))
  
  colnames(res_df)[colnames(res_df) == "GENENAME"] <- "description"
  
  res_df <- res_df[, c("symbol", "baseMean", "FoldChange",
                       "log2FoldChange", "lfcSE", "stat",
                       "pvalue", "padj", "description")]
  
  dir.create(DESEQ_OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  write.table(res_df, file=file.path(DESEQ_OUT_DIR, paste0("deseq_all_",
                                                    group1, "_vs_", group2,
                                                    ".tsv")),
              sep="\t",
              row.names = FALSE)
  
  
  res_shrunk_df <- as.data.frame(res_shrunken)
  
  res_shrunk_df$symbol <- rownames(res_shrunk_df)
  gene_info <- AnnotationDbi::select(org.Mm.eg.db,
                                     keys = res_shrunk_df$symbol,
                                     keytype = "SYMBOL",
                                     columns = c("SYMBOL", "GENENAME"))
  
  res_shrunk_df <- res_shrunk_df %>%
    left_join(gene_info, by = c("symbol" = "SYMBOL")) %>%
    mutate(FoldChange = 2^(log2FoldChange))
  
  colnames(res_shrunk_df)[colnames(res_shrunk_df) == "GENENAME"] <- "description"
  
  res_shrunk_df <- res_shrunk_df[, c("symbol", "baseMean", "FoldChange",
                       "log2FoldChange", "lfcSE",
                       "pvalue", "padj", "description")]
  
  write.table(res_shrunk_df, file=file.path(DESEQ_OUT_DIR, paste0("deseq_all_shrunken_", group1, "_vs_", group2, ".tsv")),
              sep="\t",
              row.names = FALSE)
  
  
  
  # --------------------
  # Validation PCA
  # --------------------
  res_sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 0.379)
  sig_genes <- rownames(res_sig)
  
  if (length(sig_genes) > 1) {
    dds_sub <- dds[sig_genes, ]
    vsd_sub <- varianceStabilizingTransformation(dds_sub, blind = FALSE)
    vsd_sub_subset <- vsd_sub[, vsd_sub$condition %in% c(group1, group2)]
    
    expr_mat_sub <- assay(vsd_sub_subset)
    metadata_sub <- as.data.frame(colData(vsd_sub_subset))
    
    pca_val <- pca(expr_mat_sub, metadata = metadata_sub, removeVar = 0)
    
    png(file.path(OUT_DIR, "pca", genotype, paste0("screeplotVAL_", group1, "_vs_", group2, ".png")),
        width = 1000, height = 800, res = 120)
    # Calculate individual and cumulative variance explained
    var_ind <- pca_val$variance / sum(pca_val$variance) * 100
    cum_var <- cumsum(pca_val$variance) / sum(pca_val$variance) * 100
    pc <- seq_along(cum_var)
    
    # Draw blue bars and save bar centers
    bar_centers <- barplot(var_ind,
                           names.arg = pc,
                           col = "dodgerblue",
                           xlab = "Principal component",
                           ylab = "Explained variance (%)",
                           ylim = c(0, max(cum_var) * 1.1),
                           main = paste0("SCREE plot. Genotype ", genotype, ": ", group1, " vs ", group2)
    )
    
    # Add the red cumulative variance line and points aligned with the bars
    lines(bar_centers, cum_var, type = "b", lwd = 4, col = "red", pch = 19)
    
    # Add a legend
    legend("right", legend = c("Individual variance", "Cumulative variance"),
           fill = c("dodgerblue", NA), border = NA, lty = c(NA, 1), pch = c(NA, 19), col = c("dodgerblue", "red"))
    
    dev.off()
    
    # Comprobar cuántos PCs hay disponibles
    pcs_disp <- colnames(pca_val$rotated)
    
    if (all(c("PC1", "PC2") %in% pcs_disp)) {
      # Extraer datos de PCA y metadatos
      pca_data <- data.frame(
        PC1 = pca_val$rotated[, "PC1"],
        PC2 = pca_val$rotated[, "PC2"],
        condition = pca_val$metadata$condition,
        sample = rownames(pca_val$rotated)
      )
    
      # Calcular centroides por grupo
      centroids <- aggregate(cbind(PC1, PC2) ~ condition, data = pca_data, FUN = mean)
      
      # Distancia entre centroides (por exploración)
      dist_matrix <- dist(centroids[, c("PC1","PC2")])
      
      # Armar el biplot estilo ggplot2
      biplotVAL_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
        geom_point(size = 3) +
        stat_ellipse(level = 0.95, size = 1.5) +
        geom_point(
          data = centroids,
          aes(x = PC1, y = PC2),
          shape = 4,          # Puedes elegir otro shape si lo prefieres (por ejemplo: 16 para círculo)
          size = 3,           # Tamaño más pequeño (ajústalo a gusto)
          stroke = 2,         # Grosor de línea si usas un shape con borde
          alpha = 0.6         # 0.4 significa bastante transparente; puedes ajustar entre 0 (invisible) y 1 (opaco)
        ) +
        labs(
          title = paste0("Biplot. Genotype ", genotype, ": ", group1, " vs ", group2),
          subtitle = paste0("PC1 (", round(pca_val$variance["PC1"], 2), "%) vs PC2 (", round(pca_val$variance["PC2"], 2), "%)"),
          x = paste0("PC1 (", round(pca_val$variance["PC1"], 2), "% explained variance)"),
          y = paste0("PC2 (", round(pca_val$variance["PC2"], 2), "% explained variance)"),
          color = "Condition"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(size = 22, face = "bold"),
          plot.subtitle = element_text(size = 20, face = "bold"),
          axis.title = element_text(size = 25),
          axis.text = element_text(size = 22),
          legend.title = element_text(size = 22),
          legend.text = element_text(size = 22),
          legend.position = "right"
        )
      
      ggsave(
        filename = file.path(OUT_DIR, "pca", genotype, paste0("biplotVAL_", group1, "_vs_", group2, ".png")),
        plot = biplotVAL_plot,
        width = 10, height = 8, units = "in", dpi = 120
      )
      
      # Para pairsplot, limitar el número de componentes al máximo disponible
      n_pcs <- min(5, ncol(pca_val$rotated))
      pairsplotVAL_obj <- pairsplot(
        pca_val,
        components = getComponents(pca_val, 1:n_pcs),
        colby = "condition",
        pointSize = 1.5,
        axisLabSize = 12,
        title = paste0("Pairsplot. Genotype ", genotype, ": ", group1, " vs ", group2)
      )
      ggsave(
        filename = file.path(OUT_DIR, "pca", genotype, paste0("pairsplotVAL_", group1, "_vs_", group2, ".png")),
        plot = pairsplotVAL_obj,
        width = 15, height = 15, units = "in", dpi = 120
      )
    } else {
      message(paste0(
        "PCA de validación en ", group1, " vs ", group2,
        " tiene solo ", length(pcs_disp), " componente(s) (",
        paste(pcs_disp, collapse = ", "),
        "). Se omite el biplot/pairsplot."
      ))
  }
    } else {
    message(paste0("No significant genes for PCA validation in ", group1, " vs ", group2, ". Skipping PCA validation."))
  }
  
  
  # --------------------
  # Volcano plot
  # --------------------
  res_df_clean <- res_df[!is.na(res_df$padj), ]
  res_df_clean$threshold <- "NS"
  res_df_clean$threshold[res_df_clean$padj < 0.05 & res_df_clean$log2FoldChange > 0.379] <- "Up"
  res_df_clean$threshold[res_df_clean$padj < 0.05 & res_df_clean$log2FoldChange < -0.379] <- "Down"
  
  # Filtrar solo genes significativos (up o down)
  sig_res <- res_df_clean[
    res_df_clean$padj < 0.05 &
      abs(res_df_clean$log2FoldChange) > 0.379,
  ]
  
  # Ordenar por padj y quedarte con los 10 más significativos
  top_genes <- head(sig_res[order(sig_res$padj), ], 10)
  
  dir.create(file.path(OUT_DIR, "volcano_ma_plots", genotype), showWarnings = FALSE, recursive = TRUE)

  x_min <- min(res_df_clean$log2FoldChange, na.rm = TRUE)
  x_max <- max(res_df_clean$log2FoldChange, na.rm = TRUE)
  x_range <- x_max - x_min
  
  x_lim <- c(x_min - 0.1 * x_range, x_max + 0.1 * x_range)
    
  volcano <- ggplot(res_df_clean, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(
      values = c("Down" = "blue", "NS" = "grey", "Up" = "red"),
      name = "Significance",
      labels = c("Down" = "Downregulated",
                 "NS" = "Not significant",
                 "Up" = "Upregulated")
      ) +
    geom_vline(xintercept=c(-0.379, 0.379), color="black", linetype = "dashed", size = 0.5) +
    geom_hline(yintercept=c(-log10(0.05)), color="black", linetype = "dashed", size = 0.5) +
    annotate("text",
             x = 0.379, 
             y = min(-log10(res_df_clean$padj), na.rm=TRUE)*0.97, 
             label = "log2FC = 0.379", 
             hjust = -0.1, vjust = 1.5, size = 5) +
    annotate("text", 
             x = -0.379, 
             y = min(-log10(res_df_clean$padj), na.rm=TRUE)*0.97, 
             label = "log2FC = -0.379", 
             hjust = 1.1, vjust = 1.5, size = 5) +
    annotate("text",
             x = min(res_df_clean$log2FoldChange, na.rm=TRUE),
             y = -log10(0.05), 
             label = "padj = 0.05",
             hjust = 0, vjust = 1.5, size = 5) +
    scale_x_continuous(limits = x_lim) +
    xlab("log2 Fold Change") +
    ylab("-log10 adjusted p-value") +
    ggtitle(
      paste0(
        "Volcano Plot of Differentially Expressed Genes\n",
        "Genotype ", genotype, ": ", 
        group1, " vs ", group2
      )
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 22, face = "bold"),
      axis.title = element_text(size = 25),
      axis.text = element_text(size = 22),
      legend.title = element_text(size = 22),
      legend.text = element_text(size = 22),
      legend.position = "right"
    )
  # Añade geom_text_repel solo si hay algún significativo (o los que tú quieras)
  if (nrow(sig_res) > 0) {
    volcano <- volcano +
      geom_text_repel(
        data = top_genes, 
        aes(label = symbol), 
        size = 5, 
        max.overlaps = Inf
      )
  }
  
  ggsave(
    filename = file.path(OUT_DIR, "volcano_ma_plots", genotype, paste0("volcanoplot_", group1, "_vs_", group2, ".png")),
    plot = volcano,
    width = 12, height = 8, units = "in", dpi = 120
  )
  
  
  # --------------------
  # MA plot
  # --------------------
  png(file.path(OUT_DIR, "volcano_ma_plots", genotype, paste0("MAplot_", group1, "_vs_", group2, ".png")),
      width = 1000, height = 800, res = 120)
  
  plotMA(res, alpha = 0.05,
         main = paste0("MA plot. Genotype ", genotype, ": ", group1, " vs ", group2))
  
  dev.off()
  
  # --------------------
  # GO Classification of DEGs
  # --------------------
  # Select significant genes
  sig_genes <- res_df_clean$symbol[res_df_clean$padj < 0.05 & abs(res_df_clean$log2FoldChange) > 0.379]
  
  # Crear un nombre legible para la comparación
  comp_name <- paste0(group1, "_vs_", group2)
  
  # Lista donde acumular DEGs por comparación (defínela antes del bucle)
  # deg_lists <- list()
  
  deg_lists[[comp_name]] <- sig_genes
  
  # --------------------
  # Pie plot of DEGs
  # --------------------
  # Calcular número de DEGs y no-DEGs para la comparación
  num_deg <- length(sig_genes)
  num_total <- nrow(res_df)
  num_no_deg <- num_total - num_deg
  
  # Crear data frame para el gráfico
  pie_data <- data.frame(
    Category = c("DEG", "No DEG"),
    Count = c(num_deg, num_no_deg)
  )
  
  pie_data$Percent <- round(100 * pie_data$Count / sum(pie_data$Count), 2)
  pie_data$Label <- paste0(pie_data$Count, " (", pie_data$Percent, "%)")
  
  # Pie chart con porcentajes usando ggplot2
  pie_chart <- ggplot(pie_data, aes(x = "", y = Count, fill = Category)) +
    geom_col(color = "white") +
    coord_polar(theta = "y") +
    geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 8) +
    labs(title = paste0("Percentage of DEGs (", group1, " vs ", group2, ")")) +
    theme_void() +
    scale_fill_manual(values = c("DEG" = "#E15759", "No DEG" = "#4E79A7")) +
    theme(
      plot.title = element_text(size = 22, face = "bold"),
      legend.title = element_text(size = 22),
      legend.text = element_text(size = 22)
    )
  
  # Guardar el gráfico
  dir.create(file.path(OUT_DIR, "pie_deg", genotype), showWarnings = FALSE, recursive = TRUE)
  
  ggsave(
    filename = file.path(OUT_DIR, "pie_deg", genotype, paste0("DEG_pie_", group1, "_vs_", group2, ".png")),
    plot = pie_chart,
    width = 6, height = 6, units = "in", dpi = 120
  )
  
  if (length(sig_genes) >= 1) {
    go_ann <- AnnotationDbi::select(org.Mm.eg.db,
                                    keys = sig_genes,
                                    keytype = "SYMBOL",
                                    columns = c("GO", "ONTOLOGY"))
    
    # Add readable GO term names
    go_ann$Term <- AnnotationDbi::Term(go_ann$GO)
    
    # Remove duplicate gene-GO pairs
    go_ann_unique <- go_ann %>% distinct(SYMBOL, GO, ONTOLOGY, Term)
    
    # Count unique genes per GO term and ontology
    summary_go <- go_ann_unique %>%
      group_by(ONTOLOGY, Term) %>%
      summarise(Count = n(), .groups = "drop")
    
    # For each ontology, keep top 10 terms by Count
    exclude_terms <- c("molecular_function", "biological_process", "cellular_component")
    
    top_mf <- summary_go %>%
      filter(ONTOLOGY == "MF", !(Term %in% exclude_terms)) %>%
      arrange(desc(Count))
    
    top_bp <- summary_go %>%
      filter(ONTOLOGY == "BP", !(Term %in% exclude_terms)) %>%
      arrange(desc(Count))
    
    # top_cc <- summary_go %>%
    #   filter(ONTOLOGY == "CC", !(Term %in% exclude_terms)) %>%
    #   arrange(desc(Count))
    
    # Combine results and calculate percent of all significant genes
    top_all <- bind_rows(top_bp[1:15,], top_mf[1:15,])
    top_all$ONTOLOGY <- factor(top_all$ONTOLOGY, levels = c("BP", "MF"))
    top_all$Percent <- round(100 * top_all$Count / length(sig_genes), 1)
    
    # Order 'Term' by ontology and then descending Percent for plotting
    top_all <- top_all %>%
      arrange(ONTOLOGY, desc(Percent)) %>%
      mutate(Term = factor(Term, levels = rev(unique(Term))))
    
    # Customized grouped bar plot, colored by ontology, with count/percent labels
    go_plot <- ggplot(top_all, aes(x = Percent, y = Term, fill = ONTOLOGY)) +
      geom_bar(stat = "identity") +
      geom_text(
        aes(label = paste0(Count, " (", Percent, "%)")),
        hjust = -0.1, size = 6
      ) +
      scale_fill_manual(
        name = "GO Ontology",
        values = c("BP" = "#FB64B6", "MF" = "#21BCFF")
      ) +
      xlim(0, max(top_all$Percent)*1.15) +
      labs(
        title = paste0("GO Classification of Significant Genes (", group1, " vs ", group2, ")"),
        x = "Percentage of Genes",
        y = "GO Term"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(size = 22, face = "bold"),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 22)
      )
    
    dir.create(file.path(OUT_DIR, "go_classification_deg", genotype), showWarnings = FALSE, recursive = TRUE)
    
    ggsave(
      filename = file.path(OUT_DIR, "go_classification_deg", genotype, paste0("GO_classification_", group1, "_vs_", group2, ".png")),
      plot = go_plot,
      width = 25, height = 10, units = "in", dpi = 120
    )
  } else {
    message(paste0("No significant genes for GO classification in ", group1, " vs ", group2, ". Skipping GO classification"))
  }
}

# Filtrar comparaciones sin ningún DEG para que no molesten
deg_lists_filtered <- deg_lists[sapply(deg_lists, length) > 0]

# Por ejemplo, para 2–4 comparaciones funciona muy bien
venn_plot <- ggVennDiagram(
  deg_lists_filtered,
  label = "count",        # o "percent", "both"
  label_size = 8,         # tamaño más grande
  label_alpha = 0,
  set_size = 8
  ) +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_gradient(low = "#F4FAFF", high = "#3B7CFF") +
  theme(
    text = element_text(size = 18),
    plot.title = element_text(size = 22, face = "bold")
  ) +
  labs(title = paste0("DEG Venn diagram - ", genotype))

dir.create(file.path(OUT_DIR, "venn_deg", genotype), showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = file.path(OUT_DIR, "venn_deg", genotype, "DEG_venn.png"),
  plot = venn_plot,
  width = 12, height = 8, units = "in", dpi = 120
)

dir.create(file.path("./extra/mouse/liver/", mitocarta, genotype, "venn_deg"), showWarnings = FALSE, recursive = TRUE)

saveRDS(
  deg_lists,
  file = file.path("./extra/mouse/liver/", mitocarta, genotype, "venn_deg", "DEG_lists.rds")
)


# --------------------
# Overlapping genes
# --------------------
# deg_lists <- readRDS(paste0("./extra/mouse/liver/", mitocarta, "/", genotype, "/venn_deg/DEG_lists.rds"))
# 
# # Todas las comparaciones (nombres de las listas)
# all_comps <- names(deg_lists)
# 
# # Todos los genes que aparecen en alguna comparación
# all_genes <- sort(unique(unlist(deg_lists)))
# 
# # Construir una matriz lógica gene × comparación
# presence_mat <- sapply(deg_lists, function(gs) all_genes %in% gs)
# rownames(presence_mat) <- all_genes
# 
# presence_df <- as.data.frame(presence_mat) %>%
#   mutate(Gene = rownames(.)) %>%
#   relocate(Gene)

# common_genes <- presence_df %>%
#   filter(!G1_vs_G2 & !G5_vs_G9 & G6_vs_G5 & !G14_vs_G13) %>%
#   pull(Gene)

# common_genes <- presence_df %>%
#   filter(G1_vs_G2 & G5_vs_G9) %>%
#   pull(Gene)
# 
# res_G1 <- read.csv(file.path(DESEQ_OUT_DIR, "deseq_all_G1_vs_G2.tsv"), header=TRUE, sep='\t')
# res_G5 <- read.csv(file.path(DESEQ_OUT_DIR, "deseq_all_G5_vs_G9.tsv"), header=TRUE, sep='\t')
# 
# # common_genes es un vector de símbolos
# table_common <- data.frame(symbol = common_genes)
# 
# table_common <- table_common %>%
#   left_join(res_G1 %>%
#               dplyr::select(symbol, description = description, 
#                      log2FC_G1_vs_G2 = log2FoldChange, padj_G1_vs_G2 = padj),
#             by = "symbol") %>%
#   left_join(res_G5 %>%
#               dplyr::select(symbol, 
#                      log2FC_G5_vs_G9 = log2FoldChange, padj_G5_vs_G9 = padj),
#             by = "symbol") %>%
#   mutate(
#     Direction = case_when(
#       log2FC_G1_vs_G2 > 0 & log2FC_G5_vs_G9 > 0 ~ "up in males (both)",   # adapta al coding de tus contrastes
#       log2FC_G1_vs_G2 < 0 & log2FC_G5_vs_G9 < 0 ~ "up in females (both)",
#       TRUE ~ "opposite direction"
#     )
#   )
# 
# dir.create(file.path("./extra/mouse/liver/", mitocarta, genotype, "venn_deg", "G1_vs_G2_and_G5_vs_G9"), showWarnings = FALSE, recursive = TRUE)
# 
# write.table(table_common, file=file.path("./extra/mouse/liver/", 
#                                          mitocarta, genotype, "venn_deg", 
#                                          "G1_vs_G2_and_G5_vs_G9", 
#                                          "table_common_DEGs_G1vsG2_and_G5vsG9.tsv"),
#             sep="\t",
#             row.names = FALSE)
