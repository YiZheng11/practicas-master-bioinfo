# --------------------
# Import necessary packages
# --------------------
library(ggplot2)
library(data.table)
library(dplyr)
library(pheatmap)

# --------------------
# Input variables
# --------------------
mitocarta <- "before_mitocarta"

DESEQ_DIR <- paste0("./data/05_deseq/mouse/liver/", mitocarta)
GENES_DIR <- "./data/genes_of_interest"
OUTPUT_DIR <- paste0("./figuras/mouse/liver/", mitocarta, "/log2fc_goi")

genotype <- "SKA111_EXTRA"

DESEQ_OUT_DIR <- paste0("./data/05_deseq/mouse/liver/", mitocarta, "/", genotype)

# # SKA111
# comparisons <- list(
#   c("G1", "G2"),           # MALE OLD VS FEMALE OLD
#   c("G5", "G9"),           # MALE SHAM VS FEMALE SHAM
#   c("G6", "G5"),           # MALE ORQ VS MALE SHAM
#   c("G10", "G9"),          # FEMALE OVX (1 MONTH) VS FEMALE SHAM
#   c("G14", "G13")          # FEMALE OVX (7 MONTH) VS FEMALE SHAM
# )

# SKA111
comparisons <- list(
  c("G13", "G2"),          # FEMALE SHAM (7 MONTH) VS FEMALE OLD
  c("G14", "G2"),          # FEMALE OVX (7 MONTH) VS FEMALE OLD
  c("G5", "G1"),           # MALE SHAM VS MALE OLD
  c("G6", "G1"),           # MALE ORQ VS MALE OLD
  c("G14", "G5"),          # FEMALE OVX (7 MONTH) VS MALE SHAM
  c("G6", "G9")            # MALE ORQ VS FEMALE SHAM
)

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


dir.create(file.path(OUTPUT_DIR, genotype), recursive = TRUE, showWarnings = FALSE)

# Get all genes of interest files
gene_files <- list.files(GENES_DIR, pattern = "\\.txt$", full.names = TRUE)

# Loop through each comparison group
for (comp in comparisons) {
  group1 <- comp[1]
  group2 <- comp[2]
  
  # Load DESeq2 results for this comparison
  deseq_path <- paste0(DESEQ_DIR, "/", genotype, "/deseq_all_", group1, "_vs_", group2, ".tsv")
  if (!file.exists(deseq_path)) {
    message("DESeq2 results not found for ", group1, " vs ", group2, ". Skipping.")
    next
  }
  
  df <- read.csv(deseq_path, header=TRUE, sep='\t')

  # Select columns with gene name, log2FC and padj values
  original_data <- data.frame(
    Gene = df$symbol,
    log2FC = df$log2FoldChange,
    padj = df$padj
  )

  # Loop over each genes of interest set
  for (gene_file in gene_files) {
    GOI_name <- tools::file_path_sans_ext(basename(gene_file))    # get basename of file (without extension)
    message("Processing: ", GOI_name, " for comparison ", group1, " vs ", group2)
    
    # Load GOI list
    genes_of_interest <- read.csv(gene_file, header = FALSE, sep = '\t')
    
    # Filter for genes present in current DESeq results
    data <- original_data[original_data$Gene %in% genes_of_interest[[1]],]
    data <- data %>%
      filter(!is.na(padj))   # Remove rows with padj = NA
    
    # Add regulation status
    data_filtered <- data %>%
      mutate(Regulation = ifelse(log2FC >= 0, "UP", "DOWN")) %>%
      mutate(
        sig_label = case_when(
          padj < 0.001      ~ "***",
          padj < 0.01       ~ "**",
          padj < 0.05       ~ "*",
          TRUE              ~ ""
          )
        )
    
    # If no valid genes, skip this plot
    if (nrow(data_filtered) == 0) {
      message("No genes of interest found in: ", GOI_name, " for ", group1, " vs ", group2)
      next
    }
    
    # Calculate recommended font size (larger for <20 genes, smaller otherwise)
    gene_num <- nrow(data_filtered)
    gene_font_size <- ifelse(gene_num <= 35, 20, 10)
    
    # Create and save log2FC barplot for these genes
    p <- ggplot(data_filtered, aes(Gene, log2FC, fill=Regulation)) +
      geom_bar(stat='identity') +
      geom_hline(yintercept=c(-0.38, 0.38), color="black", linetype = "dashed", size = 0.5) +
      geom_hline(yintercept = 0, color="black", size = 0.75) +
      geom_text(aes(label = sig_label),
                hjust = ifelse(data_filtered$log2FC > 0, -0.3, 1.3),
                vjust = 0.5,
                size = 5) +
      scale_fill_manual(name = "Regulation",
                        values = c("UP" = "red",
                                   "DOWN" = "blue")) +
      scale_y_continuous(
        breaks = scales::breaks_pretty(n = 8),
        limits = c(
          floor(min(data_filtered$log2FC) - 0.5),
          ceiling(max(data_filtered$log2FC) + 0.5)
        )
      ) +
      coord_flip() +
      theme_minimal(base_size = 14) + 
      theme(
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = gene_font_size),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 22)
      )
      ggtitle(paste0(GOI_name, ". Genotype ", genotype, ": ", group1, " vs ", group2))
    
    output_path <- file.path(OUTPUT_DIR, genotype,
                             paste0(group1, "_vs_", group2, "_", GOI_name, ".png"))
    
    ggsave(output_path, plot = p, width = 10, height = 15, dpi = 300)
    
    message("Log2FC figure saved in: ", output_path)
  }
}