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

group1 <- "G1_vs_G2"
group2 <- "G5_vs_G9"
sig_genes <- common_genes

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
    title = paste0("GO Classification of Common Genes between ", group1, " and ", group2, ")"),
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

dir.create(file.path("./extra/mouse/liver/", "after_mitocarta", "SKA111", "venn_deg", paste0(group1, "_and_", group2)), showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = file.path("./extra/mouse/liver/", "after_mitocarta", "SKA111", "venn_deg", paste0(group1, "_and_", group2), paste0("GO_classification_", group1, "_and_", group2, ".png")),
  plot = go_plot,
  width = 25, height = 10, units = "in", dpi = 120
)
