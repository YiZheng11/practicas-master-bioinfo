# --------------------
# Load libraries
# --------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(dunn.test)
library(car)
library(rcompanion)


# --------------------
# Input variables
# --------------------
mitocarta <- "after_mitocarta"
org <- "human"
tissue <- "liver"

OUT_DIR <- paste0("./figuras/", org, "/", tissue, "/", mitocarta, "/read_distribution")

if (org == 'mouse') {
  HTSEQ_DIR <- paste0("./data/04_htseq/", org, "/", tissue, "/", mitocarta)
  files <- list.files(path=HTSEQ_DIR, 
                      pattern="\\.tsv$", 
                      full.names=TRUE)
} else {
  HTSEQ_DIR <- paste0("./data/04_htseq/", org, "/", mitocarta)
  files <- list.files(path = HTSEQ_DIR)
}

if (tissue == 'liver') {
  tissue_n <- 3
} else if (tissue == 'heart_aa') {
  tissue_n <- 1
} else {
  tissue_n <- 2
}


# --------------------
# Function to read htseq file
# --------------------
read_htseq_mouse <- function(file) {
  # Extract sample name from file name
  sample_name <- sub(".*(G[0-9]+.*).tsv", "\\1", basename(file))
  
  # Read raw counts
  df <- read.table(file, 
                   header = FALSE, 
                   col.names = c("gene", "count"))
  
  # Remove htseq summary rows (technical counters, not genes)
  df <- df %>%
    filter(!gene %in% c("__no_feature", "__ambiguous", 
                        "__too_low_aQual", "__not_aligned", 
                        "__alignment_not_unique"))
  
  # Add sample name as a column
  df$sample <- sample_name
  
  return(df)
}

read_htseq_human <- function(file) {
  if (mitocarta == "before_mitocarta") {
    message(paste0("Reading whole GTEx human ", tissue, " file."))
    df <- read.delim(
      file,
      skip = 2,
      header = TRUE,
      sep = "\t",
      check.names = FALSE,
      row.names = 1
    )
    
    df <- df[, -c(1)]
    
    df_long <- df %>%
      as.data.frame() %>%
      tibble::rownames_to_column("gene") %>%   # gene = Ensembl ID
      pivot_longer(
        cols = -gene,
        names_to = "sample",
        values_to = "count"
      )
  } else {
    message(paste0("Reading filtered by MitoCarta GTEx human ", tissue, " file."))
    df <- read.delim(
      file,
      header = TRUE,
      sep = "\t",
      check.names = FALSE,
      row.names = 3
    )
    
    df <- df[, -c(1, 2, 3)]
    
    df_long <- df %>%
      as.data.frame() %>%
      tibble::rownames_to_column("gene") %>%   # gene = Ensembl ID
      pivot_longer(
        cols = -gene,
        names_to = "sample",
        values_to = "count"
      )
  }
}


# --------------------
# Read and combine all htseq files
# --------------------
if (org == 'mouse') {
  counts <- lapply(files, read_htseq) %>% bind_rows()
} else if (org == 'human') {
  counts <- read_htseq_human(file.path(HTSEQ_DIR, files[tissue_n]))
}

# --------------------
# Preprocessing counts for QC plots
# --------------------
if (org == 'mouse') {
  counts <- counts %>%
    mutate(
      # Add log2 transformed counts
      logcounts = log2(count +1),
      # Extract experimental condition from sample name
      condition = sub("(G[0-9]+).*", "\\1", sample)
      ) %>%
    group_by(sample) %>%
    ungroup()
} else if (org == 'human') {
  sample_ids <- counts$sample
  subject_ids <- sub("(GTEX-[A-Z0-9]+).*", "\\1", sample_ids)
  counts$subject_id <- subject_ids
  
  subject_meta_df <- readRDS(paste0("./data/gtex_v10_", tissue, "_subject_meta.rds"))
  
  counts <- left_join(counts, subject_meta_df, by = c("subject_id" = "subjectId"))
  
  counts <- counts %>%
    mutate(
      # Add log2 transformed counts
      logcounts = log2(count +1),
      # Extract experimental condition from sample name
      condition = dplyr::case_when(
        sex == "female" & ageBracket %in% c("20-29", "30-39") ~ "F_20_39",
        sex == "female" & ageBracket %in% c("40-49")          ~ "F_40_49",
        sex == "female" & ageBracket %in% c("50-59")          ~ "F_50_59",
        sex == "female" & ageBracket %in% c("60-69", "70-79") ~ "F_60_79",
        sex == "male"   & ageBracket %in% c("20-29", "30-39") ~ "M_20_39",
        sex == "male"   & ageBracket %in% c("40-49")          ~ "M_40_49",
        sex == "male"   & ageBracket %in% c("50-59")          ~ "M_50_59",
        sex == "male"   & ageBracket %in% c("60-69", "70-79") ~ "M_60_79",
        TRUE ~ NA_character_
      )
    ) %>%
    group_by(sample) %>%
    ungroup()
}

counts$condition <- factor(counts$condition,
                           levels = unique(counts$condition))


# --------------------
# Distribution plots
# --------------------
# Density plot of log2(counts + 1) by condition
density_plot <- ggplot(counts, aes(x=logcounts, color=condition)) +
  geom_density(alpha=0.6, size = 1.2) +
  labs(x="log2(counts + 1)",
       y="Density",
       title="Count distribution by condition") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 22, face = "bold"),      # título grande
    axis.title = element_text(size = 18),                     # títulos de ejes
    axis.text = element_text(size = 14),                      # números y etiquetas eje
    legend.title = element_text(size = 16),                   # título de leyenda
    legend.text = element_text(size = 14)                     # texto de leyenda
  )

dir.create(file.path(OUT_DIR), showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = file.path(OUT_DIR, "01A_density_plot.png"),
  plot = density_plot,
  width = 15, height = 8, units = "in", dpi = 120
)

# Histogram of log2(counts + 1) by condition (overlaid)
histogram <- ggplot(counts, aes(x=logcounts, fill=condition)) +
  geom_histogram(bins = 100, alpha = 0.6, position = "identity") +
  labs(x="log2(counts + 1)",
       y="Number of genes",
       title="Count distribution by condition") + 
  theme_minimal() +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

ggsave(
  filename = file.path(OUT_DIR, "01B_histogram.png"),
  plot = histogram,
  width = 15, height = 8, units = "in", dpi = 120
)

# Faceted histograms per condition
faceted_histogram <- ggplot(counts, aes(x = logcounts, colour = condition, fill = condition)) +
  geom_histogram(bins = 100, alpha = 0.6) +
  theme_bw() +
  facet_grid(. ~ condition) +
  theme(legend.position = "none") +
  labs(title = "log2(counts + 1) per condition",
       x = "log2(counts + 1)",
       y = "Number of genes") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

ggsave(
  filename = file.path(OUT_DIR, "01C_faceted_histogram.png"),
  plot = faceted_histogram,
  width = 25, height = 8, units = "in", dpi = 120
)


# --------------------
# Calculate zero and low count genes per sample
# --------------------
# Define low-count threshold (e.g. raw counts < 10)
low_threshold <- 10

per_sample_counts <- counts %>%
  group_by(subject_id) %>%
  summarize(
    total_genes   = n(),
    zero_counts   = sum(count == 0),
    low_counts    = sum(count < low_threshold),
    zero_pct      = round((zero_counts / total_genes) * 100, 2),
    low_pct       = round((low_counts / total_genes) * 100, 2)
  )

# --------------------
# Normality tests
# --------------------
# Get list of conditions
conds <- unique(counts$condition)

# QQ plots per condition (on logcounts)
# Interpretation: strong deviations from the straight line indicate non-normality
qqplots <- lapply(conds, function(cond) {
  g <- ggqqplot(counts$logcounts[counts$condition == cond],
                title = paste("QQ Plot -", cond))
  return(g)
})

qqplots_comb <- ggarrange(plotlist = qqplots, 
          ncol = 4, 
          nrow = ceiling(length(qqplots)/4))

ggsave(
  filename = file.path(OUT_DIR, "02A_qqplots.png"),
  plot = qqplots_comb,
  width = 15, height = 15, units = "in", dpi = 120
)

# Shapiro-Wilk on a random subset of genes per condition
# Note: with thousands of genes (n > 5000), any small deviation from normality becomes significant (p << 0.05)
normality_tests <- lapply(conds, function(cond) {
  vals <- counts$logcounts[counts$condition == cond]
  # Take a random subset to avoid warnings and extreme power
  if (length(vals) > 5000) vals <- sample(vals, 5000)
  res <- shapiro.test(vals)
  data.frame(
    condition = cond,
    W = res$statistic,
    p.value = res$p.value
  )
  }) %>% 
  bind_rows()

print(normality_tests)

# Evaluate variance homogeneity
# Levene's test checks if group variances are similar
levene_res <- leveneTest(logcounts ~ condition, 
                         data = counts, 
                         center = "median")
print(levene_res)


# --------------------
# Non-parametric tests
# --------------------
# Kruskal-Wallis test on logcounts by condition
# Interpretation: tests whether at least one condition differs in distribution (median ranks)
kruskal_res <- kruskal.test(logcounts ~ condition, data = counts)
print(kruskal_res)

# Dunn's test (post hoc test)
# Interpretation: pairwise multiple comparisons of mean rank sums between conditions,
# with p-values adjusted for multiple testing (here: Benjamini-Hochberg / FDR).
dunn_res <- dunn.test(x = counts$logcounts,
                      g = counts$condition,
                      method = "bh",
                      kw = FALSE,  # set TRUE if you also want to recompute Kruskal
                      list = TRUE) # return results in list form
print(dunn_res)


# --------------------
# QC visualization
# --------------------
# Calculate boxplot statistics for each condition
box_stats <- counts %>%
  group_by(condition) %>%
  summarize(
    Q1 = quantile(logcounts, 0.25),        # First quartile
    Q3 = quantile(logcounts, 0.75),        # Third quartile
    lower_whisker = min(logcounts[logcounts >= Q1 - 1.5 * (Q3 - Q1)]), # Lower whisker
    upper_whisker = max(logcounts[logcounts <= Q3 + 1.5 * (Q3 - Q1)]), # Upper whisker
    Mean = mean(logcounts),               # Mean for each condition
    .groups = "drop"
  ) %>%
  mutate(cond_x = as.numeric(factor(condition, levels = levels(counts$condition))))

whisker_width <- 0.18  # Horizontal length of whisker caps

boxplot <- ggplot(counts, aes(x = condition, y = logcounts, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +   # sin puntos negros extra, solo caja + bigotes
  stat_summary(fun = mean,
               geom = "point",
               shape = 20,
               size = 3,
               color = "red") +
  geom_segment(
    data = box_stats,
    aes(x = cond_x - whisker_width,
        xend = cond_x + whisker_width,
        y = lower_whisker,
        yend = lower_whisker),
    inherit.aes = FALSE,
    color = "black", size = 0.4
  ) +
  geom_segment(
    data = box_stats,
    aes(x = cond_x - whisker_width,
        xend = cond_x + whisker_width,
        y = upper_whisker,
        yend = upper_whisker),
    inherit.aes = FALSE,
    color = "black", size = 0.4
  ) +
  labs(title = "log2(counts + 1) by condition",
       x = "Condition",
       y = "log2(counts + 1)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )

ggsave(
  filename = file.path(OUT_DIR, "04A_boxplot.png"),
  plot = boxplot,
  width = 15, height = 8, units = "in", dpi = 120
)
# --------------------
# Save statistical test results
# --------------------
# 1. Normality tests (Shapiro-Wilk) - TSV format
write.table(normality_tests,
            file = file.path(OUT_DIR, "02B_shapiro_wilk_test.tsv"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

# 2. Levene's test for homogeneity of variance - Plain text
sink(file.path(OUT_DIR, "02C_levene_test.txt"))
cat("Levene's Test for Homogeneity of Variance\n")
cat("==========================================\n\n")
print(levene_res)
sink()

# 3. Kruskal-Wallis test - Plain text
sink(file.path(OUT_DIR, "03A_kruskal_wallis_test.txt"))
cat("Kruskal-Wallis Test\n")
cat("===================\n\n")
print(kruskal_res)
sink()

# 4. Dunn's test pairwise comparisons - TSV format (rich table)
dunn_table <- data.frame(
  Comparison = as.vector(dunn_res$comparisons),
  Z          = dunn_res$Z,
  P.unadj    = dunn_res$P,
  P.adj      = dunn_res$P.adjusted
)

# Añadir columnas Group1, Group2 y Direction
dunn_table_sep <- dunn_table %>%
  separate(Comparison, into = c("Group1", "Group2"), sep = " - ")

cond_means <- counts %>%
  group_by(condition) %>%
  summarize(mean_log = mean(logcounts), .groups = "drop")

dunn_table_final <- dunn_table_sep %>%
  left_join(cond_means, by = c("Group1" = "condition")) %>%
  dplyr::rename(mean1 = mean_log) %>%
  left_join(cond_means, by = c("Group2" = "condition")) %>%
  dplyr::rename(mean2 = mean_log) %>%
  mutate(
    Significant = ifelse(P.adj < 0.05, "*", " "),
    Direction = case_when(
      mean1 > mean2 ~ paste0(Group1, " > ", Group2),
      mean1 < mean2 ~ paste0(Group2, " > ", Group1),
      TRUE          ~ "~ equal"
    )
  ) %>%
  dplyr::select(Group1, Group2, Z, P.unadj, P.adj, Significant, Direction)

write.table(dunn_table_final,
            file = file.path(OUT_DIR, "03B_dunn_posthoc_results.tsv"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

# 5. Summary statistics per condition (from boxplot) - TSV format
box_stats_save <- box_stats %>%
  dplyr::select(-cond_x)  # Remove internal x-position variable

write.table(box_stats_save,
            file = file.path(OUT_DIR, "04B_summary_statistics_by_condition.tsv"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

# 6. Overall summary report - Plain text
sink(file.path(OUT_DIR, "05_statistical_tests_summary.txt"))
cat("EXPLORATORY STATISTICAL ANALYSIS SUMMARY\n")
cat("========================================\n\n")
cat("Dataset:", mitocarta, "\n")
cat("Number of conditions:", length(conds), "\n")
cat("Conditions:", paste(conds, collapse = ", "), "\n")
cat("Total observations:", nrow(counts), "\n\n")

cat("1. NORMALITY ASSESSMENT (Shapiro-Wilk)\n")
cat("---------------------------------------\n")
print(normality_tests)
cat("\n")
# Evaluate normality for each condition
n_reject <- sum(normality_tests$p.value < 0.05)
n_total <- nrow(normality_tests)

if (n_reject == n_total) {
  cat("Interpretation: ALL conditions reject normality (p < 0.05)\n")
  cat("All", n_total, "conditions show non-normal distributions.\n")
} else if (n_reject == 0) {
  cat("Interpretation: NO conditions reject normality (p >= 0.05)\n")
  cat("All", n_total, "conditions appear normally distributed.\n")
} else {
  cat("Interpretation:", n_reject, "out of", n_total, "conditions reject normality (p < 0.05)\n")
  cat("Conditions with normal distribution:\n")
  normal_conds <- normality_tests %>% filter(p.value >= 0.05) %>% pull(condition)
  cat(paste(normal_conds, collapse = ", "), "\n")
  cat("\nConditions with non-normal distribution:\n")
  non_normal_conds <- normality_tests %>% filter(p.value < 0.05) %>% pull(condition)
  cat(paste(non_normal_conds, collapse = ", "), "\n")
}
cat("Result saved in '02B_shapiro_wilk_test.txt'.\n")
cat("\n\n")

cat("2. HOMOGENEITY OF VARIANCE (Levene's Test)\n")
cat("-------------------------------------------\n")
print(levene_res)
cat("\n")
if (levene_res$`Pr(>F)`[1] < 0.05) {
  cat("Interpretation: Variances are NOT homogeneous across conditions (p < 0.05)\n")
} else {
  cat("Interpretation: Variances are homogeneous across conditions (p >= 0.05)\n")
}
cat("Result saved in '02C_levene_test.txt'.\n")
cat("\n\n")

cat("3. GLOBAL TEST (Kruskal-Wallis)\n")
cat("--------------------------------\n")
print(kruskal_res)
cat("\n")
if (kruskal_res$p.value < 0.05) {
  cat("Interpretation: At least one condition differs significantly (p < 0.05)\n")
} else {
  cat("Interpretation: No significant differences detected (p >= 0.05)\n")
}
cat("Result saved in '03A_kruskal_wallis_test.txt'.\n")
cat("\n\n")

cat("4. POST HOC TEST (Dunn's Test)\n")
cat("-------------------------------\n")
cat("Number of pairwise comparisons:", nrow(dunn_table), "\n")
cat("Significant comparisons (p.adj < 0.05):", sum(dunn_table$P.adj < 0.05), "\n")
cat("Non-significant comparisons (p.adj >= 0.05):", sum(dunn_table$P.adj >= 0.05), "\n")
cat("See '03B_dunn_posthoc_results.tsv' for all pairwise comparisons,\n")
cat("including Z statistic, raw and adjusted p-values, and direction\n")
cat("of the effect based on mean log2(counts + 1).\n")

sink()

cat("\n==> Statistical test results saved to:", OUT_DIR, "\n")

