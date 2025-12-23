library(dplyr)
library(ggplot2)

mitocarta <- "after_mitocarta"
tissue <- "liver"

OUT_DIR <- paste0("./figuras/human/", tissue, "/", mitocarta)

subject_meta_df <- readRDS(paste0("./data/gtex_v10_", tissue, "_subject_meta.rds"))
 
# ggplot(subject_meta_df, aes(x = ageBracket)) +
#   geom_bar() +
#   xlab("Age Bracket (years)") +
#   ylab("Number of donors") +
#   ggtitle("Age distribution of GTEx Liver donors")

df_age <- subject_meta_df %>% 
  count(ageBracket)

if (tissue == "liver") {
  tissue_name <- "Liver"
} else if (tissue == "heart_aa") {
  tissue_name <- "Heart Atrial Appendage"
} else if (tissue == "heart_lv") {
  tissue_name <- "Heart Left Ventricule"
}

distribution_plot <- ggplot(df_age, aes(x = ageBracket, y = n, group = 1)) +
  geom_col(fill = "#1874CD") +
  geom_line(color = "#FF4500", size = 1.2) +
  geom_point(color = "#FF4500", size = 2) +
  xlab("Age Bracket (years)") +
  ylab("Number of donors") +
  ggtitle(paste0("Age distribution of GTEx ", tissue_name, " donors")) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 22),
    axis.title = element_text(size = 25),
    text = element_text(size = 18),
    plot.title = element_text(size = 25, face = "bold")
  )

dir.create(file.path(OUT_DIR, "population_distribution"), showWarnings = FALSE, recursive = TRUE)

ggsave(filename = file.path(OUT_DIR, "population_distribution", "age_distribution.png"),
       plot = distribution_plot,
       width = 12, height = 8, units = "in", dpi = 120
       )
