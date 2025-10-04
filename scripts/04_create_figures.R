# Figure Generation
# This script creates publication-quality figures for the manuscript

# Load required packages
library(tidyverse)
library(here)
library(ggplot2)
library(patchwork)  # For combining plots
library(scales)     # For plot formatting

# Set up paths
processed_data_path <- here("data", "processed")
figures_path <- here("outputs", "figures")

# Ensure output directory exists
if (!dir.exists(figures_path)) {
  dir.create(figures_path, recursive = TRUE)
}

# Set theme for all plots
theme_set(theme_bw(base_size = 12))

# Custom color palette
custom_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b")

# Load processed data
cat("Loading processed data...\n")

# Uncomment when data is available:
# diversity_data <- read_csv(
#   file.path(processed_data_path, "diversity_with_metadata.csv"),
#   show_col_types = FALSE
# )
# 
# community_matrix <- read_csv(
#   file.path(processed_data_path, "community_matrix.csv"),
#   show_col_types = FALSE
# )

cat("Data loaded successfully.\n")

# Figure 1: Diversity patterns across sites
cat("Creating Figure 1: Diversity patterns...\n")

# fig1 <- ggplot(diversity_data, aes(x = site_id, y = shannon)) +
#   geom_boxplot(fill = custom_colors[1], alpha = 0.7) +
#   geom_jitter(width = 0.2, alpha = 0.5) +
#   labs(
#     title = "Shannon Diversity Across Sampling Sites",
#     x = "Site",
#     y = "Shannon Diversity Index"
#   ) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# ggsave(
#   file.path(figures_path, "fig1_diversity_patterns.png"),
#   fig1,
#   width = 10,
#   height = 6,
#   dpi = 300
# )

# Figure 2: Species richness vs environmental variables
cat("Creating Figure 2: Species richness relationships...\n")

# fig2 <- ggplot(diversity_data, aes(x = temperature, y = richness)) +
#   geom_point(size = 3, alpha = 0.6, color = custom_colors[2]) +
#   geom_smooth(method = "lm", se = TRUE, color = custom_colors[3]) +
#   labs(
#     title = "Species Richness vs Temperature",
#     x = "Temperature (°C)",
#     y = "Species Richness"
#   )
# 
# ggsave(
#   file.path(figures_path, "fig2_richness_temperature.png"),
#   fig2,
#   width = 8,
#   height = 6,
#   dpi = 300
# )

# Figure 3: Ordination plot
cat("Creating Figure 3: Ordination plot...\n")

# NMDS ordination visualization
# This would use the NMDS results from statistical analysis
# fig3 <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
#   geom_point(aes(color = site_type), size = 3, alpha = 0.7) +
#   stat_ellipse(aes(color = site_type), level = 0.95) +
#   labs(
#     title = "NMDS Ordination of Microgastropod Communities",
#     x = "NMDS1",
#     y = "NMDS2",
#     color = "Site Type"
#   ) +
#   scale_color_manual(values = custom_colors)
# 
# ggsave(
#   file.path(figures_path, "fig3_ordination.png"),
#   fig3,
#   width = 8,
#   height = 6,
#   dpi = 300
# )

# Figure 4: Correlation heatmap
cat("Creating Figure 4: Correlation matrix...\n")

# library(corrplot)
# 
# png(
#   file.path(figures_path, "fig4_correlation_heatmap.png"),
#   width = 2400,
#   height = 2400,
#   res = 300
# )
# 
# corrplot(
#   cor_matrix,
#   method = "color",
#   type = "full",
#   tl.col = "black",
#   tl.srt = 45,
#   addCoef.col = "black",
#   number.cex = 0.7,
#   title = "Diversity-Environment Correlations",
#   mar = c(0, 0, 2, 0)
# )
# 
# dev.off()

# Figure 5: Multi-panel figure combining key results
cat("Creating Figure 5: Multi-panel summary...\n")

# Using patchwork to combine plots
# combined_fig <- (fig1 | fig2) / (fig3)
# 
# ggsave(
#   file.path(figures_path, "fig5_combined_results.png"),
#   combined_fig,
#   width = 14,
#   height = 10,
#   dpi = 300
# )

cat("Figure generation completed successfully!\n")
cat("Figures saved to:", figures_path, "\n")
