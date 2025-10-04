# Statistical Analysis
# This script performs statistical analyses for microgastropods as environmental indicators

# Load required packages
library(tidyverse)
library(here)
library(vegan)     # Community ecology analyses
library(ggplot2)   # Plotting
library(corrplot)  # Correlation plots

# Set up paths
processed_data_path <- here("data", "processed")
outputs_path <- here("outputs")
figures_path <- here("outputs", "figures")
tables_path <- here("outputs", "tables")

# Ensure output directories exist
for (path in c(outputs_path, figures_path, tables_path)) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

# Load processed data
cat("Loading processed data...\n")

# Uncomment when data is available:
# community_matrix <- read_csv(
#   file.path(processed_data_path, "community_matrix.csv"),
#   show_col_types = FALSE
# ) %>%
#   column_to_rownames("sample_id")
# 
# environmental_matrix <- read_csv(
#   file.path(processed_data_path, "environmental_matrix.csv"),
#   show_col_types = FALSE
# ) %>%
#   column_to_rownames("site_id")
# 
# diversity_data <- read_csv(
#   file.path(processed_data_path, "diversity_with_metadata.csv"),
#   show_col_types = FALSE
# )

cat("Data loaded successfully.\n")

# Analysis 1: Ordination Analysis (NMDS or PCA)
cat("Performing ordination analysis...\n")

# NMDS ordination
# set.seed(123)  # For reproducibility
# nmds <- metaMDS(community_matrix, distance = "bray", k = 2, trymax = 100)
# 
# # Extract NMDS scores
# nmds_scores <- as.data.frame(scores(nmds))
# nmds_scores$site_id <- rownames(nmds_scores)

# Analysis 2: Environmental Fitting
cat("Testing environmental correlations...\n")

# Fit environmental variables to ordination
# env_fit <- envfit(nmds, environmental_matrix, permutations = 999)
# 
# # Save results
# env_fit_results <- data.frame(
#   variable = names(env_fit$vectors$pvals),
#   r2 = env_fit$vectors$r,
#   p_value = env_fit$vectors$pvals
# )

# Analysis 3: Diversity ~ Environment relationships
cat("Testing diversity-environment relationships...\n")

# Example: Shannon diversity vs environmental variables
# Linear models for each environmental variable
# lm_results <- list()
# env_vars <- names(environmental_matrix)
# 
# for (var in env_vars) {
#   formula <- as.formula(paste("shannon ~", var))
#   model <- lm(formula, data = diversity_data)
#   lm_results[[var]] <- summary(model)
# }

# Analysis 4: Indicator Species Analysis
cat("Performing indicator species analysis...\n")

# Identify species associated with environmental conditions
# This requires grouping sites by environmental conditions
# Example using multipatt from indicspecies package:
# library(indicspecies)
# 
# # Create groups based on environmental variable (e.g., pollution level)
# site_groups <- ifelse(environmental_matrix$pollution > median(environmental_matrix$pollution),
#                       "High", "Low")
# 
# # Run indicator species analysis
# indval <- multipatt(community_matrix, site_groups, func = "r.g", control = how(nperm=999))

# Analysis 5: Correlation Analysis
cat("Analyzing correlations between variables...\n")

# Correlation between diversity indices and environmental variables
# cor_matrix <- cor(
#   diversity_data[, c("richness", "shannon", "simpson", "evenness")],
#   environmental_matrix,
#   use = "complete.obs"
# )
# 
# # Statistical significance of correlations
# cor_test_results <- list()
# for (div_index in c("richness", "shannon", "simpson", "evenness")) {
#   for (env_var in env_vars) {
#     test <- cor.test(diversity_data[[div_index]], environmental_matrix[[env_var]])
#     cor_test_results[[paste(div_index, env_var, sep = "_")]] <- test
#   }
# }

# Save analysis results
cat("Saving analysis results...\n")

# write_csv(
#   env_fit_results,
#   file.path(tables_path, "environmental_fit_results.csv")
# )
# 
# write_csv(
#   as.data.frame(cor_matrix),
#   file.path(tables_path, "diversity_environment_correlations.csv")
# )

cat("Statistical analysis completed successfully!\n")
cat("Results saved to:", tables_path, "\n")
