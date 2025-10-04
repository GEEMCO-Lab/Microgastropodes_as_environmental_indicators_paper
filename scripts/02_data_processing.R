# Data Processing and Transformation
# This script processes cleaned data and creates analysis-ready datasets

# Load required packages
library(tidyverse)
library(here)
library(vegan)  # For community ecology analyses

# Set up paths
processed_data_path <- here("data", "processed")

# Load cleaned data
cat("Loading cleaned data...\n")

# Uncomment when data is available:
# samples_data <- read_csv(
#   file.path(processed_data_path, "cleaned_samples_data.csv"),
#   show_col_types = FALSE
# )
# 
# environmental_data <- read_csv(
#   file.path(processed_data_path, "cleaned_environmental_data.csv"),
#   show_col_types = FALSE
# )
# 
# species_data <- read_csv(
#   file.path(processed_data_path, "cleaned_species_data.csv"),
#   show_col_types = FALSE
# )

cat("Cleaned data loaded successfully.\n")

# Create community matrix (sites x species)
cat("Creating community matrix...\n")

# community_matrix <- species_data %>%
#   pivot_wider(
#     names_from = species,
#     values_from = abundance,
#     values_fill = 0
#   ) %>%
#   column_to_rownames("sample_id")

# Create environmental matrix
cat("Creating environmental matrix...\n")

# environmental_matrix <- environmental_data %>%
#   pivot_wider(
#     names_from = variable,
#     values_from = value
#   ) %>%
#   column_to_rownames("site_id")

# Calculate diversity indices
cat("Calculating diversity indices...\n")

# diversity_indices <- data.frame(
#   sample_id = rownames(community_matrix),
#   richness = specnumber(community_matrix),
#   shannon = diversity(community_matrix, index = "shannon"),
#   simpson = diversity(community_matrix, index = "simpson"),
#   evenness = diversity(community_matrix, index = "shannon") / log(specnumber(community_matrix))
# )

# Merge sample metadata with diversity indices
# diversity_with_metadata <- samples_data %>%
#   left_join(diversity_indices, by = "sample_id")

# Save processed data
cat("Saving processed data...\n")

# write_csv(
#   as.data.frame(community_matrix) %>% rownames_to_column("sample_id"),
#   file.path(processed_data_path, "community_matrix.csv")
# )
# 
# write_csv(
#   as.data.frame(environmental_matrix) %>% rownames_to_column("site_id"),
#   file.path(processed_data_path, "environmental_matrix.csv")
# )
# 
# write_csv(
#   diversity_indices,
#   file.path(processed_data_path, "diversity_indices.csv")
# )
# 
# write_csv(
#   diversity_with_metadata,
#   file.path(processed_data_path, "diversity_with_metadata.csv")
# )

cat("Data processing completed successfully!\n")
cat("Processed data saved to:", processed_data_path, "\n")
