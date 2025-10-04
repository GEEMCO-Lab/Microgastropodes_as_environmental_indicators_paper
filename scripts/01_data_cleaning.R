# Data Cleaning and Validation
# This script loads raw data, performs quality checks, and creates cleaned datasets

# Load required packages
library(tidyverse)
library(here)

# Set up paths
raw_data_path <- here("data", "raw")
processed_data_path <- here("data", "processed")

# Ensure output directory exists
if (!dir.exists(processed_data_path)) {
  dir.create(processed_data_path, recursive = TRUE)
}

# Function to load and validate data
load_and_validate <- function(filename, required_cols = NULL) {
  filepath <- file.path(raw_data_path, filename)
  
  if (!file.exists(filepath)) {
    stop(paste("File not found:", filepath))
  }
  
  data <- read_csv(filepath, show_col_types = FALSE)
  
  # Check for required columns if specified
  if (!is.null(required_cols)) {
    missing_cols <- setdiff(required_cols, names(data))
    if (length(missing_cols) > 0) {
      stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
    }
  }
  
  return(data)
}

# Main data loading and cleaning
cat("Loading raw data...\n")

# Load microgastropod samples data
# Uncomment and modify when data is available:
# samples_data <- load_and_validate(
#   "microgastropods_samples.csv",
#   required_cols = c("sample_id", "site_id", "date", "latitude", "longitude")
# )

# Load environmental variables
# environmental_data <- load_and_validate(
#   "environmental_variables.csv",
#   required_cols = c("site_id", "variable", "value")
# )

# Load species abundance data
# species_data <- load_and_validate(
#   "species_abundance.csv",
#   required_cols = c("sample_id", "species", "abundance")
# )

cat("Raw data loaded successfully.\n")

# Data cleaning steps
cat("Performing data cleaning...\n")

# 1. Remove duplicates
# samples_data_clean <- samples_data %>%
#   distinct()

# 2. Handle missing values
# Check for NA values and decide on strategy (remove, impute, etc.)
# samples_data_clean <- samples_data_clean %>%
#   filter(!is.na(sample_id))

# 3. Standardize formats
# Convert dates to Date format
# samples_data_clean <- samples_data_clean %>%
#   mutate(date = as.Date(date))

# 4. Validate ranges
# Check that numeric values are within expected ranges
# environmental_data_clean <- environmental_data %>%
#   filter(value >= 0)

# 5. Remove outliers if necessary
# Use appropriate methods based on your data

cat("Data cleaning completed.\n")

# Save cleaned data
cat("Saving cleaned data...\n")

# write_csv(
#   samples_data_clean,
#   file.path(processed_data_path, "cleaned_samples_data.csv")
# )
# 
# write_csv(
#   environmental_data_clean,
#   file.path(processed_data_path, "cleaned_environmental_data.csv")
# )
# 
# write_csv(
#   species_data_clean,
#   file.path(processed_data_path, "cleaned_species_data.csv")
# )

cat("Data cleaning script completed successfully!\n")
cat("Cleaned data saved to:", processed_data_path, "\n")
