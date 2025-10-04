# Project Setup Script
# Run this script once when you first clone the repository

cat("==============================================\n")
cat("Setting up Microgastropods Research Project\n")
cat("==============================================\n\n")

# Function to check if package is installed
check_package <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    return(FALSE)
  }
  return(TRUE)
}

# List of required packages
required_packages <- c(
  "tidyverse",
  "here",
  "vegan",
  "ggplot2",
  "patchwork",
  "corrplot",
  "scales",
  "knitr",
  "rmarkdown"
)

# Check which packages are missing
cat("Checking for required packages...\n")
missing_packages <- c()
for (pkg in required_packages) {
  if (!check_package(pkg)) {
    missing_packages <- c(missing_packages, pkg)
  }
}

# Install missing packages
if (length(missing_packages) > 0) {
  cat("\nThe following packages need to be installed:\n")
  cat(paste("-", missing_packages, collapse = "\n"), "\n\n")
  
  response <- readline(prompt = "Would you like to install them now? (y/n): ")
  
  if (tolower(response) == "y") {
    cat("\nInstalling packages...\n")
    install.packages(missing_packages, repos = "https://cran.r-project.org")
    cat("\nPackages installed successfully!\n")
  } else {
    cat("\nSkipping package installation.\n")
    cat("You can install them later using:\n")
    cat("install.packages(c(", paste0('"', missing_packages, '"', collapse = ", "), "))\n")
  }
} else {
  cat("✓ All required packages are already installed!\n")
}

# Check/create directory structure
cat("\nChecking directory structure...\n")
library(here)

dirs_to_check <- c(
  here("data", "raw"),
  here("data", "processed"),
  here("scripts"),
  here("outputs", "figures"),
  here("outputs", "tables")
)

for (dir in dirs_to_check) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat("✓ Created:", dir, "\n")
  } else {
    cat("✓ Exists:", dir, "\n")
  }
}

# Display next steps
cat("\n==============================================\n")
cat("Setup Complete!\n")
cat("==============================================\n\n")

cat("Next steps:\n")
cat("1. Place your raw data files in data/raw/\n")
cat("   - See data/raw/README.md for expected file formats\n")
cat("   - Use the *_TEMPLATE.csv files as guides\n\n")
cat("2. Run the analysis:\n")
cat("   source('scripts/run_all.R')\n\n")
cat("3. Generate the manuscript:\n")
cat("   rmarkdown::render('manuscript.Rmd')\n")
cat("   # or\n")
cat("   quarto::quarto_render('manuscript.qmd')\n\n")

cat("For more information, see README.md\n")
cat("==============================================\n")
