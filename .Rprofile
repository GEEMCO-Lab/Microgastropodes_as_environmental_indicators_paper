# .Rprofile for Microgastropods Research Project
# This file is automatically loaded when R starts in this project

# Display welcome message
cat("\n")
cat("================================================\n")
cat("Microgastropods as Environmental Indicators\n")
cat("================================================\n")
cat("\n")
cat("Quick Start:\n")
cat("  source('setup.R')           # First-time setup\n")
cat("  source('scripts/run_all.R') # Run all analyses\n")
cat("\n")
cat("See QUICKSTART.md for detailed instructions\n")
cat("================================================\n")
cat("\n")

# Set options for better display
options(
  width = 80,
  repos = c(CRAN = "https://cran.r-project.org"),
  stringsAsFactors = FALSE,
  scipen = 999  # Avoid scientific notation
)

# Auto-load commonly used packages
.First <- function() {
  if (interactive()) {
    # Try to load here package for path management
    if (require(here, quietly = TRUE)) {
      cat("✓ Project root:", here::here(), "\n\n")
    }
  }
}
