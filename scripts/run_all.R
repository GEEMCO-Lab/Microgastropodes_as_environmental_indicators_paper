# Master Script - Run All Analyses
# This script executes all analysis scripts in the correct order

# Load required package
library(here)

# Print session info for reproducibility
cat("==============================================\n")
cat("Running Master Analysis Script\n")
cat("==============================================\n\n")

cat("Session Info:\n")
print(sessionInfo())
cat("\n")

# Record start time
start_time <- Sys.time()

# Function to run script and report status
run_script <- function(script_name) {
  cat("\n==============================================\n")
  cat("Running:", script_name, "\n")
  cat("==============================================\n")
  
  script_path <- here("scripts", script_name)
  
  if (!file.exists(script_path)) {
    stop(paste("Script not found:", script_path))
  }
  
  tryCatch({
    source(script_path)
    cat("\n✓ Successfully completed:", script_name, "\n")
  }, error = function(e) {
    cat("\n✗ Error in", script_name, ":", conditionMessage(e), "\n")
    stop(paste("Stopping execution due to error in", script_name))
  })
}

# Run scripts in order
scripts <- c(
  "01_data_cleaning.R",
  "02_data_processing.R",
  "03_statistical_analysis.R",
  "04_create_figures.R"
)

for (script in scripts) {
  run_script(script)
}

# Calculate and report total runtime
end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")

cat("\n==============================================\n")
cat("All analyses completed successfully!\n")
cat("Total runtime:", round(runtime, 2), "minutes\n")
cat("==============================================\n")
