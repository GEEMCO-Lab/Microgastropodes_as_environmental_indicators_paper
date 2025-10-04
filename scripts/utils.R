# Utility Functions
# Common functions used across multiple analysis scripts

# Function to check and install packages
check_and_install <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      message(paste("Installing package:", pkg))
      install.packages(pkg, repos = "https://cran.r-project.org")
      library(pkg, character.only = TRUE)
    }
  }
}

# Function to create output directories
create_output_dirs <- function(base_path = here::here()) {
  dirs <- c(
    file.path(base_path, "data", "raw"),
    file.path(base_path, "data", "processed"),
    file.path(base_path, "outputs", "figures"),
    file.path(base_path, "outputs", "tables")
  )
  
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      message(paste("Created directory:", dir))
    }
  }
}

# Function to save plot with consistent settings
save_plot <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  ggsave(
    filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  message(paste("Saved plot:", filename))
}

# Function to export table to multiple formats
export_table <- function(data, base_filename, output_dir) {
  # CSV
  readr::write_csv(data, file.path(output_dir, paste0(base_filename, ".csv")))
  
  # Excel (if openxlsx is available)
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    openxlsx::write.xlsx(data, file.path(output_dir, paste0(base_filename, ".xlsx")))
  }
  
  message(paste("Saved table:", base_filename))
}

# Function to create summary statistics table
create_summary_table <- function(data, group_var = NULL) {
  if (is.null(group_var)) {
    summary_stats <- data %>%
      summarise(across(where(is.numeric), list(
        mean = ~mean(.x, na.rm = TRUE),
        sd = ~sd(.x, na.rm = TRUE),
        median = ~median(.x, na.rm = TRUE),
        min = ~min(.x, na.rm = TRUE),
        max = ~max(.x, na.rm = TRUE),
        n = ~sum(!is.na(.x))
      )))
  } else {
    summary_stats <- data %>%
      group_by(across(all_of(group_var))) %>%
      summarise(across(where(is.numeric), list(
        mean = ~mean(.x, na.rm = TRUE),
        sd = ~sd(.x, na.rm = TRUE),
        median = ~median(.x, na.rm = TRUE),
        min = ~min(.x, na.rm = TRUE),
        max = ~max(.x, na.rm = TRUE),
        n = ~sum(!is.na(.x))
      )), .groups = "drop")
  }
  
  return(summary_stats)
}

# Function to check data quality
check_data_quality <- function(data, data_name = "Dataset") {
  cat("\nData Quality Check for:", data_name, "\n")
  cat("Dimensions:", nrow(data), "rows x", ncol(data), "columns\n")
  cat("\nMissing values:\n")
  print(colSums(is.na(data)))
  cat("\nDuplicate rows:", sum(duplicated(data)), "\n")
  cat("\nColumn types:\n")
  print(sapply(data, class))
  cat("\n")
}

# Function to standardize column names
standardize_names <- function(data) {
  names(data) <- names(data) %>%
    tolower() %>%
    gsub(" ", "_", .) %>%
    gsub("[^a-z0-9_]", "", .)
  return(data)
}
