# R Script: microgastropodes_environmental_main_analysis.R
# Author: Eduardo HC Galvao - eduardohcgalvao@gmail.com
# Version: 6.0 09/10/2025
# Date: 30/09/2025 (dd/mm/yyyy) #nolint
# Description: Analysis of microgastropod community in relation to environmental
# variables across wet and dry seasons. This script performs VIF-based variable
# selection, distance-based redundancy analysis (dbRDA), NMDS ordination with
# environmental fitting, PERMANOVA and ANOSIM tests, diversity comparison:
# shannon, simpson, richness, pielou evenness and margalef richness and
# SIMPER analysis to identify species driving seasonal differences.
# Additionally, it includes GAM models relating diversity indices to
# environmental variables.

# Input: environmental_data.csv, environmental_data_average_macronutrients.csv,
# microgastropodes_abundance_matrix.csv

# If one the following packages are not installed run:
# install.packages("package") #nolint
library(ggplot2) # For plotting
library(dplyr) # For data manipulation
library(tidyr) # For data reshaping
library(reshape2) # For data manipulation
library(car) # For VIF analysis and additional statistical tests
library(vegan)     # For dbRDA, PERMANOVA, diversity indices
library(cowplot)   # For combining plots
library(mgcv) # For GAM models

set.seed(123) # For reproducibility

# Analysis workflow:
# Data loading and preparation
# Environmental variable barplots
# VIF analysis for multicollinearity
# dbRDA (distance-based redundancy analysis)
# NMDS ordination with environmental fitting
# PERMANOVA tests
# ANOSIM analysis
# Diversity indices analysis (Shannon, Simpson, Richness, Evenness, Margalef)
# SIMPER analysis
# Models for diversity indices vs environmental variables

###########################################################################
######################### Load and prepare data ###########################
###########################################################################

# Load abiotic and biotic data
env_macro_data <- read.csv("C:/Users/Kelmo/Desktop/analises_eduardo/dados/environmental_data_average_macronutrients.csv")
env_var_data <- read.csv("C:/Users/Kelmo/Desktop/analises_eduardo/dados/environmental_variables_data.csv")
community_data <- read.csv("C:/Users/Kelmo/Desktop/analises_eduardo/dados/microgastropodes_abundance_matrix.csv")

# Transform biotic data from wide to long format
# Column names: species, then sampling points with season codes (A=Wet, O=Dry)
species_names <- community_data[,1]
abundance_data <- community_data[,-1]

# Create sampling point information
sampling_points <- colnames(abundance_data)
site_codes <- c("C5", "C6", "C7", "C8", "C9", "C10", "N1", "N2", "N3", "N4", "S1", "S2", "S3", "S4", "S5")

# Create species matrix for community analysis
# Rows = sampling points, columns = species
community_matrix <- t(abundance_data)
colnames(community_matrix) <- species_names
rownames(community_matrix) <- sampling_points

# Extract season and site information from column names
season_codes <- substr(sampling_points, nchar(sampling_points), nchar(sampling_points))
site_names <- substr(sampling_points, 1, nchar(sampling_points)-1)

# Convert season codes to readable names (A = Wet, O = Dry)
seasons <- ifelse(season_codes == "A", "Wet", "Dry")

# Create site-season dataframe
site_season_data <- data.frame(
  sampling_points = sampling_points,
  site = site_names,
  season = seasons,
  stringsAsFactors = FALSE
)

# Add environment mapping for macronutrients merge
# Map sites to environments based on the data structure
site_to_environment <- data.frame(
  site = c("N1", "N2", "N3", "N4", "C5", "C6", "C7", "C8", "C9", "C10", "S1", "S2", "S3", "S4", "S5"),
  environment = c(rep("Serinhaem", 4), rep("Sorojo", 6), rep("Marau", 5)),
  stringsAsFactors = FALSE
)

# Add environment to site_season_data
site_season_data <- merge(site_season_data, site_to_environment, by = "site", all.x = TRUE)

# Merge both biotic and abiotic tables
# Verify data types
env_var_data$sampling_points <- as.character(env_var_data$sampling_points)
env_var_data$season <- as.factor(env_var_data$season)
env_var_data$temperature <- as.numeric(env_var_data$temperature)
env_var_data$salinity <- as.numeric(env_var_data$salinity)
env_var_data$pH <- as.numeric(env_var_data$pH)
env_var_data$dissolved_o2 <- as.numeric(env_var_data$dissolved_o2)
env_var_data$suspended_solids <- as.numeric(env_var_data$suspended_solids)
env_var_data$sediment_particle_size <- as.numeric(env_var_data$sediment_particle_size)

# Convert season in env_var_data to match format (Wet/Dry with capital letters)
env_var_data$season <- as.character(env_var_data$season)
env_var_data$season <- ifelse(tolower(env_var_data$season) %in% c("wet"), "Wet",
                             ifelse(tolower(env_var_data$season) %in% c("dry"), "Dry", 
                                   env_var_data$season))

# Merge environmental variables data with site-season information
combined_data <- merge(site_season_data, env_var_data, 
                      by.x = c("site", "season"), by.y = c("sampling_points", "season"), 
                      suffixes = c("", "_env_var"))

# Merge macronutrients data with combined data by environment and season
combined_data <- merge(combined_data, env_macro_data, 
                      by = c("environment", "season"), 
                      suffixes = c("", "_macro"))

# Ensure the order matches the community matrix
combined_data <- combined_data[match(site_season_data$sampling_points, combined_data$sampling_points),]

# Create environmental matrix for analyses (including macronutrients)
env_matrix <- combined_data[, c("temperature", "salinity", "pH", "dissolved_o2", "suspended_solids", 
                               "sediment_particle_size", "Mg", "K", "Ca", "P", "Al", "V", "B", "Fe", "Cu")]
rownames(env_matrix) <- combined_data$sampling_points

################################################################################
####################### Environmental Variables Barplots #######################
################################################################################

# Prepare data for environmental barplots (including macronutrients)
env_plot_data <- combined_data %>%
  select(sampling_points, site, season, temperature, salinity, pH, dissolved_o2, 
         suspended_solids, sediment_particle_size, Mg, K, Ca, P, Al, V, B, Fe, Cu) %>%
  pivot_longer(cols = c(temperature, salinity, pH, dissolved_o2, suspended_solids, sediment_particle_size,
                       Mg, K, Ca, P, Al, V, B, Fe, Cu),
               names_to = "variable", values_to = "value") %>%
  mutate(
    variable_label = case_when(
      variable == "temperature" ~ "Temperature (°C)",
      variable == "salinity" ~ "Salinity",
      variable == "pH" ~ "pH",
      variable == "dissolved_o2" ~ "Dissolved O2 (mg/L)",
      variable == "suspended_solids" ~ "Suspended solids (mg/L)",
      variable == "sediment_particle_size" ~ "Sediment particle size (φ)",
      variable == "Mg" ~ "Magnesium (mg/L)",
      variable == "K" ~ "Potassium (mg/L)",
      variable == "Ca" ~ "Calcium (mg/L)",
      variable == "P" ~ "Phosphorus (mg/L)",
      variable == "Al" ~ "Aluminum (mg/L)",
      variable == "V" ~ "Vanadium (mg/L)",
      variable == "B" ~ "Boron (mg/L)",
      variable == "Fe" ~ "Iron (mg/L)",
      variable == "Cu" ~ "Copper (mg/L)",
      TRUE ~ variable
    ),
    variable_label = factor(variable_label, 
                           levels = c("Temperature (°C)", "Salinity", "pH", 
                                    "Dissolved O2 (mg/L)", "Suspended solids (mg/L)", "Sediment particle size (φ)",
                                    "Magnesium (mg/L)", "Potassium (mg/L)", "Calcium (mg/L)", "Phosphorus (mg/L)",
                                    "Aluminum (mg/L)", "Vanadium (mg/L)", "Boron (mg/L)", "Iron (mg/L)", "Copper (mg/L)")),
    # Set the site order as requested: N1, N2, N3, N4, C5, C6, C7, C8, C9, C10, S1, S2, S3, S4, S5
    site = factor(site, levels = c("N1", "N2", "N3", "N4", "C5", "C6", "C7", "C8", "C9", "C10", 
                                  "S1", "S2", "S3", "S4", "S5"))
  )

# Create individual environmental plots for environmental variables

# Temperature plot
temp_data <- env_plot_data %>% filter(variable_label == "Temperature (°C)")
temp_plot <- ggplot(temp_data, aes(x = site, y = value, fill = season)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  scale_fill_manual(values = c("Wet" = "#72309e", "Dry" = "#9ec4e3"), name = "") +
  labs(x = "", y = "Temperature (°C)") +
  theme_classic() +
  coord_cartesian(ylim = c(25, 30)) +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank()
  )

# Salinity plot
sal_data <- env_plot_data %>% filter(variable_label == "Salinity")
sal_plot <- ggplot(sal_data, aes(x = site, y = value, fill = season)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  scale_fill_manual(values = c("Wet" = "#72309e", "Dry" = "#9ec4e3"), name = "") +
  labs(x = "", y = "Salinity") +
  theme_classic() +
  coord_cartesian(ylim = c(5, 40)) +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank()
  )

# pH plot (with adjusted y-axis to show variation clearly)
ph_data <- env_plot_data %>% filter(variable_label == "pH")
ph_plot <- ggplot(ph_data, aes(x = site, y = value, fill = season)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  scale_fill_manual(values = c("Wet" = "#72309e", "Dry" = "#9ec4e3"), name = "") +
  labs(x = "", y = "pH") +
  theme_classic() +
  coord_cartesian(ylim = c(6, 9)) +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank()
  )

# Dissolved O2 plot
do_data <- env_plot_data %>% filter(variable_label == "Dissolved O2 (mg/L)")
do_plot <- ggplot(do_data, aes(x = site, y = value, fill = season)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  scale_fill_manual(values = c("Wet" = "#72309e", "Dry" = "#9ec4e3"), name = "") +
  labs(x = "", y = "Dissolved O2 (mg/L)") +
  theme_classic() +
  coord_cartesian(ylim = c(2, 6)) +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank()
  )

# Suspended solids plot
ss_data <- env_plot_data %>% filter(variable_label == "Suspended solids (mg/L)")
ss_plot <- ggplot(ss_data, aes(x = site, y = value, fill = season)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  scale_fill_manual(values = c("Wet" = "#72309e", "Dry" = "#9ec4e3"), name = "") +
  labs(x = "", y = "Suspended solids (mg/L)") +
  theme_classic() +
  coord_cartesian(ylim = c(15, 40)) +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank()
  )

# Sediment particle size plot
ps_data <- env_plot_data %>% filter(variable_label == "Sediment particle size (φ)")
ps_plot <- ggplot(ps_data, aes(x = site, y = value, fill = season)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  scale_fill_manual(values = c("Wet" = "#72309e", "Dry" = "#9ec4e3"), name = "") +
  labs(x = "", y = "Sediment particle size (φ)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank()
  )

# Combine plots with cowplot and add shared legend
env_plots_grid <- plot_grid(
  temp_plot, do_plot, sal_plot,
  ss_plot, ph_plot, ps_plot,
  ncol = 2, nrow = 3, align = "hv",
  labels = c("A", "B", "C", "D", "E", "F"),
  label_size = 16, label_fontface = "bold"
)

# Create shared legend
legend <- get_legend(
  temp_plot +
    scale_fill_manual(values = c("Wet" = "#72309e", "Dry" = "#9ec4e3"), name = "") +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 14, face = "bold"),
          legend.key.size = unit(1, "cm"))
)

# Combine plots with legend at bottom center
env_plots_grid_leg <- plot_grid(
  env_plots_grid,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.1)
)

# Save the combined plot
#ggsave("environmental_variables_barplot.png", env_plots_grid_leg, width = 14, height = 12, dpi = 1200)

###########################################################################
############## Macronutrients Descriptive Statistics by Environment ######
###########################################################################

# Calculate mean and sd for macronutrients by environment
macro_nutrients <- c("Mg", "K", "Ca", "P", "Al", "V", "B", "Fe", "Cu")
environments <- c("Serinhaem", "Sorojo", "Marau")

# Create empty dataframe to store results
macro_stats_by_env <- data.frame(
  Environment = character(),
  Variable = character(),
  Mean = numeric(),
  SD = numeric(),
  N = numeric(),
  stringsAsFactors = FALSE
)

# Calculate statistics for each environment and macronutrient
for (env in environments) {
  env_data <- combined_data[combined_data$environment == env, ]
  
  for (nutrient in macro_nutrients) {
    if (nutrient %in% colnames(env_data)) {
      nutrient_values <- env_data[[nutrient]]
      nutrient_values <- nutrient_values[!is.na(nutrient_values)]
      
      if (length(nutrient_values) > 0) {
        macro_stats_by_env <- rbind(macro_stats_by_env, data.frame(
          Environment = env,
          Variable = nutrient,
          Mean = round(mean(nutrient_values), 3),
          SD = round(sd(nutrient_values), 3),
          N = length(nutrient_values)
        ))
      }
    }
  }
}

# Display results
#print(macro_stats_by_env)

# Create a more readable table format
macro_wide <- reshape(macro_stats_by_env[, c("Environment", "Variable", "Mean", "SD")], 
                     idvar = "Variable", 
                     timevar = "Environment", 
                     direction = "wide")

#print(macro_wide)

# Save to CSV
#write.csv(macro_stats_by_env, "macronutrients_stats_by_environment.csv", row.names = FALSE)
#write.csv(macro_wide, "macronutrients_summary_table.csv", row.names = FALSE)

# Create summary table with Mean ± SD format
macro_summary <- data.frame(Variable = macro_nutrients)

for (env in environments) {
  env_stats <- macro_stats_by_env[macro_stats_by_env$Environment == env, ]
  mean_sd_values <- paste0(env_stats$Mean, " ± ", env_stats$SD)
  names(mean_sd_values) <- env_stats$Variable
  
  # Match the order of variables
  ordered_values <- mean_sd_values[macro_nutrients]
  macro_summary[[env]] <- ordered_values
}
#print(macro_summary)

#write.csv(macro_summary, "macronutrients_mean_sd_summary.csv", row.names = FALSE)

###########################################################################
############## Macronutrients Summary by Season ##########################
###########################################################################

# Calculate mean and sd for macronutrients by season only (Wet vs Dry)
seasons <- c("Wet", "Dry")

# Create empty dataframe to store seasonal results
macro_stats_by_season <- data.frame(
  Season = character(),
  Variable = character(),
  Mean = numeric(),
  SD = numeric(),
  N = numeric(),
  stringsAsFactors = FALSE
)

# Calculate statistics for each season and macronutrient
for (season in seasons) {
  season_data <- combined_data[combined_data$season == season, ]
  
  for (nutrient in macro_nutrients) {
    if (nutrient %in% colnames(season_data)) {
      nutrient_values <- season_data[[nutrient]]
      nutrient_values <- nutrient_values[!is.na(nutrient_values)]
      
      if (length(nutrient_values) > 0) {
        macro_stats_by_season <- rbind(macro_stats_by_season, data.frame(
          Season = season,
          Variable = nutrient,
          Mean = round(mean(nutrient_values), 3),
          SD = round(sd(nutrient_values), 3),
          N = length(nutrient_values)
        ))
      }
    }
  }
}

# Display seasonal results
#print(macro_stats_by_season)

# Save seasonal results
#write.csv(macro_stats_by_season, "macronutrients_seasonal_stats.csv", row.names = FALSE)

# Create a seasonal summary table with Mean ± SD format
macro_seasonal_summary <- data.frame(Variable = macro_nutrients)

for (season in seasons) {
  season_stats <- macro_stats_by_season[macro_stats_by_season$Season == season, ]
  mean_sd_values <- paste0(season_stats$Mean, " ± ", season_stats$SD)
  names(mean_sd_values) <- season_stats$Variable
  
  # Match the order of variables
  ordered_values <- mean_sd_values[macro_nutrients]
  macro_seasonal_summary[[season]] <- ordered_values
}

#print(macro_seasonal_summary)

# Create a comparison table showing differences between seasons
macro_comparison <- data.frame(
  Variable = macro_nutrients,
  Wet_Mean = numeric(length(macro_nutrients)),
  Dry_Mean = numeric(length(macro_nutrients)),
  Difference = numeric(length(macro_nutrients)),
  Percent_Change = numeric(length(macro_nutrients))
)

for (i in 1:length(macro_nutrients)) {
  nutrient <- macro_nutrients[i]
  wet_mean <- macro_stats_by_season[macro_stats_by_season$Season == "Wet" & 
                                   macro_stats_by_season$Variable == nutrient, "Mean"]
  dry_mean <- macro_stats_by_season[macro_stats_by_season$Season == "Dry" & 
                                   macro_stats_by_season$Variable == nutrient, "Mean"]
  
  if (length(wet_mean) > 0 && length(dry_mean) > 0) {
    macro_comparison$Wet_Mean[i] <- wet_mean
    macro_comparison$Dry_Mean[i] <- dry_mean
    macro_comparison$Difference[i] <- round(wet_mean - dry_mean, 3)
    macro_comparison$Percent_Change[i] <- round(((wet_mean - dry_mean) / dry_mean) * 100, 1)
  }
}

#print(macro_comparison)

# Save seasonal results
#write.csv(macro_stats_by_season, "macronutrients_seasonal_stats.csv", row.names = FALSE)
#write.csv(macro_seasonal_summary, "macronutrients_seasonal_summary.csv", row.names = FALSE)
#write.csv(macro_comparison, "macronutrients_seasonal_comparison.csv", row.names = FALSE)

###########################################################################
############## Macronutrients Heatmap Visualization ######################
###########################################################################

# Calculate mean values for each environment × season combination
macro_heatmap_data <- data.frame(
  Environment_Season = character(),
  Variable = character(),
  Mean_Value = numeric(),
  stringsAsFactors = FALSE
)

# Get unique environment × season combinations
env_season_combinations <- unique(combined_data[, c("environment", "season")])
env_season_combinations <- env_season_combinations[complete.cases(env_season_combinations), ]

# Calculate means for each combination
for (i in 1:nrow(env_season_combinations)) {
  env <- env_season_combinations$environment[i]
  season <- env_season_combinations$season[i]
  env_season_label <- paste0(env, "_", season)
  
  # Filter data for this environment × season combination
  subset_data <- combined_data[combined_data$environment == env & combined_data$season == season, ]
  
  for (nutrient in macro_nutrients) {
    if (nutrient %in% colnames(subset_data)) {
      nutrient_values <- subset_data[[nutrient]]
      nutrient_values <- nutrient_values[!is.na(nutrient_values)]
      
      if (length(nutrient_values) > 0) {
        macro_heatmap_data <- rbind(macro_heatmap_data, data.frame(
          Environment_Season = env_season_label,
          Variable = nutrient,
          Mean_Value = mean(nutrient_values)
        ))
      }
    }
  }
}

# Convert to wide format for heatmap
heatmap_matrix <- reshape(macro_heatmap_data, 
                         idvar = "Variable", 
                         timevar = "Environment_Season", 
                         direction = "wide")

# Clean column names (remove "Mean_Value." prefix)
colnames(heatmap_matrix) <- gsub("Mean_Value.", "", colnames(heatmap_matrix))

# Set row names to variables
rownames(heatmap_matrix) <- heatmap_matrix$Variable
heatmap_matrix$Variable <- NULL

# Convert to matrix for heatmap
heatmap_matrix <- as.matrix(heatmap_matrix)

# Create the heatmap using ggplot2
# First, convert matrix back to long format for ggplot
heatmap_long <- reshape2::melt(heatmap_matrix)
colnames(heatmap_long) <- c("Nutrient", "Environment_Season", "Mean_Value")

# Create environment and season columns for better labeling
heatmap_long$Environment <- gsub("_.*", "", heatmap_long$Environment_Season)
heatmap_long$Season <- gsub(".*_", "", heatmap_long$Environment_Season)

# Create the heatmap plot
macro_heatmap <- ggplot(heatmap_long, aes(x = Environment_Season, y = Nutrient, fill = Mean_Value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", 
                       midpoint = median(heatmap_long$Mean_Value, na.rm = TRUE),
                       name = "Mean\nConcentration") +
  geom_text(aes(label = round(Mean_Value, 2)), color = "black", size = 3.5, fontface = "bold") +
  labs(title = "Macronutrients Concentration Heatmap",
       subtitle = "Mean values by Environment × Season combinations",
       x = "Environment × Season", 
       y = "Macronutrients") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold.italic"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  coord_fixed(ratio = 0.8)

# Display the heatmap
macro_heatmap

# Alternative heatmap with better color scheme for publication
macro_heatmap_alt <- ggplot(heatmap_long, aes(x = Environment_Season, y = Nutrient, fill = Mean_Value)) +
  geom_tile(color = "white", linewidth = 1) +
  scale_fill_viridis_c(name = "Mean\nConcentration", option = "cividis") +
  geom_text(aes(label = round(Mean_Value, 2)), color = "white", size = 3.5, fontface = "bold") +
  labs(title = "Macronutrients Concentration Heatmap",
       subtitle = "Mean values by Environment × Season combinations",
       x = "Environment × Season", 
       y = "Macronutrients") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray40"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, angle = 90, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 16, face = "bold.italic"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "white", fill = NA, linewidth = 1)
  ) +
  coord_fixed(ratio = 0.8)

# Display the alternative heatmap
macro_heatmap_alt

# Save heatmap data and plots
#write.csv(macro_heatmap_data, "macronutrients_heatmap_data.csv", row.names = FALSE)
#write.csv(heatmap_matrix, "macronutrients_heatmap_matrix.csv", row.names = TRUE)
#ggsave("macronutrients_heatmap.png", macro_heatmap, width = 12, height = 8, dpi = 1200)
#ggsave("macronutrients_heatmap_viridis.png", macro_heatmap_alt, width = 12, height = 8, dpi = 1200)

################################################################
######################### VIF Analysis #########################
################################################################

# Check multicollinearity of environmental variables using VIF
# VIF > 7 indicates multicollinearity
env_complete <- na.omit(env_matrix)
community_complete <- community_matrix[rownames(env_complete), ]

# Calculate VIF for each environmental variable
vif_results <- data.frame(Variable = character(), VIF = numeric(), stringsAsFactors = FALSE)

for (i in 1:ncol(env_complete)) {
  # Create a model with all other variables as predictors
  other_vars <- env_complete[, -i, drop = FALSE]
  if (ncol(other_vars) > 0 && nrow(other_vars) > ncol(other_vars)) {
    target_var <- env_complete[, i]
    vif_model <- lm(target_var ~ ., data = as.data.frame(other_vars))
    r_squared <- summary(vif_model)$r.squared
    vif_value <- 1 / (1 - r_squared)
    vif_results <- rbind(vif_results, data.frame(Variable = colnames(env_complete)[i], VIF = vif_value))
  }
}

#print(vif_results)

# Convert environmental variable names for better visualization in VIF plot
vif_results$Variable_Label <- case_when(
  vif_results$Variable == "salinity" ~ "Salinity",
  vif_results$Variable == "temperature" ~ "Temperature", 
  vif_results$Variable == "pH" ~ "pH",
  vif_results$Variable == "dissolved_o2" ~ "Dissolved O2",
  vif_results$Variable == "sediment_particle_size" ~ "Sediment particle size",
  vif_results$Variable == "suspended_solids" ~ "Suspended solids",
  vif_results$Variable == "Mg" ~ "Magnesium",
  vif_results$Variable == "K" ~ "Potassium",
  vif_results$Variable == "Ca" ~ "Calcium",
  vif_results$Variable == "P" ~ "Phosphorus",
  vif_results$Variable == "Al" ~ "Aluminum",
  vif_results$Variable == "V" ~ "Vanadium",
  vif_results$Variable == "B" ~ "Boron",
  vif_results$Variable == "Fe" ~ "Iron",
  vif_results$Variable == "Cu" ~ "Copper",
  TRUE ~ vif_results$Variable  # Keep original name for any other variables
)

# Barplot of VIF values
vif_plot <- ggplot(vif_results, aes(x = reorder(Variable_Label, VIF), y = VIF)) +
  geom_bar(stat = "identity", fill = "white", colour = "black") +
  coord_flip() +
  labs(title = "Variance Inflation Factor (VIF) for Environmental Variables",
       x = "Environmental Variables", y = "VIF") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),
    axis.title.y = element_blank(),
  ) +
  geom_hline(yintercept = 7, linetype = "dashed", color = "black") +
  annotate("text", x = Inf, y = 7, label = "VIF = 7", vjust = -0.5, hjust = 1.1, color = "red", size = 3)

# Save VIF results to CSV
#write.csv(vif_results, "vif_results.csv", row.names = FALSE)
#ggsave("vif_plot.png", vif_plot, width = 8, height = 6, dpi = 1200)

# Filter variables with VIF < 7
low_vif_vars <- vif_results$Variable[vif_results$VIF < 7]
env_filtered <- env_complete[, low_vif_vars, drop = FALSE]

###############################################################################
################# Transform community data to fourth root #####################
###############################################################################

# Apply fourth root transformation to community data
community_transformed <- community_complete^0.25

###############################################################################
################# Distance-based Redundancy Analysis (dbRDA) #################
###############################################################################

# Calculate distance matrix for community data (Bray-Curtis dissimilarity)
community_dist <- vegdist(community_transformed, method = "bray")

# Perform dbRDA with VIF-filtered environmental variables
dbrda_result <- dbrda(community_dist ~ ., data = as.data.frame(env_filtered))
#print(summary(dbrda_result))

# Test significance of dbRDA axes
dbrda_significance <- anova(dbrda_result, permutations = 9999)
#print(dbrda_significance)

# Test significance of individual environmental variables
dbrda_var_significance <- anova(dbrda_result, by = "terms", permutations = 9999)
#print(dbrda_var_significance)

# Save dbRDA results to CSV
dbrda_summary <- summary(dbrda_result)
#write.csv(as.data.frame(dbrda_summary$cont$importance), "dbrda_eigenvalues.csv", row.names = TRUE)
#write.csv(as.data.frame(dbrda_var_significance), "dbrda_variable_significance.csv", row.names = TRUE)
#write.csv(as.data.frame(dbrda_significance), "dbrda_overall_significance.csv", row.names = TRUE)

# Create dbRDA biplot
# Extract dbRDA scores for plotting
site_scores <- as.data.frame(scores(dbrda_result, display = "sites"))
species_scores <- as.data.frame(scores(dbrda_result, display = "species"))
env_scores <- as.data.frame(scores(dbrda_result, display = "bp"))

# Add season information (needed for plotting)
season_info <- combined_data$season[match(rownames(community_complete), combined_data$sampling_points)]

# Ensure season_info has capitalized values
season_info <- ifelse(season_info == "wet", "Wet",
                     ifelse(season_info == "dry", "Dry", season_info))

# Add labels and season information
site_scores$site_label <- rownames(site_scores)
site_scores$season <- season_info
species_scores$species_label <- rownames(species_scores)
env_scores$env_label <- rownames(env_scores)

# Convert environmental variable names for better visualization
env_scores$env_label <- case_when(
  env_scores$env_label == "salinity" ~ "Salinity",
  env_scores$env_label == "temperature" ~ "Temperature", 
  env_scores$env_label == "pH" ~ "pH",
  env_scores$env_label == "dissolved_o2" ~ "Dissolved O2",
  env_scores$env_label == "sediment_particle_size" ~ "Sediment particle size",
  env_scores$env_label == "suspended_solids" ~ "Suspended solids",
  env_scores$env_label == "Mg" ~ "Magnesium",
  env_scores$env_label == "K" ~ "Potassium",
  env_scores$env_label == "Ca" ~ "Calcium",
  env_scores$env_label == "P" ~ "Phosphorus",
  env_scores$env_label == "Al" ~ "Aluminum",
  env_scores$env_label == "V" ~ "Vanadium",
  env_scores$env_label == "B" ~ "Boron",
  env_scores$env_label == "Fe" ~ "Iron",
  env_scores$env_label == "Cu" ~ "Copper",
  TRUE ~ env_scores$env_label  # Keep original name for any other variables
)

dbrda_plot <- ggplot() +
  geom_segment(data = env_scores, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2), 
               color = "darkgreen", arrow = arrow(length = unit(1, "cm")), 
               linewidth = 1.5, alpha = 0.7) +
  geom_text(data = env_scores, aes(x = dbRDA1, y = dbRDA2, label = env_label), 
            color = "darkgreen", size = 5, fontface = "bold", 
            vjust = -1.0, hjust = 0.5, alpha = 0.7) +
  geom_point(data = site_scores, aes(x = dbRDA1, y = dbRDA2, color = season), 
             size = 5, alpha = 0.7) +
  stat_ellipse(data = site_scores, aes(x = dbRDA1, y = dbRDA2, color = season), 
               level = 0.95, linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_text(data = site_scores, aes(x = dbRDA1, y = dbRDA2, label = site_label), 
            size = 5, vjust = -0.5, hjust = 0.5) +
  scale_color_manual(values = c("Wet" = "#72309e", "Dry" = "#9ec4e3"), 
                     name = "Season") +
  labs(title = "dbRDA Biplot - Sites and Environmental Variables",
       subtitle = paste0("Overall dbRDA significance: p = ", round(dbrda_significance$`Pr(>F)`[1], 3)),
       x = "dbRDA1", y = "dbRDA2") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)

# Save the plot
#ggsave("dbrda_biplot_ggplot.png", dbrda_plot, width = 12, height = 8, dpi = 1200)

################################################################################
################## NMDS Ordination with Environmental Fitting ##################
################################################################################

# Perform NMDS with Bray-Curtis dissimilarity
nmds_result <- metaMDS(community_transformed, distance = "bray", k = 2, trymax = 100)
#print(nmds_result)

# Fit environmental vectors to NMDS ordination
env_fit <- envfit(nmds_result, env_filtered, permutations = 9999)
#print(env_fit)

# Create NMDS plot with environmental vectors
# Extract NMDS scores for plotting
nmds_site_scores <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_species_scores <- as.data.frame(scores(nmds_result, display = "species"))

# Add labels and season information
nmds_site_scores$site_label <- rownames(nmds_site_scores)
nmds_site_scores$season <- season_info
nmds_species_scores$species_label <- rownames(nmds_species_scores)

# Extract environmental vectors for plotting
env_vectors <- as.data.frame(scores(env_fit, display = "vectors"))
env_vectors$env_label <- rownames(env_vectors)

# Convert environmental variable names for better visualization
env_vectors$env_label <- case_when(
  env_vectors$env_label == "salinity" ~ "Salinity",
  env_vectors$env_label == "temperature" ~ "Temperature", 
  env_vectors$env_label == "pH" ~ "pH",
  env_vectors$env_label == "dissolved_o2" ~ "Dissolved O2",
  env_vectors$env_label == "sediment_particle_size" ~ "Sediment particle size",
  env_vectors$env_label == "suspended_solids" ~ "Suspended solids",
  env_vectors$env_label == "Mg" ~ "Magnesium",
  env_vectors$env_label == "K" ~ "Potassium",
  env_vectors$env_label == "Ca" ~ "Calcium",
  env_vectors$env_label == "P" ~ "Phosphorus",
  env_vectors$env_label == "Al" ~ "Aluminum",
  env_vectors$env_label == "V" ~ "Vanadium",
  env_vectors$env_label == "B" ~ "Boron",
  env_vectors$env_label == "Fe" ~ "Iron",
  env_vectors$env_label == "Cu" ~ "Copper",
  TRUE ~ env_vectors$env_label  # Keep original name for any other variables
)

# Create NMDS plot with environmental vectors and PERMANOVA results
# Extract PERMANOVA combined results for subtitle
permanova_r2 <- round(permanova_combined$R2[1], 3)  # Total R² explained
permanova_p <- round(permanova_combined$`Pr(>F)`[1], 3)       # Overall p-value

# Format p-value for display (show full p-value without rounding)
permanova_p_text <- if(permanova_p < 0.001) "p < 0.001" else paste("p =", permanova_p)

nmds_plot <- ggplot() +
  geom_point(data = nmds_site_scores, aes(x = NMDS1, y = NMDS2, color = season),
             size = 4, alpha = 0.8) +
  stat_ellipse(data = nmds_site_scores, aes(x = NMDS1, y = NMDS2, color = season),
               level = 0.95, linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_text(data = nmds_site_scores, aes(x = NMDS1, y = NMDS2, label = site_label),
            size = 5, vjust = -0.5, hjust = 0.5) +
  #geom_segment(data = env_vectors, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
  #             color = "darkgreen", arrow = arrow(length = unit(0.5, "cm")),
  #             linewidth = 1.2, alpha = 0.7) +
  #geom_text(data = env_vectors, aes(x = NMDS1, y = NMDS2, label = env_label),
  #          color = "darkgreen", size = 4, fontface = "bold", 
  #          vjust = -0.5, hjust = 0.5, alpha = 0.7) +
  scale_color_manual(values = c("Wet" = "#72309e", "Dry" = "#9ec4e3"),
                     name = "Season") +
  labs(title = "NMDS Ordination",
       subtitle = paste0("Stress = ", round(nmds_result$stress, 3), 
                        " | PERMANOVA: R² = ", permanova_r2, ", ", permanova_p_text),
       x = "NMDS1", y = "NMDS2") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),
  ) +
  # Add origin lines
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5, color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, color = "gray50")

# Save the plot
#ggsave("nmds_plot_ggplot_PERMANOVA.png", nmds_plot, width = 12, height = 8, dpi = 1200)

######################################################################
######################### PERMANOVA Analysis #########################
######################################################################

# Prepare data for PERMANOVA
season_vector <- season_info
site_vector <- combined_data$site[match(rownames(community_transformed), combined_data$sampling_points)]
environment_vector <- combined_data$environment[match(rownames(community_complete), combined_data$sampling_points)]

# Create distance matrix
community_dist <- vegdist(community_transformed, method = "bray")

# PERMANOVA testing season effects
permanova_season <- adonis2(community_dist ~ season_vector, permutations = 9999)
#print(permanova_season)

# PERMANOVA testing environment type effects
permanova_environment <- adonis2(community_dist ~ environment_vector, permutations = 9999)
#print(permanova_environment)

# Combined PERMANOVA with season and environment
permanova_combined <- adonis2(community_dist ~ season_vector + environment_vector, permutations = 9999)
#print(permanova_combined)

# Test interaction between season and environment
permanova_interaction <- adonis2(community_dist ~ season_vector * environment_vector, permutations = 9999)
#print(permanova_interaction)

# PERMADISPER to test homogeneity of dispersions
betadisper_season <- betadisper(community_dist, season_vector)
permutest_season <- permutest(betadisper_season, permutations = 9999)
#print(permutest_season)

# PERMADISPER to test homogeneity of dispersions
betadisper_environment <- betadisper(community_dist, environment_vector)
permutest_environment <- permutest(betadisper_environment, permutations = 9999)
#print(permutest_environment)

# Save PERMANOVA results to CSV
#write.csv(as.data.frame(permanova_season), "permanova_season_results.csv", row.names = TRUE)
#write.csv(as.data.frame(permanova_environment), "permanova_environmental_results.csv", row.names = TRUE)
#write.csv(as.data.frame(permanova_combined), "permanova_combined_results.csv", row.names = TRUE)
#write.csv(as.data.frame(permanova_interaction), "permanova_interaction_results.csv", row.names = TRUE)

# Save PERMADISPER results to CSV
#write.csv(as.data.frame(permutest_season$tab), "permadisper_season_results.csv", row.names = TRUE)
#write.csv(as.data.frame(permutest_environment$tab), "permadisper_environment_results.csv", row.names = TRUE)

###########################################################################
############################# ANOSIM Analysis #############################
###########################################################################

# ANOSIM testing season effects
anosim_season <- anosim(community_dist, season_vector, permutations = 9999)
#print(anosim_season)

# ANOSIM testing environment effects (Serinhaem, Sorojo, Marau)
anosim_environment <- anosim(community_dist, environment_vector, permutations = 9999)
#print(anosim_environment)

# Save ANOSIM results to CSV
anosim_season_results <- data.frame(
  Analysis = "Season_Effects",
  R_statistic = anosim_season$statistic,
  p_value = anosim_season$signif,
  stringsAsFactors = FALSE
)

anosim_environment_results <- data.frame(
  Analysis = "Environment_Effects", 
  R_statistic = anosim_environment$statistic,
  p_value = anosim_environment$signif,
  stringsAsFactors = FALSE
)

# Combine both ANOSIM results into a single dataframe
anosim_combined_results <- rbind(anosim_season_results, anosim_environment_results)

# Save ANOSIM results
#write.csv(anosim_combined_results, "anosim_combined_results.csv", row.names = FALSE)

##############################################################################
######################### Diversity Indices Analysis #########################
##############################################################################

# Calculate diversity indices
shannon_diversity <- diversity(community_complete, index = "shannon")
simpson_diversity <- diversity(community_complete, index = "simpson")
richness <- specnumber(community_complete)
evenness <- shannon_diversity / log(richness)
# Margalef richness index
margalef_richness <- (richness - 1) / log(nrow(community_complete))

# Create diversity dataframe
diversity_data <- data.frame(
  sampling_points = rownames(community_complete),
  season = season_vector,
  site = site_vector,
  shannon = shannon_diversity,
  simpson = simpson_diversity,
  richness = richness,
  evenness = evenness,
  margalef = margalef_richness
)
#print(summary(diversity_data))

# Calculate the mean and sd value of each index per season
diversity_summary <- diversity_data %>%
  group_by(season) %>%
  summarise(
    mean_shannon = mean(shannon, na.rm = TRUE),
    sd_shannon = sd(shannon, na.rm = TRUE),
    mean_simpson = mean(simpson, na.rm = TRUE),
    sd_simpson = sd(simpson, na.rm = TRUE),
    mean_richness = mean(richness, na.rm = TRUE),
    sd_richness = sd(richness, na.rm = TRUE),
    mean_evenness = mean(evenness, na.rm = TRUE),
    sd_evenness = sd(evenness, na.rm = TRUE),
    mean_margalef = mean(margalef, na.rm = TRUE),
    sd_margalef = sd(margalef, na.rm = TRUE)
  )

# Save diversity indices to CSV
#write.csv(diversity_data, "diversity_indices.csv", row.names = FALSE)

# Test normality of ANOVA residuals for choosing appropriate statistical test
# If residuals are normal (p > 0.05), use ANOVA; if not (p < 0.05), use Kruskal-Wallis

# Shannon diversity
shannon_aov <- aov(shannon ~ season, data = diversity_data)
shapiro_shannon <- shapiro.test(residuals(shannon_aov))
if (shapiro_shannon$p.value > 0.05) {
  shannon_test <- shannon_aov
  shannon_p_value <- summary(shannon_test)[[1]][["Pr(>F)"]][1]
  shannon_test_type <- "ANOVA"
  #print(summary(shannon_test))
} else {
  shannon_test <- kruskal.test(shannon ~ season, data = diversity_data)
  shannon_p_value <- shannon_test$p.value
  shannon_test_type <- "Kruskal-Wallis"
  #print(shannon_test)
}

# Simpson diversity
simpson_aov <- aov(simpson ~ season, data = diversity_data)
shapiro_simpson <- shapiro.test(residuals(simpson_aov))
if (shapiro_simpson$p.value > 0.05) {
  simpson_test <- simpson_aov
  simpson_p_value <- summary(simpson_test)[[1]][["Pr(>F)"]][1]
  simpson_test_type <- "ANOVA"
  #print(summary(simpson_test))
} else {
  simpson_test <- kruskal.test(simpson ~ season, data = diversity_data)
  simpson_p_value <- simpson_test$p.value
  simpson_test_type <- "Kruskal-Wallis"
  #print(simpson_test)
}

# Species richness
richness_aov <- aov(richness ~ season, data = diversity_data)
shapiro_richness <- shapiro.test(residuals(richness_aov))
if (shapiro_richness$p.value > 0.05) {
  richness_test <- richness_aov
  richness_p_value <- summary(richness_test)[[1]][["Pr(>F)"]][1]
  richness_test_type <- "ANOVA"
  #print(summary(richness_test))
} else {
  richness_test <- kruskal.test(richness ~ season, data = diversity_data)
  richness_p_value <- richness_test$p.value
  richness_test_type <- "Kruskal-Wallis"
  #print(richness_test)
}

# Evenness
evenness_aov <- aov(evenness ~ season, data = diversity_data)
shapiro_evenness <- shapiro.test(residuals(evenness_aov))
if (shapiro_evenness$p.value > 0.05) {
  evenness_test <- evenness_aov
  evenness_p_value <- summary(evenness_test)[[1]][["Pr(>F)"]][1]
  evenness_test_type <- "ANOVA"
  #print(summary(evenness_test))
} else {
  evenness_test <- kruskal.test(evenness ~ season, data = diversity_data)
  evenness_p_value <- evenness_test$p.value
  evenness_test_type <- "Kruskal-Wallis"
  #print(evenness_test)
}

# Marlgalef richness
margalef_aov <- aov(margalef ~ season, data = diversity_data)
shapiro_margalef <- shapiro.test(residuals(margalef_aov))
if (shapiro_margalef$p.value > 0.05) {
  margalef_test <- margalef_aov
  margalef_p_value <- summary(margalef_test)[[1]][["Pr(>F)"]][1]
  margalef_test_type <- "ANOVA"
  #print(summary(margalef_test))
} else {
  margalef_test <- kruskal.test(margalef ~ season, data = diversity_data)
  margalef_p_value <- margalef_test$p.value
  margalef_test_type <- "Kruskal-Wallis"
  #print(margalef_test)
}

# Save shapiro diversity results to CSV
shapiro_results <- data.frame(
  Diversity_Index = c("Shannon", "Simpson", "Richness", "Evenness", "Margalef"),
  Shapiro_p_value = c(shapiro_shannon$p.value, shapiro_simpson$p.value, 
                      shapiro_richness$p.value, shapiro_evenness$p.value, 
                      shapiro_margalef$p.value),
  Test_Type = c(shannon_test_type, simpson_test_type, richness_test_type, 
                evenness_test_type, margalef_test_type)
)
#write.csv(shapiro_results, "diversity_shapiro_results.csv", row.names = FALSE)

# Save diversity indices to CSV significance results
diversity_significance <- data.frame(
  Diversity_Index = c("Shannon", "Simpson", "Richness", "Evenness", "Margalef"),
  p_value = c(shannon_p_value, simpson_p_value, richness_p_value, evenness_p_value, margalef_p_value),
  Test_Type = c(shannon_test_type, simpson_test_type, richness_test_type, evenness_test_type, margalef_test_type)
)
#write.csv(diversity_significance, "diversity_significance_results.csv", row.names = FALSE)

# Create diversity plots using ggplot2
# Reshape data for ggplot2
diversity_long <- diversity_data %>%
  select(sampling_points, season, site, shannon, simpson, richness, evenness, margalef) %>%
  pivot_longer(cols = c(shannon, simpson, richness, evenness, margalef), 
               names_to = "diversity_index", 
               values_to = "value") %>%
  mutate(
    diversity_index = factor(diversity_index, 
                           levels = c("shannon", "simpson", "richness", "evenness", "margalef"),
                           labels = c("Shannon Diversity", "Simpson Diversity", 
                                    "Species Richness", "Evenness", "Margalef Richness"))
  )

# Create individual plots for each diversity index with p-values
# Format p-values for display
# Format p-values for display
format_p_value <- function(p_val) {
  if (p_val < 0.001) return("p < 0.001")
  if (p_val < 0.01) return(paste("p =", round(p_val, 3)))
  if (p_val < 0.05) return(paste("p =", round(p_val, 3)))
  return(paste("p =", round(p_val, 3)))
}

# Calculate boxplot statistics with mean as middle line for all diversity indices
shannon_stats <- diversity_data %>%
  group_by(season) %>%
  summarise(
    ymin = min(shannon, na.rm = TRUE),
    lower = quantile(shannon, 0.25, na.rm = TRUE),
    middle = mean(shannon, na.rm = TRUE),   # mean instead of median
    upper = quantile(shannon, 0.75, na.rm = TRUE),
    ymax = max(shannon, na.rm = TRUE)
  )

simpson_stats <- diversity_data %>%
  group_by(season) %>%
  summarise(
    ymin = min(simpson, na.rm = TRUE),
    lower = quantile(simpson, 0.25, na.rm = TRUE),
    middle = mean(simpson, na.rm = TRUE),   # mean instead of median
    upper = quantile(simpson, 0.75, na.rm = TRUE),
    ymax = max(simpson, na.rm = TRUE)
  )

richness_stats <- diversity_data %>%
  group_by(season) %>%
  summarise(
    ymin = min(richness, na.rm = TRUE),
    lower = quantile(richness, 0.25, na.rm = TRUE),
    middle = mean(richness, na.rm = TRUE),   # mean instead of median
    upper = quantile(richness, 0.75, na.rm = TRUE),
    ymax = max(richness, na.rm = TRUE)
  )

evenness_stats <- diversity_data %>%
  group_by(season) %>%
  summarise(
    ymin = min(evenness, na.rm = TRUE),
    lower = quantile(evenness, 0.25, na.rm = TRUE),
    middle = mean(evenness, na.rm = TRUE),   # mean instead of median
    upper = quantile(evenness, 0.75, na.rm = TRUE),
    ymax = max(evenness, na.rm = TRUE)
  )

margalef_stats <- diversity_data %>%
  group_by(season) %>%
  summarise(
    ymin = min(margalef, na.rm = TRUE),
    lower = quantile(margalef, 0.25, na.rm = TRUE),
    middle = mean(margalef, na.rm = TRUE),   # mean instead of median
    upper = quantile(margalef, 0.75, na.rm = TRUE),
    ymax = max(margalef, na.rm = TRUE)
  )

shannon_plot <- ggplot(diversity_data, aes(x = season, fill = season)) +
  geom_boxplot(
    data = shannon_stats,
    aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax),
    stat = "identity", alpha = 0.8, outlier.shape = 16, outlier.size = 2
  ) +
  geom_jitter(aes(y = shannon), width = 0.2, alpha = 0.6, size = 2) +
  scale_fill_manual(values = c("Wet" = "#72309e", "Dry" = "#9ec4e3")) +
  labs(x = "Season", y = "Shannon Index") +
  annotate("text", x = Inf, y = Inf, label = format_p_value(shannon_p_value), 
           hjust = 1.1, vjust = 1.5, size = 5, color = "black") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

simpson_plot <- ggplot(diversity_data, aes(x = season, fill = season)) +
  geom_boxplot(
    data = simpson_stats,
    aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax),
    stat = "identity", alpha = 0.8, outlier.shape = 16, outlier.size = 2
  ) +
  geom_jitter(aes(y = simpson), width = 0.2, alpha = 0.6, size = 2) +
  scale_fill_manual(values = c("Wet" = "#72309e", "Dry" = "#9ec4e3")) +
  labs(x = "Season", y = "Simpson Index") +
  annotate("text", x = Inf, y = Inf, label = format_p_value(simpson_p_value), 
           hjust = 1.1, vjust = 1.5, size = 5, color = "black") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

richness_plot <- ggplot(diversity_data, aes(x = season, fill = season)) +
  geom_boxplot(
    data = richness_stats,
    aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax),
    stat = "identity", alpha = 0.8, outlier.shape = 16, outlier.size = 2
  ) +
  geom_jitter(aes(y = richness), width = 0.2, alpha = 0.6, size = 2) +
  scale_fill_manual(values = c("Wet" = "#72309e", "Dry" = "#9ec4e3")) +
  labs(x = "Season", y = "Number of Species") +
  annotate("text", x = Inf, y = Inf, label = format_p_value(richness_p_value), 
           hjust = 1.1, vjust = 1.5, size = 5, color = "black") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

evenness_plot <- ggplot(diversity_data, aes(x = season, fill = season)) +
  geom_boxplot(
    data = evenness_stats,
    aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax),
    stat = "identity", alpha = 0.8, outlier.shape = 16, outlier.size = 2
  ) +
  geom_jitter(aes(y = evenness), width = 0.2, alpha = 0.6, size = 2) +
  scale_fill_manual(values = c("Wet" = "#72309e", "Dry" = "#9ec4e3")) +
  labs(x = "Season", y = "Pielou's Evenness") +
  annotate("text", x = Inf, y = Inf, label = format_p_value(evenness_p_value), 
           hjust = 1.1, vjust = 1.5, size = 5, color = "black") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

margalef_plot <- ggplot(diversity_data, aes(x = season, fill = season)) +
  geom_boxplot(
    data = margalef_stats,
    aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax),
    stat = "identity", alpha = 0.8, outlier.shape = 16, outlier.size = 2
  ) +
  geom_jitter(aes(y = margalef), width = 0.2, alpha = 0.6, size = 2) +
  scale_fill_manual(values = c("Wet" = "#72309e", "Dry" = "#9ec4e3")) +
  labs(x = "Season", y = "Margalef Richness") +
  annotate("text", x = Inf, y = Inf, label = format_p_value(margalef_p_value), 
           hjust = 1.1, vjust = 1.5, size = 5, color = "black") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

diversity_combined <- plot_grid(shannon_plot, simpson_plot, richness_plot, evenness_plot,
                               ncol = 2, nrow = 2, align = "hv",
                               labels = c("A", "B", "C", "D"),
                               label_size = 16, label_fontface = "bold")

diversity_combined_margalef <- plot_grid(shannon_plot, margalef_plot, richness_plot, evenness_plot,
                               ncol = 2, nrow = 2, align = "hv",
                               labels = c("A", "B", "C", "D"),
                               label_size = 16, label_fontface = "bold")

# Save the combined plot
#ggsave("diversity_boxplots_simpson.png", diversity_combined, width = 12, height = 10, dpi = 1200)
#ggsave("diversity_boxplots_margalef.png", diversity_combined_margalef, width = 12, height = 10, dpi = 1200)

# Save individual plots
#ggsave("shannon_boxplot_ggplot.png", shannon_plot, width = 10, height = 8, dpi = 1200)
#ggsave("simpson_boxplot_ggplot.png", simpson_plot, width = 10, height = 8, dpi = 1200)
#ggsave("richness_boxplot_ggplot.png", richness_plot, width = 10, height = 8, dpi = 1200)
#ggsave("evenness_boxplot_ggplot.png", evenness_plot, width = 10, height = 8, dpi = 1200)
#ggsave("margalef_boxplot_ggplot.png", margalef_plot, width = 10, height = 8, dpi = 1200)

###################################################################
######################### SIMPER Analysis #########################
###################################################################

# SIMPER analysis between seasons
simper_seasons <- simper(community_complete, season_vector, permutations = 9999)
#print(summary(simper_seasons))

# Extract top contributing species for each comparison
simper_summary <- summary(simper_seasons, ordered = TRUE)

# Save SIMPER results to CSV
comparison_name <- names(simper_summary)[1]
top_species <- simper_summary[[1]]
#write.csv(top_species, "simper_top_species.csv", row.names = TRUE)

# Set number of species to display in barplot
n_species <- 10  # Change this number to display more or fewer species

# Prepare data for barplot
top_n <- head(top_species, n_species)
simper_plot_data <- data.frame(
  species = factor(rownames(top_n), levels = rownames(top_n)[order(top_n$average)]),
  contribution = top_n$average,
  cumulative = top_n$cumsum,
  p_value = round(top_n$p, 3)
)

# Create barplot
simper_barplot <- ggplot(simper_plot_data, aes(x = species, y = contribution)) +
  geom_bar(stat = "identity", aes(fill = p_value < 0.05), color = "gray", alpha = 0.8, linewidth = 0.3) +
  scale_fill_manual(values = c("TRUE" = "gray", "FALSE" = "white")) +
  geom_text(aes(label = paste0("p = ", p_value)),
            hjust = 1.2, vjust = 0.5, size = 4, fontface = "bold") +
  coord_flip() +
  labs(title = paste0("Top ", n_species, " Species Contributing to Seasonal Differences"),
       subtitle = "SIMPER Analysis Results",
       x = "Species", 
       y = "Average Contribution (%)") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "italic"),
    axis.text.x = element_text(size = 14, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(color = "white", fill = NA, linewidth = 0.5),
    legend.position = "none"
  ) +
  expand_limits(x = max(simper_plot_data$contribution) * 1.15)

# Save the plot
#ggsave("simper_barplot.png", simper_barplot, width = 10, height = 8, dpi = 1200)

# Create lollipop plot
simper_lollipop <- ggplot(simper_plot_data, aes(x = species, y = contribution)) +
  geom_segment(aes(x = species, xend = species, y = 0, yend = contribution), color = "gray40", linewidth = 3) +
  geom_point(aes(fill = p_value < 0.05), size = 6, color = "black", shape = 21, stroke = 1) +
  scale_fill_manual(values = c("TRUE" = "white", "FALSE" = "black")) +
  geom_text(aes(label = paste0("p = ", p_value)),
            hjust = -0.3, vjust = 0.5, size = 5, fontface = "bold") +
  coord_flip() +
  labs(title = paste0("Top ", n_species, " Species Contributing to Seasonal Differences"),
       subtitle = "SIMPER Analysis Results",
       x = "Species", 
       y = "Average Contribution (%)") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 15, face = "italic"),
    axis.text.x = element_text(size = 14, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(color = "white", fill = NA, linewidth = 0.5),
    legend.position = "none"
  ) +
  expand_limits(y = max(simper_plot_data$contribution) * 1.15)

# Save the plot
#ggsave("simper_lollipop_plot.png", simper_lollipop, width = 10, height = 8, dpi = 1200)

################################################################################
############################# Diversity Modeling ###############################
################################################################################

# Prepare data for model selection
model_data <- data.frame(
  shannon = diversity_data$shannon,
  richness = diversity_data$richness,
  temperature = env_filtered$temperature,
  salinity = env_filtered$salinity,
  dissolved_o2 = env_filtered$dissolved_o2,
  suspended_solids = env_filtered$suspended_solids,
  sediment_particle_size = env_filtered$sediment_particle_size,
  pH = env_complete$pH
)

################################################################################
# Dissolved O2 GAM

# GAM model of Shannon ~ Dissolved O2
shannon_o2_gam <- gam(shannon ~ s(dissolved_o2), data = model_data)

# GAM model of richness ~ Dissolved O2
richness_o2_gam <- gam(richness ~ s(dissolved_o2), data = model_data)

# Plots for Dissolved O2 models
# Extract statistics for Shannon ~ Dissolved O2
shannon_o2_stats <- summary(shannon_o2_gam)
shannon_o2_p <- round(shannon_o2_stats$s.table[1, "p-value"], 3)
shannon_o2_dev <- round(shannon_o2_stats$dev.expl * 100, 1)
shannon_o2_p_text <- ifelse(shannon_o2_p < 0.001, "p < 0.001", paste("p =", shannon_o2_p))
shannon_o2_annotation <- paste0(shannon_o2_p_text, "\nDev. = ", shannon_o2_dev, "%")

shannon_o2_gam_plot <- ggplot(model_data, aes(x = dissolved_o2, y = shannon)) +
  geom_point(size = 3) +
  geom_smooth(method = "gam", colour = "red", linewidth = 3) +
  labs(x = "Dissolved O2 (mg/L)", y = "Shannon Diversity") +
  annotate("text", x = Inf, y = Inf, label = shannon_o2_annotation, 
           hjust = 1.1, vjust = 1.1, size = 5, fontface = "bold") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank())

# Extract statistics for Richness ~ Dissolved O2
richness_o2_stats <- summary(richness_o2_gam)
richness_o2_p <- round(richness_o2_stats$s.table[1, "p-value"], 3)
richness_o2_dev <- round(richness_o2_stats$dev.expl * 100, 1)
richness_o2_p_text <- ifelse(richness_o2_p < 0.001, "p < 0.001", paste("p =", richness_o2_p))
richness_o2_annotation <- paste0(richness_o2_p_text, "\nDev. = ", richness_o2_dev, "%")

richness_o2_gam_plot <- ggplot(model_data, aes(x = dissolved_o2, y = richness)) +
  geom_point(size = 3) +
  geom_smooth(method = "gam", colour = "red", linewidth = 3) +
  labs(x = "Dissolved O2 (mg/L)", y = "Species Richness") +
  annotate("text", x = Inf, y = Inf, label = richness_o2_annotation, 
           hjust = 1.1, vjust = 1.1, size = 5, fontface = "bold") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16, face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank())

# Combine plots
o2_gam_combined_plot <- plot_grid(shannon_o2_gam_plot, richness_o2_gam_plot,
                               ncol = 1, nrow = 2, align = "hv",
                               label_size = 16, label_fontface = "bold")

################################################################################
# Temperature GAM

# GAM model of Shannon ~ Temperature
shannon_temp_gam <- gam(shannon ~ s(temperature), data = model_data)

# GAM model of richness ~ Temperature
richness_temp_gam <- gam(richness ~ s(temperature), data = model_data)

# Plots for Temperature models
# Extract statistics for Shannon ~ Temperature
shannon_temp_stats <- summary(shannon_temp_gam)
shannon_temp_p <- round(shannon_temp_stats$s.table[1, "p-value"], 3)
shannon_temp_dev <- round(shannon_temp_stats$dev.expl * 100, 1)
shannon_temp_p_text <- ifelse(shannon_temp_p < 0.001, "p < 0.001", paste("p =", shannon_temp_p))
shannon_temp_annotation <- paste0(shannon_temp_p_text, "\nDev. = ", shannon_temp_dev, "%")

shannon_temp_gam_plot <- ggplot(model_data, aes(x = temperature, y = shannon)) +
  geom_point(size = 3) +
  geom_smooth(method = "gam", colour = "red", linewidth = 3) +
  labs(x = "Temperature (°C)", y = "Shannon Diversity") +
  annotate("text", x = Inf, y = Inf, label = shannon_temp_annotation, 
           hjust = 1.1, vjust = 1.1, size = 5, fontface = "bold") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank())

# Extract statistics for Richness ~ Temperature
richness_temp_stats <- summary(richness_temp_gam)
richness_temp_p <- round(richness_temp_stats$s.table[1, "p-value"], 3)
richness_temp_dev <- round(richness_temp_stats$dev.expl * 100, 1)
richness_temp_p_text <- ifelse(richness_temp_p < 0.001, "p < 0.001", paste("p =", richness_temp_p))
richness_temp_annotation <- paste0(richness_temp_p_text, "\nDev. = ", richness_temp_dev, "%")

richness_temp_gam_plot <- ggplot(model_data, aes(x = temperature, y = richness)) +
  geom_point(size = 3) +
  geom_smooth(method = "gam", colour = "red", linewidth = 3) +
  labs(x = "Temperature (°C)", y = "Species Richness") +
  annotate("text", x = Inf, y = Inf, label = richness_temp_annotation, 
           hjust = 1.1, vjust = 1.1, size = 5, fontface = "bold") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16, face = "bold"),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank())

# Combine plots
temp_combined_plot <- plot_grid(shannon_temp_gam_plot, richness_temp_gam_plot,
                               ncol = 1, nrow = 2, align = "hv",
                               label_size = 16, label_fontface = "bold")

################################################################################
# Salinity GAM

# GAM model of Shannon ~ Salinity
shannon_salinity_gam <- gam(shannon ~ s(salinity), data = model_data)

# GAM model of Richness ~ Salinity
richness_salinity_gam <- gam(richness ~ s(salinity), data = model_data)

# Plots for Salinity models
# Extract statistics for Shannon ~ Salinity
shannon_salinity_stats <- summary(shannon_salinity_gam)
shannon_salinity_p <- round(shannon_salinity_stats$s.table[1, "p-value"], 3)
shannon_salinity_dev <- round(shannon_salinity_stats$dev.expl * 100, 1)
shannon_salinity_p_text <- ifelse(shannon_salinity_p < 0.001, "p < 0.001", paste("p =", shannon_salinity_p))
shannon_salinity_annotation <- paste0(shannon_salinity_p_text, "\nDev. = ", shannon_salinity_dev, "%")

shannon_salinity_gam_plot <- ggplot(model_data, aes(x = salinity, y = shannon)) +
  geom_point(size = 3) +
  geom_smooth(method = "gam", colour = "red", linewidth = 3) +
  labs(x = "Salinity (ppt)", y = "Shannon Diversity") +
  annotate("text", x = Inf, y = Inf, label = shannon_salinity_annotation, 
           hjust = 1.1, vjust = 1.1, size = 5, fontface = "bold") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16, face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank())

# Extract statistics for Richness ~ Salinity
richness_salinity_stats <- summary(richness_salinity_gam)
richness_salinity_p <- round(richness_salinity_stats$s.table[1, "p-value"], 3)
richness_salinity_dev <- round(richness_salinity_stats$dev.expl * 100, 1)
richness_salinity_p_text <- ifelse(richness_salinity_p < 0.001, "p < 0.001", paste("p =", richness_salinity_p))
richness_salinity_annotation <- paste0(richness_salinity_p_text, "\nDev. = ", richness_salinity_dev, "%")

richness_salinity_gam_plot <- ggplot(model_data, aes(x = salinity, y = richness)) +
  geom_point(size = 3) +
  geom_smooth(method = "gam", colour = "red", linewidth = 3) +
  labs(x = "Salinity (ppt)", y = "Species Richness") +
  annotate("text", x = Inf, y = Inf, label = richness_salinity_annotation, 
           hjust = 1.1, vjust = 1.1, size = 5, fontface = "bold") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16, face = "bold"),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank())

# Combine plots
salinity_combined_plot <- plot_grid(shannon_salinity_gam_plot, richness_salinity_gam_plot,
                               ncol = 1, nrow = 2, align = "hv",
                               label_size = 16, label_fontface = "bold")

################################################################################
# Suspended Solids

# GAM model of Shannon ~ Suspended Solids
shannon_ss_gam <- gam(shannon ~ s(suspended_solids), data = model_data)

# GAM model of Richness ~ Suspended Solids
richness_ss_gam <- gam(richness ~ s(suspended_solids), data = model_data)

# Plots for Suspended Solids models
# Extract statistics for Shannon ~ Suspended Solids
shannon_ss_stats <- summary(shannon_ss_gam)
shannon_ss_p <- round(shannon_ss_stats$s.table[1, "p-value"], 3)
shannon_ss_dev <- round(shannon_ss_stats$dev.expl * 100, 1)
shannon_ss_p_text <- ifelse(shannon_ss_p < 0.001, "p < 0.001", paste("p =", shannon_ss_p))
shannon_ss_annotation <- paste0(shannon_ss_p_text, "\nDev. = ", shannon_ss_dev, "%")

shannon_ss_gam_plot <- ggplot(model_data, aes(x = suspended_solids, y = shannon)) +
  geom_point(size = 3) +
  geom_smooth(method = "gam", colour = "red", linewidth = 3) +
  labs(x = "Suspended Solids (mg/L)", y = "Shannon Diversity") +
  annotate("text", x = Inf, y = Inf, label = shannon_ss_annotation, 
           hjust = 1.1, vjust = 1.1, size = 5, fontface = "bold") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16, face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank())

# Extract statistics for Richness ~ Suspended Solids
richness_ss_stats <- summary(richness_ss_gam)
richness_ss_p <- round(richness_ss_stats$s.table[1, "p-value"], 3)
richness_ss_dev <- round(richness_ss_stats$dev.expl * 100, 1)
richness_ss_p_text <- ifelse(richness_ss_p < 0.001, "p < 0.001", paste("p =", richness_ss_p))
richness_ss_annotation <- paste0(richness_ss_p_text, "\nDev. = ", richness_ss_dev, "%")

richness_ss_gam_plot <- ggplot(model_data, aes(x = suspended_solids, y = richness)) +
  geom_point(size = 3) +
  geom_smooth(method = "gam", colour = "red", linewidth = 3) +
  labs(x = "Suspended Solids (mg/L)", y = "Species Richness") +
  annotate("text", x = Inf, y = Inf, label = richness_ss_annotation, 
           hjust = 1.1, vjust = 1.1, size = 5, fontface = "bold") +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16, face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank())

# Combine plots
ss_combined_plot <- plot_grid(shannon_ss_gam_plot, richness_ss_gam_plot,
                               ncol = 1, nrow = 2, align = "hv",
                               label_size = 16, label_fontface = "bold")

################################################################################
# Sediment Particle Size

# GAM model of Shannon ~ Sediment Particle Size
shannon_sps_gam <- gam(shannon ~ s(sediment_particle_size), data = model_data)

# GAM model of Richness ~ Sediment Particle Size
richness_sps_gam <- gam(richness ~ s(sediment_particle_size), data = model_data)

# Plots for Sediment Particle Size models
# Extract statistics for Shannon ~ Sediment Particle Size
shannon_sps_stats <- summary(shannon_sps_gam)
shannon_sps_p <- round(shannon_sps_stats$s.table[1, "p-value"], 3)
shannon_sps_dev <- round(shannon_sps_stats$dev.expl * 100, 1)
shannon_sps_p_text <- ifelse(shannon_sps_p < 0.001, "p < 0.001", paste("p =", shannon_sps_p))
shannon_sps_annotation <- paste0(shannon_sps_p_text, "\nDev. = ", shannon_sps_dev, "%")

shannon_sps_gam_plot <- ggplot(model_data, aes(x = sediment_particle_size, y = shannon)) +
  geom_point(size = 3) +
  geom_smooth(method = "gam", colour = "red", linewidth = 3) +
  labs(x = "Sediment Particle Size (φ)", y = "Shannon Diversity") +
  annotate("text", x = Inf, y = Inf, label = shannon_sps_annotation, 
           hjust = 1.1, vjust = 1.1, size = 5, fontface = "bold") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16, face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank())

# Extract statistics for Richness ~ Sediment Particle Size
richness_sps_stats <- summary(richness_sps_gam)
richness_sps_p <- round(richness_sps_stats$s.table[1, "p-value"], 3)
richness_sps_dev <- round(richness_sps_stats$dev.expl * 100, 1)
richness_sps_p_text <- ifelse(richness_sps_p < 0.001, "p < 0.001", paste("p =", richness_sps_p))
richness_sps_annotation <- paste0(richness_sps_p_text, "\nDev. = ", richness_sps_dev, "%")

richness_sps_gam_plot <- ggplot(model_data, aes(x = sediment_particle_size, y = richness)) +
  geom_point(size = 3) +
  geom_smooth(method = "gam", colour = "red", linewidth = 3) +
  labs(x = "Sediment Particle Size (φ)", y = "Species Richness") +
  annotate("text", x = Inf, y = Inf, label = richness_sps_annotation, 
           hjust = 1.1, vjust = 1.1, size = 5, fontface = "bold") +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16, face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank())

# Combine plots
sps_combined_plot <- plot_grid(shannon_sps_gam_plot, richness_sps_gam_plot,
                               ncol = 1, nrow = 2, align = "hv",
                               label_size = 16, label_fontface = "bold")

################################################################################
# pH

# GAM model of Shannon ~ pH
shannon_ph_gam <- gam(shannon ~ s(pH), data = model_data)

# GAM model of Richness ~ pH
richness_ph_gam <- gam(richness ~ s(pH), data = model_data)

# Plots for pH models
# Extract statistics for Shannon ~ pH
shannon_ph_stats <- summary(shannon_ph_gam)
shannon_ph_p <- round(shannon_ph_stats$s.table[1, "p-value"], 3)
shannon_ph_dev <- round(shannon_ph_stats$dev.expl * 100, 1)
shannon_ph_p_text <- ifelse(shannon_ph_p < 0.001, "p < 0.001", paste("p =", shannon_ph_p))
shannon_ph_annotation <- paste0(shannon_ph_p_text, "\nDev. = ", shannon_ph_dev, "%")

shannon_ph_gam_plot <- ggplot(model_data, aes(x = pH, y = shannon)) +
  geom_point(size = 3) +
  geom_smooth(method = "gam", colour = "red", linewidth = 3) +
  labs(x = "pH", y = "Shannon Diversity") +
  annotate("text", x = Inf, y = Inf, label = shannon_ph_annotation, 
           hjust = 1.1, vjust = 1.1, size = 5, fontface = "bold") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16, face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank())

# Extract statistics for Richness ~ pH
richness_ph_stats <- summary(richness_ph_gam)
richness_ph_p <- round(richness_ph_stats$s.table[1, "p-value"], 3)
richness_ph_dev <- round(richness_ph_stats$dev.expl * 100, 1)
richness_ph_p_text <- ifelse(richness_ph_p < 0.001, "p < 0.001", paste("p =", richness_ph_p))
richness_ph_annotation <- paste0(richness_ph_p_text, "\nDev. = ", richness_ph_dev, "%")

richness_ph_gam_plot <- ggplot(model_data, aes(x = pH, y = richness)) +
  geom_point(size = 3) +
  geom_smooth(method = "gam", colour = "red", linewidth = 3) +
  labs(x = "pH", y = "Species Richness") +
  annotate("text", x = Inf, y = Inf, label = richness_ph_annotation, 
           hjust = 1.1, vjust = 1.1, size = 5, fontface = "bold") +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16, face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank())

# Combine plots
ph_combined_plot <- plot_grid(shannon_ph_gam_plot, richness_ph_gam_plot,
                               ncol = 1, nrow = 2, align = "hv",
                               label_size = 16, label_fontface = "bold")

################################################################################

# Combine all gam plots into one figure
all_gam_combined_plot <- plot_grid(o2_gam_combined_plot, temp_combined_plot, 
                                   salinity_combined_plot, ss_combined_plot,
                                   sps_combined_plot, ph_combined_plot,
                               ncol = 6, nrow = 1, align = "hv",
                               label_size = 16, label_fontface = "bold")

# Save the combined plot
#ggsave("all_gam_models.png", all_gam_combined_plot, width = 18, height = 10, dpi = 1200)

################################################################################
###################### Comprehensive GAM Results Summary #######################
################################################################################

# Function to extract comprehensive GAM statistics
extract_gam_stats <- function(gam_model, model_name) {
  gam_sum <- summary(gam_model)
  
  # Extract smooth terms statistics
  smooth_terms <- gam_sum$s.table
  
  # Get the first (and usually only) smooth term
  if (nrow(smooth_terms) > 0) {
    edf <- round(smooth_terms[1, "edf"], 3)
    ref_df <- round(smooth_terms[1, "Ref.df"], 3)
    f_value <- round(smooth_terms[1, "F"], 3)
    p_value <- smooth_terms[1, "p-value"]
    
    # Format p-value
    if (p_value < 0.001) {
      p_formatted <- "< 0.001"
    } else {
      p_formatted <- round(p_value, 3)
    }
  } else {
    edf <- NA
    ref_df <- NA
    f_value <- NA
    p_formatted <- NA
  }
  
  # Model fit statistics
  r_squared <- round(gam_sum$r.sq, 3)
  adj_r_squared <- round(gam_sum$r.sq, 3)  # For GAM, this is the same
  deviance_explained <- round(gam_sum$dev.expl * 100, 1)
  n_obs <- nrow(gam_model$model)
  
  # AIC and model comparison
  aic_value <- round(AIC(gam_model), 2)
  
  # Residual statistics
  residual_df <- gam_sum$residual.df
  scale_est <- round(gam_sum$scale, 4)
  
  return(data.frame(
    Model = model_name,
    N = n_obs,
    R_squared = r_squared,
    Dev_Explained_Percent = deviance_explained,
    EDF = edf,
    Ref_DF = ref_df,
    F_value = f_value,
    P_value = p_formatted,
    AIC = aic_value,
    Residual_DF = residual_df,
    Scale_Parameter = scale_est,
    stringsAsFactors = FALSE
  ))
}

# Apply to all GAM models
gam_models_list <- list(
  shannon_o2_gam = "Shannon ~ s(Dissolved O2)",
  richness_o2_gam = "Richness ~ s(Dissolved O2)",
  shannon_temp_gam = "Shannon ~ s(Temperature)",
  richness_temp_gam = "Richness ~ s(Temperature)",
  shannon_salinity_gam = "Shannon ~ s(Salinity)",
  richness_salinity_gam = "Richness ~ s(Salinity)",
  shannon_ss_gam = "Shannon ~ s(Suspended Solids)",
  richness_ss_gam = "Richness ~ s(Suspended Solids)",
  shannon_sps_gam = "Shannon ~ s(Sediment Particle Size)",
  richness_sps_gam = "Richness ~ s(Sediment Particle Size)",
  shannon_ph_gam = "Shannon ~ s(pH)",
  richness_ph_gam = "Richness ~ s(pH)"
)

# Create comprehensive summary table
gam_comprehensive_summary <- data.frame()

for (model_name in names(gam_models_list)) {
  if (exists(model_name)) {
    model_obj <- get(model_name)
    model_label <- gam_models_list[[model_name]]
    stats <- extract_gam_stats(model_obj, model_label)
    gam_comprehensive_summary <- rbind(gam_comprehensive_summary, stats)
  }
}

# Create results formatted for manuscript reporting
manuscript_results <- data.frame(
  Model = gam_comprehensive_summary$Model,
  Statistical_Result = paste0(
    "F(", gam_comprehensive_summary$EDF, ",", gam_comprehensive_summary$Residual_DF, ") = ",
    gam_comprehensive_summary$F_value, ", p ", 
    ifelse(gam_comprehensive_summary$P_value == "< 0.001", "< 0.001", 
           paste("=", gam_comprehensive_summary$P_value))
  ),
  Effect_Size = paste0(
    "R² = ", gam_comprehensive_summary$R_squared,
    " (", gam_comprehensive_summary$Dev_Explained_Percent, "% deviance explained)"
  ),
  Sample_Size = paste0("n = ", gam_comprehensive_summary$N),
  Model_Selection = paste0("AIC = ", gam_comprehensive_summary$AIC)
)

# Identify significant models
significant_models <- gam_comprehensive_summary[gam_comprehensive_summary$P_value != "NS" & 
                                               gam_comprehensive_summary$P_value != "> 0.05" &
                                               !is.na(gam_comprehensive_summary$P_value), ]
# Add R-squared and Deviance Explained columns to significant models
significant_models$R_squared <- gam_comprehensive_summary$R_squared[gam_comprehensive_summary$Model %in% significant_models$Model]
significant_models$Dev_Explained_Percent <- gam_comprehensive_summary$Dev_Explained_Percent[gam_comprehensive_summary$Model %in% significant_models$Model]

# Save comprehensive results
#write.csv(gam_comprehensive_summary, "gam_comprehensive_results.csv", row.names = FALSE)
#write.csv(manuscript_results, "gam_manuscript_results.csv", row.names = FALSE)
#write.csv(significant_models, "gam_significant_models.csv", row.names = FALSE)