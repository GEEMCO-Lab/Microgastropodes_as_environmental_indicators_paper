# R Script: microgastropodes_environmental_main_analysis.R
# Author: Eduardo HC Galvao - eduardohcgalvao@gmail.com
# Version: 3.5 04/10/2025
# Date: 30/09/2025 (dd/mm/yyyy)
# Description: Comprehensive multivariate analysis of microgastropod community 
# structure in relation to environmental variables across wet and dry seasons.
# Performs VIF-based variable selection, canonical correspondence analysis (CCA),
# NMDS ordination with environmental fitting, PERMANOVA tests, diversity indices
# comparison, and SIMPER analysis to identify key species driving seasonal differences.

# Input: environmental_data.csv, environmental_data_average_macronutrients.csv,
# microgastropodes_abundance_matrix.csv
# Output: CCA_biplot.png, NMDS_plot.png, PERMANOVA_results.csv,
# diversity_indices.csv, SIMPER_results.csv

# If one the following packages are not installed run:
# install.packages("package")
library(ggplot2) # For plotting
library(dplyr) # For data manipulation
library(tidyr) # For data reshaping
library(reshape2) # For data manipulation
library(car) # For VIF analysis and additional statistical tests
library(vegan)     # For CCA, PERMANOVA, diversity indices
library(cowplot)   # For combining plots

set.seed(123) # For reproducibility

# Analysis workflow:
# 1. Load and prepare data
# 2. VIF analysis to check multicollinearity of environmental variables
# 3. Canonical Correspondence Analysis (CCA) with VIF filtered variables
# 4. Test significance of CCA axes and environmental variables
# 5. NMDS ordination with environmental fitting
# 6. PERMANOVA to test effects of environmental variables and seasons
# 7. Kruskal-Wallis or ANOVA for diversity indices across seasons
# 8. SIMPER analysis to identify species contributing to differences

###########################################################################
######################### Load and prepare data ###########################
###########################################################################

# Load abiotic and biotic data
env_macro_data <- read.csv("C:/Users/Usuário/Desktop/Profissional/lab_geemco/Microgastropodes_as_environmental_indicators_paper/data/environmental_data_average_macronutrients.csv")
env_var_data <- read.csv("C:/Users/Usuário/Desktop/Profissional/lab_geemco/Microgastropodes_as_environmental_indicators_paper/data/environmental_variables_data.csv")
community_data <- read.csv("C:/Users/Usuário/Desktop/Profissional/lab_geemco/Microgastropodes_as_environmental_indicators_paper/data/microgastropodes_abundance_matrix.csv")

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

# Create individual environmental plots for macronutrients (Mg, K, Ca, P, Al, V, B, Fe, Cu)

# Magnesium plot
mg_data <- env_plot_data %>% filter(variable_label == "Magnesium (mg/L)")
mg_plot <- ggplot(mg_data, aes(x = site, y = value, fill = season)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  scale_fill_manual(values = c("Wet" = "#72309e", "Dry" = "#9ec4e3"), name = "") +
  labs(x = "", y = "Magnesium (mg/L)") +
  theme_classic() +
  #coord_cartesian(ylim = c(50, 120)) +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank()
  )

# Potassium plot
k_data <- env_plot_data %>% filter(variable_label == "Potassium (mg/L)")
k_plot <- ggplot(k_data, aes(x = site, y = value, fill = season)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  scale_fill_manual(values = c("Wet" = "#72309e", "Dry" = "#9ec4e3"), name = "") +
  labs(x = "", y = "Potassium (mg/L)") +
  theme_classic() +
  #coord_cartesian(ylim = c(5, 20)) +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank()
  )

# Calcium plot
ca_data <- env_plot_data %>% filter(variable_label == "Calcium (mg/L)")
ca_plot <- ggplot(ca_data, aes(x = site, y = value, fill = season)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  scale_fill_manual(values = c("Wet" = "#72309e", "Dry" = "#9ec4e3"), name = "") +
  labs(x = "", y = "Calcium (mg/L)") +
  theme_classic() +
  #coord_cartesian(ylim = c(20, 80)) +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank()
  )

# Phosphorus plot
p_data <- env_plot_data %>% filter(variable_label == "Phosphorus (mg/L)")
p_plot <- ggplot(p_data, aes(x = site, y = value, fill = season)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  scale_fill_manual(values = c("Wet" = "#72309e", "Dry" = "#9ec4e3"), name = "") +
  labs(x = "", y = "Phosphorus (mg/L)") +
  theme_classic() +
  #coord_cartesian(ylim = c(0.5, 2)) +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank()
  )

# Aluminum plot
al_data <- env_plot_data %>% filter(variable_label == "Aluminum (mg/L)")
al_plot <- ggplot(al_data, aes(x = site, y = value, fill = season)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  scale_fill_manual(values = c("Wet" = "#72309e", "Dry" = "#9ec4e3"), name = "") +
  labs(x = "", y = "Aluminum (mg/L)") +
  theme_classic() +
  #coord_cartesian(ylim = c(0.5, 2)) +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank()
  )

# Vanadium plot
v_data <- env_plot_data %>% filter(variable_label == "Vanadium (mg/L)")
v_plot <- ggplot(v_data, aes(x = site, y = value, fill = season)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  scale_fill_manual(values = c("Wet" = "#72309e", "Dry" = "#9ec4e3"), name = "") +
  labs(x = "", y = "Vanadium (mg/L)") +
  theme_classic() +
  #coord_cartesian(ylim = c(0.1, 0.5)) +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank()
  )

# Boron plot
b_data <- env_plot_data %>% filter(variable_label == "Boron (mg/L)")
b_plot <- ggplot(b_data, aes(x = site, y = value, fill = season)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  scale_fill_manual(values = c("Wet" = "#72309e", "Dry" = "#9ec4e3"), name = "") +
  labs(x = "", y = "Boron (mg/L)") +
  theme_classic() +
  #coord_cartesian(ylim = c(0.5, 2)) +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank()
  )

# Iron plot
fe_data <- env_plot_data %>% filter(variable_label == "Iron (mg/L)")
fe_plot <- ggplot(fe_data, aes(x = site, y = value, fill = season)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  scale_fill_manual(values = c("Wet" = "#72309e", "Dry" = "#9ec4e3"), name = "") +
  labs(x = "", y = "Iron (mg/L)") +
  theme_classic() +
  #coord_cartesian(ylim = c(0.5, 2)) +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank()
  )

# Copper plot
cu_data <- env_plot_data %>% filter(variable_label == "Copper (mg/L)")
cu_plot <- ggplot(cu_data, aes(x = site, y = value, fill = season)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  scale_fill_manual(values = c("Wet" = "#72309e", "Dry" = "#9ec4e3"), name = "") +
  labs(x = "", y = "Copper (mg/L)") +
  theme_classic() +
  #coord_cartesian(ylim = c(0.1, 0.5)) +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank()
  )

# Combine macronutrient plots with cowplot and add shared legend
macro_plots_grid <- plot_grid(
  mg_plot, k_plot, ca_plot,
  p_plot, al_plot, v_plot,
  b_plot, fe_plot, cu_plot,
  ncol = 3, nrow = 3, align = "hv",
  labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
  label_size = 16, label_fontface = "bold"
)

# Create shared legend
legend_macro <- get_legend(
  mg_plot +
    scale_fill_manual(values = c("Wet" = "#72309e", "Dry" = "#9ec4e3"), name = "") +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 14, face = "bold"),
          legend.key.size = unit(1, "cm"))
)

# Combine plots with legend at bottom center
macro_plots_grid_leg <- plot_grid(
  macro_plots_grid,
  legend_macro,
  ncol = 1,
  rel_heights = c(1, 0.1)
)

# Save the combined macronutrient plot
#ggsave("macronutrients_barplot.png", macro_plots_grid_leg, width = 18, height = 14, dpi = 1200)

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

print(vif_results)

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
vif_plot

# Save VIF results to CSV
#write.csv(vif_results, "vif_results.csv", row.names = FALSE)
#ggsave("vif_plot.png", vif_plot, width = 8, height = 6, dpi = 1200)

# Filter variables with VIF < 7
low_vif_vars <- vif_results$Variable[vif_results$VIF < 7]
env_filtered <- env_complete[, low_vif_vars, drop = FALSE]

print(paste("Variables retained after VIF filtering:", paste(low_vif_vars, collapse = ", ")))
removed_vars <- vif_results$Variable[vif_results$VIF >= 7]
if (length(removed_vars) > 0) {
  print(paste("Variables removed due to high VIF:", paste(removed_vars, collapse = ", ")))
} else {
  print("No variables removed due to high VIF.")
}

###############################################################################
################### Canonical Correspondence Analysis (CCA) ###################
###############################################################################

# Perform CCA with VIF-filtered environmental variables
cca_result <- cca(community_complete ~ ., data = as.data.frame(env_filtered))
print(summary(cca_result))

# Test significance of CCA axes
cca_significance <- anova(cca_result, permutations = 9999)
print(cca_significance)

# Test significance of individual environmental variables
cca_var_significance <- anova(cca_result, by = "terms", permutations = 9999)
print(cca_var_significance)

# Save CCA results to CSV
cca_summary <- summary(cca_result)
#write.csv(as.data.frame(cca_summary$cont$importance), "cca_eigenvalues.csv", row.names = TRUE)
#write.csv(as.data.frame(cca_summary$biplot), "cca_biplot_scores.csv", row.names = TRUE)
#write.csv(as.data.frame(cca_var_significance), "cca_variable_significance.csv", row.names = TRUE)
#write.csv(as.data.frame(cca_significance), "cca_overall_significance.csv", row.names = TRUE)

# Create CCA biplot
# Extract CCA scores for plotting
site_scores <- as.data.frame(scores(cca_result, display = "sites"))
species_scores <- as.data.frame(scores(cca_result, display = "species"))
env_scores <- as.data.frame(scores(cca_result, display = "bp"))

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

cca_plot <- ggplot() +
  geom_segment(data = env_scores, aes(x = 0, y = 0, xend = CCA1, yend = CCA2), 
               color = "darkgreen", arrow = arrow(length = unit(0.5, "cm")), 
               linewidth = 1.2, alpha = 0.7) +
  geom_text(data = env_scores, aes(x = CCA1, y = CCA2, label = env_label), 
            color = "darkgreen", size = 3.5, fontface = "bold", 
            vjust = -1.0, hjust = 0.5, alpha = 0.7) +
  geom_point(data = site_scores, aes(x = CCA1, y = CCA2, color = season), 
             size = 3, alpha = 0.7) +
  geom_text(data = site_scores, aes(x = CCA1, y = CCA2, label = site_label), 
            size = 4, vjust = -0.5, hjust = 0.5) +
  geom_point(data = species_scores, aes(x = CCA1, y = CCA2), 
             color = "#000000ff", size = 1.2, alpha = 0.6) +
  geom_text(data = species_scores, aes(x = CCA1, y = CCA2, label = species_label), 
            color = "#000000ff", size = 3.0, alpha = 0.8, check_overlap = TRUE, vjust = -0.5) +
  scale_color_manual(values = c("Wet" = "#72309e", "Dry" = "#9ec4e3"), 
                     name = "Season") +
  labs(title = "CCA Biplot - Species and Environmental Variables",
       x = "CCA1", y = "CCA2") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)

# Save the plot
#ggsave("cca_biplot_ggplot.png", cca_plot, width = 12, height = 8, dpi = 1200)

################################################################################
################## NMDS Ordination with Environmental Fitting ##################
################################################################################

# Perform NMDS with Bray-Curtis dissimilarity
nmds_result <- metaMDS(community_complete, distance = "bray", k = 2, trymax = 100)
print(nmds_result)

# Fit environmental vectors to NMDS ordination
env_fit <- envfit(nmds_result, env_filtered, permutations = 9999)
print(env_fit)

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
permanova_p <- permanova_combined$`Pr(>F)`[1]       # Overall p-value

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
ggsave("nmds_plot_ggplot_PERMANOVA_combined.png", nmds_plot, width = 12, height = 8, dpi = 1200)

######################################################################
######################### PERMANOVA Analysis #########################
######################################################################

# Prepare data for PERMANOVA
season_vector <- season_info
site_vector <- combined_data$site[match(rownames(community_complete), combined_data$sampling_points)]

# Create distance matrix
community_dist <- vegdist(community_complete, method = "bray")

# PERMANOVA testing season effects
permanova_season <- adonis2(community_dist ~ season_vector, permutations = 9999)
print(permanova_season)

# PERMANOVA testing environmental variables effects
permanova_env <- adonis2(community_dist ~ ., data = as.data.frame(env_filtered), permutations = 9999)
print(permanova_env)

# Combined PERMANOVA with season and environmental variables
combined_env_data <- cbind(as.data.frame(env_filtered), season = season_vector)
permanova_combined <- adonis2(community_dist ~ ., data = combined_env_data, permutations = 9999)
print(permanova_combined)

# Test each environmental variable separately
# Use only the variables that passed VIF filtering (available in combined_env_data)
# Get the variable names from env_filtered (excluding 'season' which we add separately)
vif_filtered_vars <- colnames(combined_env_data)[colnames(combined_env_data) != "season"]

# Create formula dynamically using only available variables
permanova_formula <- as.formula(paste("community_dist ~", paste(vif_filtered_vars, collapse = " + "), "+ season"))

permanova_terms <- adonis2(
  permanova_formula,
  data = combined_env_data,
  permutations = 9999,
  by = "margin"
)
print(permanova_terms)

# PERMADISPER to test homogeneity of dispersions (assumption check)
betadisper_season <- betadisper(community_dist, season_vector)
permutest_betadisper <- permutest(betadisper_season, permutations = 9999)
print(permutest_betadisper)

# Save PERMANOVA results to CSV
#write.csv(as.data.frame(permanova_season), "permanova_season_results.csv", row.names = TRUE)
#write.csv(as.data.frame(permanova_env), "permanova_environmental_results.csv", row.names = TRUE)
#write.csv(as.data.frame(permanova_combined), "permanova_combined_results.csv", row.names = TRUE)
#write.csv(as.data.frame(permanova_terms), "permanova_terms_results.csv", row.names = TRUE)

###########################################################################
############################# ANOSIM Analysis #############################
###########################################################################

# ANOSIM testing season effects
anosim_season <- anosim(community_dist, season_vector, permutations = 9999)
print(anosim_season)

# ANOSIM testing environment effects (Serinhaem, Sorojo, Marau)
environment_vector <- combined_data$environment[match(rownames(community_complete), combined_data$sampling_points)]
anosim_environment <- anosim(community_dist, environment_vector, permutations = 9999)
print(anosim_environment)

# ANOSIM testing environmental variables effects
anosim_env_results <- data.frame(
  Variable = character(),
  R_statistic = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Get VIF-filtered variable names
vif_filtered_vars <- colnames(env_filtered)

for (var in vif_filtered_vars) {
  # Get the environmental variable values
  env_values <- env_filtered[, var]
  
  # Create high/low groups based on median split
  median_val <- median(env_values, na.rm = TRUE)
  env_groups <- ifelse(env_values >= median_val, "High", "Low")
  
  # Run ANOSIM
  anosim_result <- anosim(community_dist, env_groups, permutations = 9999)
  
  # Store results
  anosim_env_results <- rbind(anosim_env_results, 
                             data.frame(Variable = var,
                                       R_statistic = anosim_result$statistic,
                                       p_value = anosim_result$signif))
  
  # Print individual results
  cat("\n", var, ":\n")
  cat("R statistic:", round(anosim_result$statistic, 4), "\n")
  cat("p-value:", anosim_result$signif, "\n")
  cat("Groups: High (>=", round(median_val, 3), ") vs Low (<", round(median_val, 3), ")\n")
}

# Print summary table
print(anosim_env_results)

# Save ANOSIM results to CSV
#write.csv(anosim_env_results, "anosim_environmental_results.csv", row.names = FALSE)

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
print(summary(diversity_data))

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
  print("Shannon Diversity - ANOVA (residuals normally distributed):")
  print(summary(shannon_test))
} else {
  shannon_test <- kruskal.test(shannon ~ season, data = diversity_data)
  shannon_p_value <- shannon_test$p.value
  shannon_test_type <- "Kruskal-Wallis"
  print("Shannon Diversity - Kruskal-Wallis (residuals not normally distributed):")
  print(shannon_test)
}
print(paste("Shannon residuals normality p-value:", round(shapiro_shannon$p.value, 4)))

# Simpson diversity
simpson_aov <- aov(simpson ~ season, data = diversity_data)
shapiro_simpson <- shapiro.test(residuals(simpson_aov))
if (shapiro_simpson$p.value > 0.05) {
  simpson_test <- simpson_aov
  simpson_p_value <- summary(simpson_test)[[1]][["Pr(>F)"]][1]
  simpson_test_type <- "ANOVA"
  print("Simpson Diversity - ANOVA (residuals normally distributed):")
  print(summary(simpson_test))
} else {
  simpson_test <- kruskal.test(simpson ~ season, data = diversity_data)
  simpson_p_value <- simpson_test$p.value
  simpson_test_type <- "Kruskal-Wallis"
  print("Simpson Diversity - Kruskal-Wallis (residuals not normally distributed):")
  print(simpson_test)
}
print(paste("Simpson residuals normality p-value:", round(shapiro_simpson$p.value, 4)))

# Species richness
richness_aov <- aov(richness ~ season, data = diversity_data)
shapiro_richness <- shapiro.test(residuals(richness_aov))
if (shapiro_richness$p.value > 0.05) {
  richness_test <- richness_aov
  richness_p_value <- summary(richness_test)[[1]][["Pr(>F)"]][1]
  richness_test_type <- "ANOVA"
  print("Species Richness - ANOVA (residuals normally distributed):")
  print(summary(richness_test))
} else {
  richness_test <- kruskal.test(richness ~ season, data = diversity_data)
  richness_p_value <- richness_test$p.value
  richness_test_type <- "Kruskal-Wallis"
  print("Species Richness - Kruskal-Wallis (residuals not normally distributed):")
  print(richness_test)
}
print(paste("Richness residuals normality p-value:", round(shapiro_richness$p.value, 4)))

# Evenness
evenness_aov <- aov(evenness ~ season, data = diversity_data)
shapiro_evenness <- shapiro.test(residuals(evenness_aov))
if (shapiro_evenness$p.value > 0.05) {
  evenness_test <- evenness_aov
  evenness_p_value <- summary(evenness_test)[[1]][["Pr(>F)"]][1]
  evenness_test_type <- "ANOVA"
  print("Evenness - ANOVA (residuals normally distributed):")
  print(summary(evenness_test))
} else {
  evenness_test <- kruskal.test(evenness ~ season, data = diversity_data)
  evenness_p_value <- evenness_test$p.value
  evenness_test_type <- "Kruskal-Wallis"
  print("Evenness - Kruskal-Wallis (residuals not normally distributed):")
  print(evenness_test)
}

# Marlgalef richness
margalef_aov <- aov(margalef ~ season, data = diversity_data)
shapiro_margalef <- shapiro.test(residuals(margalef_aov))
if (shapiro_margalef$p.value > 0.05) {
  margalef_test <- margalef_aov
  margalef_p_value <- summary(margalef_test)[[1]][["Pr(>F)"]][1]
  margalef_test_type <- "ANOVA"
  print("Margalef Richness - ANOVA (residuals normally distributed):")
  print(summary(margalef_test))
} else {
  margalef_test <- kruskal.test(margalef ~ season, data = diversity_data)
  margalef_p_value <- margalef_test$p.value
  margalef_test_type <- "Kruskal-Wallis"
  print("Margalef Richness - Kruskal-Wallis (residuals not normally distributed):")
  print(margalef_test)
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
print(summary(simper_seasons))

# Extract top contributing species for each comparison
simper_summary <- summary(simper_seasons, ordered = TRUE)

# Save SIMPER results to CSV
if (length(simper_summary) > 0) {
  # Get the first comparison (usually most important)
  comparison_name <- names(simper_summary)[1]
  top_species <- simper_summary[[1]]
  
  #write.csv(top_species, "simper_top_species.csv", row.names = TRUE)
  print(paste("Top species contributing to differences saved to simper_top_species.csv"))
  
  # Print top 10 species
  print("Top 10 species contributing to seasonal differences:")
  print(head(top_species, 10))
}

# Create SIMPER plots for top contributing species (lollipop charts)
if (exists("top_species") && nrow(top_species) > 0) {
  
  # Function to create lollipop chart
  create_simper_lollipop <- function(n_species) {
    top_n <- head(top_species, n_species)
    
    # Prepare data
    simper_plot_data <- data.frame(
      species = factor(rownames(top_n), levels = rev(rownames(top_n))),
      contribution = top_n$average,
      cumulative = top_n$cumsum,
      p_value = round(top_n$p, 4)  # Add p-values rounded to 4 decimal places
    )
    
    # Create lollipop chart
    lollipop_plot <- ggplot(simper_plot_data, aes(x = species, y = contribution)) +
      geom_segment(aes(x = species, xend = species, y = 0, yend = contribution), 
                   color = "black", linewidth = 1, alpha = 0.8) +
      geom_point(color = "black", size = 4, alpha = 0.9) +
      geom_text(aes(label = paste0("  p = ", p_value)), 
                hjust = -0.1, size = 5, fontface = "bold") +
      coord_flip() +
      labs(title = paste0("Top ", n_species, " Species Contributing to Seasonal Differences"),
           subtitle = "SIMPER Analysis Results",
           x = "Species", 
           y = "Average Contribution (%)") +
      theme_classic() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5, color = "black"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "italic"),
        axis.text.x = element_text(size = 14),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
      ) +
      expand_limits(y = max(simper_plot_data$contribution) * 1.15)
    
    return(lollipop_plot)
  }
  
  # Create lollipop charts for different numbers of species
  # Top 10 species
  if (nrow(top_species) >= 10) {
    simper_lollipop_10 <- create_simper_lollipop(10)
    ggsave("simper_lollipop_top10_with_p_values.png", simper_lollipop_10, width = 10, height = 8, dpi = 1200)
  }
  
  # Top 25 species
  if (nrow(top_species) >= 25) {
    simper_lollipop_25 <- create_simper_lollipop(25)
    ggsave("simper_lollipop_top25_with_p_values.png", simper_lollipop_25, width = 12, height = 14, dpi = 1200)
  }
  
  # Top 50 species
  if (nrow(top_species) >= 50) {
    simper_lollipop_50 <- create_simper_lollipop(50)
    ggsave("simper_lollipop_top50_with_p_values.png", simper_lollipop_50, width = 14, height = 20, dpi = 1200)
  }
}
