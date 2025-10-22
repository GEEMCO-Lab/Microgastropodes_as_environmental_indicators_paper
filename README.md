# Microgastropods as Benthic Environmental Indicators

This repository contains the data analysis and code supporting the paper:  
**“Microgastropods as benthic environmental indicators in remote tropical areas: detecting spatial and seasonal differences in the little-studied Camamu Bay, Bahia, Brazil.”**

The study evaluates **microgastropods (<10 mm)** as **bioindicators** because of the high diversity and ecological variations across **three tropical river basins** in **Camamu Bay, Bahia (Brazil)** during **wet and dry seasons**.  
It explores how environmental variables influence community structure, diversity, and composition in benthic ecosystems of remote tropical estuaries.

---

## Repository Structure

```text
├── data/
│   ├── env_data_microgastropodes.csv           # Environmental variables (pH, salinity, DO, etc.)
│   ├── community_data_microgastropodes.csv     # Species abundance matrix
│
├── results/
│   ├── *.csv                                   # Statistical summaries, ANOVA, PERMANOVA, VIF, GAM outputs
│   ├── *.png                                   # Visualization panels (diversity, environmental, NMDS, GAM)
│
├── microgastropodes_env_analyses_eduhcgalvao.R # Main R analysis script
│
└── README.md                                   # Repository documentation
```

---

## Workflow Overview

The main script **`microgastropodes_env_analyses_eduhcgalvao.R`** performs the following analytical workflow:

1. **Data Import and Cleaning**  
   Loads and merges environmental and community datasets for 15 sampling sites.

2. **Environmental Analyses**  
   Tests spatial and seasonal differences in environmental variables such as:
   - pH, temperature, salinity, dissolved oxygen
   - Suspended solids and sediment particle size

3. **Diversity Analyses**  
   Calculates and compares:
   - Shannon–Wiener diversity
   - Number of species
   - Margalef richness
   - Pielou’s evenness  
   across seasons and river basins.

4. **Community Composition**  
   - Tests assemblage differences using **PERMANOVA** (Bray-Curtis dissimilarity).  
   - Visualizes results via **NMDS** ordination with 95% confidence ellipses.

5. **Environmental–Diversity Relationships**  
   - Fits **Generalized Additive Models (GAMs)** to explore how environmental gradients affect diversity and richness.  
   - Extracts deviance explained and significance values for each variable.

6. **Assemblage–Environment Relationships**  
   - Applies **distance-based redundancy analysis (dbRDA)** to assess environmental influence on community composition after multicollinearity filtering (VIF < 7).

7. **Species Contribution (SIMPER)**  
   - Identifies top 10 taxa contributing most to spatial and seasonal dissimilarities in community composition.

---

## R Packages Used

| Package | Functionality |
|----------|----------------|
| **dplyr**, **tidyr**, **reshape2** | Data cleaning and transformation |
| **ggplot2**, **cowplot** | Visualization and figure arrangement |
| **car** | Variance Inflation Factor (VIF) analysis and statistical tests |
| **vegan** | Ecological and multivariate analyses (PERMANOVA, NMDS, dbRDA, SIMPER, diversity) |
| **mgcv** | Generalized Additive Models (GAMs) |

---

## Reproducibility

All analyses were conducted in **R version 4.3.3**.

To reproduce the workflow:

```r
# Set working directory to repository root
setwd("path/to/Microgastropodes_as_environmental_indicators")

# Run analysis
source("microgastropodes_env_analyses_eduhcgalvao.R")
