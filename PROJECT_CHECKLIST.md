# Project Checklist

Use this checklist to track your progress through the project setup and analysis.

## Initial Setup

- [ ] Clone the repository
- [ ] Open the R project file (.Rproj)
- [ ] Run `source("setup.R")` to install packages
- [ ] Review QUICKSTART.md and README.md
- [ ] Review data templates in `data/raw/`

## Data Preparation

- [ ] Collect raw data
- [ ] Format data according to templates
- [ ] Place data files in `data/raw/`
  - [ ] microgastropods_samples.csv
  - [ ] environmental_variables.csv
  - [ ] species_abundance.csv
  - [ ] site_metadata.csv (optional)
- [ ] Document data sources in README
- [ ] Verify data formats (UTF-8, proper date formats, etc.)

## Running Analyses

- [ ] Run data cleaning: `source("scripts/01_data_cleaning.R")`
- [ ] Review cleaned data in `data/processed/`
- [ ] Run data processing: `source("scripts/02_data_processing.R")`
- [ ] Review processed data files
- [ ] Run statistical analyses: `source("scripts/03_statistical_analysis.R")`
- [ ] Review results in `outputs/tables/`
- [ ] Generate figures: `source("scripts/04_create_figures.R")`
- [ ] Review figures in `outputs/figures/`
- [ ] Or run all at once: `source("scripts/run_all.R")`

## Manuscript Development

- [ ] Update manuscript.Rmd or manuscript.qmd with results
- [ ] Add study area description
- [ ] Add methods details
- [ ] Insert generated figures
- [ ] Insert results tables
- [ ] Write discussion
- [ ] Add references to references.bib
- [ ] Generate manuscript: `rmarkdown::render("manuscript.Rmd")`
- [ ] Review generated document

## Quality Control

- [ ] Check all figures are publication quality
- [ ] Verify all statistical results
- [ ] Ensure reproducibility (run from clean environment)
- [ ] Check for missing citations
- [ ] Proofread manuscript
- [ ] Verify all data sources are documented

## Version Control

- [ ] Commit changes regularly
- [ ] Write descriptive commit messages
- [ ] Push changes to GitHub
- [ ] Update .gitignore if needed
- [ ] Tag release versions

## Publication Preparation

- [ ] Finalize manuscript
- [ ] Prepare supplementary materials
- [ ] Create data dictionary if needed
- [ ] Add DOI for data (e.g., Zenodo)
- [ ] Update CITATION.cff with final details
- [ ] Archive final version

## Sharing and Collaboration

- [ ] Update README with final information
- [ ] Ensure all scripts are well-commented
- [ ] Test reproducibility on clean system
- [ ] Share repository link
- [ ] Respond to feedback and issues

## Notes

Use this space to track any project-specific notes or deviations from the standard workflow:

---
