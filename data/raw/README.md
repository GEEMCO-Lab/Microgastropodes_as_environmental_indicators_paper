# Raw Data

This folder contains the raw, unprocessed data used in the analysis.

## Data Files

### Expected Files:
- `microgastropods_samples.csv` - Sample collection data including site information, GPS coordinates, sampling dates
- `environmental_variables.csv` - Environmental measurements (water quality, substrate type, etc.)
- `species_abundance.csv` - Species identification and abundance data
- `site_metadata.csv` - Site characteristics and environmental conditions

## Data Collection

*Add information about data collection methods, dates, locations, and any protocols followed.*

## Data Format

All data files should be in CSV format with:
- UTF-8 encoding
- Column headers in the first row
- Missing values represented as NA
- Date format: YYYY-MM-DD

## Notes

- Do not modify files in this folder
- All data processing should be done in R scripts that read from this folder and write to `data/processed/`
- Original data sources and citations should be documented here or in the main README
