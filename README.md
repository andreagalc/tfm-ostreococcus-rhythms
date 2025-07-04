# tfm-ostreococcus-rhythms

This repository contains the R script used for the analysis of proteomic rhythmicity in *Ostreococcus tauri* under short-day (SD) and long-day (LD) photoperiods, as part of my MSc Thesis in the “Máster Universitario en Análisis de Datos Ómicos y Biología de Sistemas”.

## Project Summary

The aim of the analysis is to identify and characterize rhythmic protein abundance profiles in *Ostreococcus tauri* using time-course SWATH-MS data under seasonal diel cycles. Rhythmicity is assessed using the RAIN algorithm, and clustering is performed via PCA and HCPC to reveal temporal expression patterns.

## Repository contents

- `script_final.R`: Full R script including:
  - Data preprocessing and normalization using NormalyzerDE
  - Missing data imputation
  - Rhythmicity analysis with the RAIN algorithm
  - Visualization: boxplots, PCA, dendrograms, and Venn diagrams

## 📊 Input data (not included here)

To run this script, place the following files in your working directory:
- `report.pg_matrix_SD.tsv`: raw SWATH-MS data for SD conditions
- `report.pg_matrix_LD.tsv`: raw SWATH-MS data for LD conditions

These data are available at:
- **PRIDE**: [PXD046992](https://www.ebi.ac.uk/pride/archive/projects/PXD046992)

## How to run

1. Install required packages:

```r
install.packages(c("FactoMineR", "factoextra", "readr", "NormalyzerDE", "rain", "VennDiagram"))
````
2. Open `script_final.R` in RStudio.

3. Set the working directory where the input `.tsv` files are located.

4. Run the script step by step or all at once.

## Output

The script will generate:
- Normalized protein abundance matrices for SD and LD
- Lists of rhythmic proteins (24h and 12h periods)
- Plots:
  - Boxplots before and after normalization
  - PCA and HCPC cluster dendrograms
  - Barplots of rhythmic vs non-rhythmic proteins
  - Venn diagram comparing LD and SD rhythmic proteins

## Author

**Andrea García Alcaide**  
MSc in Omics Data Analysis and Systems Biology  
University of Seville & Universidad Internacional de Andalucía  

