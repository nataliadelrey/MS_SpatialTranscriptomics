# Transcriptome analysis with spatial resolution for the study of sex differences in multiple sclerosis

## Master's Thesis in Bioinformatics

This repository contains all the code and resources used for the spatial transcriptomic analysis of multiple sclerosis (MS) tissue samples and healthy control tissue, focusing on sex differences and using the workflow integrated in the [`Giotto`](https://drieslab.github.io/Giotto_website/) package.

 

## Data

### **Source:**

-   The raw data comes from the study by [Lerma-Martin et al., 2024](https://doi.org/10.1038/s41593-024-01796-z), generated using 10x Genomics Visium technology.

### **Acess:**

-   Download the data from GEO using the provided Bash scripts (`.sh`). These scripts will automatically organize the files into the output format required by Space Ranger.

 

## Bioinformatics Analysis

1.  Run the Bash scripts to download and organize the raw data.

2.  Execute the R scripts in numerical order for data processing, bioinformatics analysis and results visualization.

3.  Scripts generate all figures and plots included in the study.

 

## Requeriments

The analysis was performed in R (version 4.4.2) with the following package versions:

**Table 1.** R packages and versions used in the analysis.

| Package                | Version |
|------------------------|---------|
| `dplyr`                | 1.1.4   |
| `EnhancedVolcano`      | 1.24.0  |
| `ggbeeswarm`           | 0.7.2   |
| `ggplot2`              | 3.5.2   |
| `ggvenn`               | 0.1.10  |
| `Giotto`               | 4.2.1   |
| `Matrix`               | 1.7-3   |
| `patchwork`            | 1.3.0   |
| `scater`               | 1.34.1  |
| `scDblFinder`          | 1.20.2  |
| `SingleCellExperiment` | 1.28.1  |
| `tibble`               | 3.2.1   |
| `tidyr`                | 1.3.1   |
| `zellkonverter`        | 1.16.0  |
