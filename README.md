# NULISAseqR

**NULISAseqR** is an R package designed for comprehensive analysis of proteomic data from the **NULISAseq** platform (Alamar Biosciences). NULISAseq is a high-multiplexed proteomic assay that uses nucleic acid-linked immunosandwich assays to measure hundreds of proteins simultaneously.

## Features

- **Data Import**: Read and parse NULISAseq XML output files
- **Quality Control**: Automated QC reporting with metrics and visualizations
- **Normalization**: Intra-plate and inter-plate normalization using internal controls
- **Statistical Analysis**: Differential expression testing and predictive modeling
- **Visualization**: Heatmaps, PCA plots, volcano plots, and more
- **Reporting**: Generate comprehensive HTML QC reports

## Installation

### System Requirements

- R version 4.4 or higher
- Installation time: < 5 minutes

### Install from GitHub

Run the following R code (typical installation time < 5 minutes).

```r
install.packages("devtools")
devtools::install_github(
  "Alamar-Biosciences/NULISAseqR",
  ref = "main"
)
```

Load the package:

```r
library(NULISAseqR)
```

## Resources

### Complete User Guide

For detailed tutorials, workflows, and examples, visit the [**NULISAseqR User Guide**](https://vignettes.nulisaseqr.alamarbio.com/) 

### Function Reference

Browse all functions organized by category at [Reference](https://nulisaseqr.alamarbio.com/reference/index.html)

### Getting Help

- **Function documentation**: Type `?function_name` in R (e.g., `?loadNULISAseq`)
- **GitHub Issues**: [Report bugs or request features](https://github.com/Alamar-Biosciences/NULISAseqR/issues)
