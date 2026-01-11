# NULISAseqR <a href="https://nulisaseqr.alamarbio.com/reference/index.html"><img src="man/figures/logo.png" align="right" height="138" alt="NULISAseqR logo"/></a>

**NULISAseqR** is an R package designed for comprehensive analysis of proteomic data from the **NULISAseq** platform (Alamar Biosciences). NULISAseq is a high-multiplexed proteomic assay that uses nucleic acid-linked immunosandwich assays to measure hundreds of proteins simultaneously.

## Features
- **Data Import**: Read and parse NULISAseq XML output files
- **Quality Control**: QC metrics and visualizations
- **Normalization**: Intra-plate and inter-plate normalization using internal controls
- **Statistical Analysis**: Differential abundance testing and outcome modeling
- **Visualization**: Heatmaps, PCA plots, volcano plots, and more
- **Reporting**: Generate comprehensive automated HTML QC reports with metrics and figures

## Installation

### System Requirements
- R version 4.4 or higher
- Installation time: < 5 minutes

### macOS Users
**Important**: Before installing NULISAseqR on macOS, you need to install Xcode Command Line Tools:

1. Install Xcode from the App Store (if not already installed)
2. Open Terminal and run:
```bash
sudo xcodebuild -license accept
```

### Windows Users
**Important**: Before installing NULISAseqR on Windows, you need to install Rtools:

1. Go to [https://cran.r-project.org/bin/windows/Rtools/](https://cran.r-project.org/bin/windows/Rtools/)
2. Download and install the appropriate Rtools version for your R installation
3. Follow the installation instructions on the website

### Installation Steps
```r
# 1. Install devtools
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

# 2. Install BiocManager for Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# 3. Install ComplexHeatmap from Bioconductor
BiocManager::install("ComplexHeatmap")

# 4. Install ggalt from CRAN snapshot
install.packages('ggalt', repos = "http://packagemanager.posit.co/cran/2025-08-02")

# 5. Install PCAtools (Alamar fork)
devtools::install_github('Alamar-Biosciences/PCAtools')

# 6. Install NULISAseqR
devtools::install_github('Alamar-Biosciences/NULISAseqR')
```

> **Note**: If you encounter issues installing packages from source, you may need to restart R (Cmd/Ctrl + Shift + F10 in RStudio) between major installation steps.

Load the package:
```r
library(NULISAseqR)
```

## Resources

### Documentation
- [**Function Reference**](https://nulisaseqr.alamarbio.com/reference/) - Complete documentation of all package functions
- [**User Guide**](https://nulisaseqr.alamarbio.com/user-guide/) - Detailed tutorials, workflows, and examples
- [**News**](https://nulisaseqr.alamarbio.com/news/) - Release notes and changelog

### Getting Help
- **Function documentation**: Type `?function_name` in R (e.g., `?importNULISAseq`)
- **GitHub Issues**: [Report bugs or request features](https://github.com/Alamar-Biosciences/NULISAseqR/issues)
