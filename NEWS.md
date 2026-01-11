# NULISAseqR 1.4.0 (2026-01-11)

## Overview

Version 1.4.0 represents a major expansion of the NULISAseqR package, introducing new analytical capabilities, enhanced visualization tools, and improved data processing functions.

## New Features

### Documentation & Installation
* Added pkgdown website for package documentation
* Added comprehensive vignette covering data loading, QC, differential expression, visualization, and case studies
* Added MacOS and Windows installation instructions to documentation

### Data Import & Export
* **importNULISAseq()** - New streamlined function to import NULISAseq data from multiple XML files with improved error handling and validation
* **getXMLVersion()** - Retrieve XML version information from NULISAseq files
* Added support for XML v1.3.0 format compatibility in `loadNULISAseq()` which accommodates absolute quantification (AQ) NULISAseq assay panels

### Quality Control
* **render_QC_report()** - Generate automated quality control reports in HTML format
* **detectability_summary()** - Summarize detectability across multiple runs and sample matrix types

### Statistical Analysis
* **permutation_anova()** - Perform permutation-based ANOVA testing for robust statistical inference

#### Single-Protein Prediction Models
Four new functions enable using single-target NPQ as a predictor in covariate-adjusted linear and logistic regression models:

* **lmNULISAseq_predict()** - Predictions for continuous outcomes using linear regression models
* **lmerNULISAseq_predict()** - Predictions for continuous outcomes from linear mixed-effects models for hierarchical data
* **glmNULISAseq_predict()** - Predictions for binary/count outcomes using generalized linear models
* **glmerNULISAseq_predict()** - Predictions for binary/count outcomes using generalized linear mixed-effects models for hierarchical data

### Visualization Suite

#### Heatmaps
* **generate_heatmap()** - Create publication-quality protein abundance heatmaps with ComplexHeatmap integration, supporting clustering, annotations, and custom color schemes
* **QCplateHeatmap()** - Plate-level quality control heatmaps for identifying spatial patterns

#### Sample & Target Visualization
* **sampleBoxplot()** - Boxplots showing sample distributions
* **sampleQCplot()** - Comprehensive sample quality control plots with multiple metrics
* **plot_plateLayout()** - Visual representation of plate layouts for experimental design

#### Dimensionality Reduction
* **generate_pca()** - Principal component analysis with biplot generation and customizable aesthetics

### NULISAseq Absolute Quantification (AQ) Analysis 
* **targetQCplot()** - Target-level QC visualizations for AQ assay performance monitoring
* **CV_AQ()** - Calculate intra-plate and inter-plate coefficient of variation for AQ runs, with automatic handling of values outside the dynamic range
* **CV_AQ_Hist()** - Visualize CV distributions for quality control monitoring
* **quantifiability()** - Calculate quantifiability metrics across multiple runs
* **quantHist()** - Histogram plots of quantifiability distributions
* Added `withinDR` matrix to AQ output for dynamic range filtering


## Enhanced Functions

### Data Import/Export
* **readNULISAseq()** - Improved XML parsing for better compatibility across file versions, enhanced error handling and validation
* **writeNULISAseq()** - Refactored to utilize `importNULISAseq` function
* **writeUpdatedXML()** - Now loads XML internally instead of requiring pre-loaded data

### Visualization
* **volcanoPlot()** - Major enhancements including:
  - Dual plotting mode: plot both unadjusted and FDR-adjusted p-values simultaneously with color coding (light colors for unadjusted significance, darker colors for FDR significance)
  - Flexible p-value input: accepts either single vector or named list with 'unadj' and 'fdr' p-values
  - Fold change thresholds: `upper_log2FC_threshold` and `lower_log2FC_threshold` parameters for labeling targets based on effect size
  - Automatic axis label adjustment based on p-value type (unadjusted vs FDR-adjusted)
  - Enhanced customization options for colors, fonts, and plot dimensions
* **targetBoxplot()** - Better handling of normalized vs. unnormalized data
* Target detectability boxplots margins adjusted to prevent clipping of long target names
* **alamarColorPalette()** - Expanded color palette options


### QC Report Improvements
* Added batch effect QC section with configurable significance thresholds
* Added ICC (Intraclass Correlation Coefficient) scatterplot for batch effect assessment
* Added combined detectability calculation across plates when common sample matrix types exist
* Both unnormalized and normalized sample boxplots now shown on all reports to better enable internal control QC checking
* Improved QC plot formatting and spacing
* Intra-plate CV "Overall" column renamed to "Average" for clarity
* Detectability summary tables now show denominator (total targets per plate)
* Interactive reactable tables allow row sorting by clicking on column names

### Other Improvements
* License updated to GPL-3

## Bug Fixes

### Data Processing
* Fixed sample tag search to only search within Data subnodes, preventing false matches
* Fixed NA handling in `writeNULISAseq` (empty cells vs "NA" string)
* Fixed numeric covariate detection in linear model functions to properly identify NA values stored as strings

### Quality Control
* Fixed subsetting errors when data contains only 1 row or column (added `drop=FALSE`)
* Fixed Sample QC percentage values (now multiplied by 100)
* Fixed Target QC and CV table display issues

### Statistical Models
* Fixed `drop=FALSE` placement bug in predict model functions to prevent errors when dataset contains only 1 row or column

### Visualization
* Fixed special characters displaying as "&#124;" in tables

### Data Management
* Fixed lazy-load database corruption by disabling LazyData
* In QC report, fixed overall detectability weight calculation when sample matrices differ between plates

## Breaking Changes

* Version number updated from 1.2.0 to 1.4.0 (skipping 1.3.0 as standalone release)
* Some function parameters may have changed order or names - please review documentation

## Getting Help

* **Documentation**: [https://nulisaseqr.alamarbio.com](https://nulisaseqr.alamarbio.com)
* **Issues**: [https://github.com/Alamar-Biosciences/NULISAseqR/issues](https://github.com/Alamar-Biosciences/NULISAseqR/issues)
* **Support**: Contact Alamar Biosciences Bioinformatics Team for assistance
* **Vignettes**: Comprehensive guides covering data import, QC, visualization, and statistical analysis

---

**Full Changelog**: [https://github.com/Alamar-Biosciences/NULISAseqR/compare/main...1.4](https://github.com/Alamar-Biosciences/NULISAseqR/compare/main...1.4)
