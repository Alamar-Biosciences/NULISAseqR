# NULISAseqR 1.5.2 (2026-08-24)

## Changes

### New Features
* **XML-driven QC thresholds** - QC thresholds and criteria can now be defined directly in the panel XML via `<QCThresholds>` (and processing parameters via `<Parameters>`), overriding the hardcoded defaults in `QC.R` on a per-flag basis (issue #670). Flags absent from the XML fall back to their hardcoded defaults. New helper `mergeQCCriteria()` overlays XML thresholds onto the `QC*Criteria()` defaults (overriding only value-like fields; flag set, ordering, and `thresholdNames` stay code-defined). Supports per-matrix sample Detectability keys (`Detectability_PLASMA` → `Detectability.plasma`), `appliesWhen="TAP"` gating, and `operator="none"` informational-only rows (displayed, never flag). A new `forceDefaultQC` argument to `loadNULISAseq()` forces the hardcoded defaults and ignores the XML.
* **loadNULISAseq()** - Added optional AQ outlier-removal controls: `AQ_NC_outlier_removal` and `AQ_IPC_outlier_removal` (both default `FALSE`), plus `AQ_IPC_outlier_threshold` (default `3`) and `AQ_IPC_mad_floor` (default `0.3`). When enabled, a maximum of one NC and/or one IPC/CAL outlier is removed per-target during absolute quantification. Default behavior is unchanged. These options only take effect when AQ is computed via `NULISAseqAQ::applyAQ`; they are ignored in the XML-embedded fallback path.
* **mergeNULISAseq()** - CNS and Inflammation panel version variants (e.g. CNS Disease Panel 220 v1 and v2) can now be merged (#3263). `extract_assay_names()` also now treats any unreported or blank `<Assay>` value as the literal name "NULISAseq" rather than "unknown, compatible with anything."
* **insertCovariatesXML()** - Added opt-in attribute-name sanitization for covariate names written into XML.

### Performance
* **lmNULISAseq() / lmerNULISAseq()** - Hoisted invariant per-target merge out of the differential-expression fit loop, reducing redundant work on each target iteration.

### Changes
* **Detectability QC no longer flags** - The hardcoded defaults for all three Detectability checks (sample-level per-matrix, plate-level, target-level `Target_Detectability`) are now threshold=0 and `operator="none"`: they display but never produce a pass/warning verdict. This is a permanent, non-overridable setting - `mergeQCCriteria()` explicitly excludes Detectability from panel-XML overrides, since many existing panel XML files carry a `<QCThresholds>` entry for it that mirrors the old hardcoded default and would otherwise silently re-enable flagging.

### Bug Fixes
* **writeUpdatedXML()** - No longer crashes with "subscript out of bounds" when the input XML has `ReadCount` nodes for a target or sample barcode that was dropped from `data$targets`/`data$samples` (e.g. a hidden target, `modifiers="hide"`). `.annotate_readcounts()` looked up barcodes with `[[` on a plain named vector, which errors on an unmatched name instead of returning `NA`; it now uses single-bracket indexing so unmatched barcodes are skipped like any other unresolved lookup.
* **mergeQCCriteria()** - Now returns `NULL` immediately when `defaults` is `NULL` (e.g. `QCTargetCriteria(advancedQC=FALSE)` on a plain RQ run), instead of proceeding to warn that every panel-XML `<TargetQC>` threshold key has "no matching computation." Previously, any plain RQ run with panel-XML target thresholds warned on every `QCFlagTarget()` recompute (issue #702).
* **loadNULISAseq()** - The AQ result structure is now consistent between the `NULISAseqAQ` and XML-embedded fallback modes: the fallback path includes the `blank_outlier_table` and `IPC_outlier_table` elements (as `NULL`) that `applyAQ()` returns.
* **loadNULISAseq()** - Outlier-removal arguments are now passed to `NULISAseqAQ::applyAQ` only when the installed version accepts them (filtered against `formals()`), preventing an "unused arguments" error and aborted AQ processing when an older `NULISAseqAQ` is installed.
* **QC Report skeleton (batch effect)** - Replaced regex/substring matching of target and plate names with exact/bounded matching in the batch-effect section (issue #655), fixing incorrect batch-effect warnings and per-plate "Warning/Pass" verdicts. Site A: per-target significance now tests exact membership in each plate's `sig_targets` list instead of `grepl(target, ...)` (which matched substrings like `CCL1`→`CCL14` and regex metacharacters like `LTA|LTB`). Site B: per-plate pairwise fraction now matches plate names with a bounded pattern instead of `grep(plate, ...)` (which matched `Plate_1`→`Plate_10`), while still matching emmeans' parenthesized factor levels.
* **QC Report skeleton (batch effect)** - The pairwise post-hoc test now drops targets whose mixed model failed to fit (`NULL`) before calling `emmeans::emmeans()` (issue #655), which previously aborted the entire report with "Can't handle an object of class NULL".
* **QC Report skeleton (read summary)** - Zero-read / IC-zero sample removal now prunes `IC_normed$normData` alongside `Data`/`samples`/`SampleNames` (issue #655), keeping columns aligned so the IC-normalized control CV% is computed on the correct wells.
* **QC Report skeleton (inter-plate normalization)** - Zero-read / IC-zero sample removal now also prunes the excluded well names from `IPC`/`NC`/`SC`/`Bridge`/`Calibrator` (issue #697). Previously, if the excluded well was itself an IPC/NC/SC/Bridge control (e.g. an IPC well with zero internal-control reads), its name lingered in those lists after its column was dropped from `Data`, causing `interPlateNorm()` to fail with "subscript out of bounds" and aborting QC report generation entirely.
* **interPlateNorm()** - Column subsets used for IPC and intensity normalization now use `drop=FALSE`, and the function now errors clearly ("No IPC wells remain..." / "No samples remain for intensity normalization...") when a plate's IPC list or resolved intensity-normalization sample set (which covers NC and Bridge-based normalization) is empty, rather than silently treating a missing control as "no normalization needed." It also now validates that IPC well names/indices actually exist in the data, raising a specific "not found in the data" error instead of the generic "subscript out of bounds." Previously, a plate left with exactly one control well after exclusion crashed with the opaque `dim(X) must have a positive length`, and a plate left with zero wells of a given type silently skipped that normalization step (e.g. IPC factors defaulted to 1) instead of failing - producing a report that looked normalized but wasn't. This also closes the same hole for `loadNULISAseq(excludeSamples=...)`, which routes through the same function. Also fixed `IPC_method='mean'`, which previously errored on any input due to an argument passed positionally into the wrong parameter.
* **lod()** - Now requires at least 2 NC/blank wells and uses `drop=FALSE` when subsetting them, instead of crashing with `dim(X) must have a positive length` when exactly one blank well remains after sample exclusion, or silently computing a meaningless single-sample LOD.
* **QC Report skeleton (sample QC / Bridge labeling)** - Sample-type relabeling for SC/Bridge wells now checks that the plate actually has wells of that type before matching by name; previously, if a plate's Bridge list became empty (all Bridge wells excluded), `grepl("", ...)` matched every sample name and mislabeled all samples on that plate as "Bridge". Also gave the sample-QC dotchart's label index (`inds`) an explicit default so it can no longer be left undefined when neither the SC nor Bridge relabeling block runs for a plate.
* **readNULISAseq()** - Now builds `aboveLOD` and detectability matrices for csv_long (Counts Report) input, which previously lacked these fields.
* **readNULISAseq()** - Hidden-target reads are now excluded from the RunSummary `ParseableMatch` count, which previously inflated the reported match rate.
* **readNULISAseq()** - Sample matching is now consistent when hide targets are dropped; previously, dropping hidden targets could desynchronize sample column alignment.
* **loadNULISAseq()** - IC (internal control) targets are now kept in sync with target exclusion. Previously, excluding a target that was also the IC could cause a subscript crash. The function now warns on partial IC loss and names the target on full loss.
* **readNULISAseq()** - Numeric IC indices are now resolved to target names before exclusion tracking, preventing misidentification when targets are reordered.
* **sampleBoxplot()** - Guards the IC reference line when the IC target is absent from the normalized matrix, instead of crashing. `filter_run_data()` now uses `drop=FALSE` when subsetting samples/targets to preserve data.frame structure.
* **generate_heatmap()** - Default `row_split` is now gated on `cluster_rows`; previously, passing `cluster_rows=FALSE` with the default `row_split` could produce an error (#669).
* **QCFlagTarget() / QCFlagPlate()** - Guards against dimension-drop when only a single SC or IPC well remains after exclusion, which previously collapsed the matrix to a vector and crashed downstream computations.
* **interCV()** - Now handles the empty replicate-set case (no replicates found) without crashing.
* **insertCovariatesXML()** - Covariate attributes are now written correctly when `sampleName` is plate-suffixed; previously, the match was not anchored to the `_` delimiter, causing mismatches (#644).
* **QC Report** - Now guards against empty or missing template paths before calling `render()`, instead of producing an opaque error (#691).
* **QC Report** - CV Average row now computes correctly, and table columns sort numerically instead of lexicographically (#646). Zero-padded sample IDs are preserved during `render_table` coercion instead of being silently converted to integers.
* **QC Report** - AQ Run QC explanations now use "CAL" terminology consistently instead of "IPC" for calibrator-based absolute quantification.
---

# NULISAseqR 1.5.1 (2026-05-04)

## Changes

### New Features
* **insertCovariatesXML()** - Added `sanitize_names` argument (default `FALSE`). When `TRUE`, covariate column names are passed through `clean_covariate_names(case = "lower_camel")` before being written as XML attributes, matching NAS's covariate-upload sanitation so attribute names stay stable across NAS export → NULISAseqR import → NAS re-import (follow-up to Alamar-Biosciences/NULISA-Analysis-Software#3199). Default `FALSE` preserves the historical write-verbatim contract for direct callers. Errors with a clear message if sanitation would collide two original columns onto the same attribute name (e.g. `sampleName` + `"Sample Name"` both targeting `sampleName`)

### Bug Fixes
* **format_wide_to_long()** - Renamed conflicting user covariate columns that share reserved names (e.g. `PlateID`) by appending `_covar` suffix; emits a message when renaming occurs to inform users of the change
* **detectability_summary()** - Replaced `do.call(rbind)` with `dplyr::bind_rows()` when aggregating samples across plates, fixing crashes when importing plates with different XML schemas (e.g. mismatched metadata columns)
* **mergeNULISAseq()** - Replaced `do.call(rbind)` with `dplyr::bind_rows()` for `RunSummary` aggregation; drops `NA`-named slots from each plate's `RunSummary` list to eliminate spurious `NA.` columns when plates have no balancers in the XML
* **QC Report skeleton** - Fixed sample boxplot axis ticks and NPQ axis label; corrected instrument subsetting for plate effect test when some plates are excluded; removed ANOVA instrument effect test section from internal report output
* **processXML()** - Used `is.na()` instead of `is.null()` for XML attribute checks, fixing silent failures when attributes return `NA` rather than `NULL`
* **insertCovariatesXML()** - Fixed silent skip of barcodes whose `covariates$sampleName` is plate-suffixed (e.g. `LowQC_Plasma_begin_<plate>.xml`), which previously produced `NA` for the injected covariate on every QC sample on re-import (issue Alamar-Biosciences/NULISA-Analysis-Software#3199). Match now anchors on either exact equality or the `_<plate>.xml` suffix used by NAS, avoiding both the prior silent-skip and the substring over-match it could regress into. Empty or whitespace-only barcode text is now skipped explicitly

### Documentation
* **User Guide** - Updated high-abundance and rare case target description: added information regarding the new Neuro 220 panel, clarified language for high abundance and `noDetectability` target behavior

### Infrastructure
* **DESCRIPTION** - Added `XML` and `fields` to Imports, enabling automatic installation of dependencies

---

# NULISAseqR 1.5.0 (2026-03-14)

## Changes

### New Features
* **loadNULISAseq()** - Now accepts a pre-built list structure in addition to file paths, enabling reprocessing of data with sample exclusions without re-parsing XML files
* **get_reverse_curve_targets()** - New exported helper to identify reverse curve targets (Curve_Quant starting with "R")
* **get_noDetectability_targets()** - New exported helper to identify targets with the XML `noDetectability` modifier

### Enhancements
* **Reverse curve & noDetectability target handling** - Reverse curve targets are now fully excluded from detectability and labeled "High Abundance"; rare-case targets with XML `noDetectability` modifier have individual detectability computed but are excluded from summary statistics (mean, sd, median, min, max, # detectable targets)
* **detectability_summary()** - Added `exclude_targets` parameter; consolidated and centralized detectability formatting logic to avoid duplication; "High Abundance" label now shown only for PLASMA/SERUM matrix types; detectability set to NA for non-plasma/serum sample types; returns numeric columns by default (`format=FALSE`) to preserve downstream computation
* **writeNULISAseq()** - IC target(s) now placed at the bottom rows of the RQ data sheet when `include_IC_counts = TRUE`
* **Target name sorting** - Applied case-insensitive sorting (`tolower`) in `quantifiability()` and the QC report skeleton to ensure consistent ordering across platforms
* **Batch effect QC** - Revised batch effect messaging; non-RC `noDetectability` targets kept in batch effect assessment; guarded against edge cases
* **DESCRIPTION** - Minimum R version now declared (required for native pipe usage in `lmNULISAseq.R`); removed `LazyData` field; added `withr` to Suggests

### Bug Fixes
* **Batch effect table** - Fixed crash on pagination and PCA legend truncation in QC report
* **detectability_summary()** - Fixed `apply()` dimension drop in detectability output table; fixed `rowSums` NA handling and guarded against empty target sets in aggregation
* **QCFlagTarget** - Fixed incorrect exclusion of non-RC `noDetectability` targets from detectability calculations
* **Failed_Targets Run QC** - Non-RC `noDetectability` targets now correctly included in Failed_Targets for CV criterion only
* **Well position** - Corrected zero-padding for well position values
* **loadNULISAseq()** - Fixed handling of AQ projects with list input; properly recreates `numericCovariates` for list inputs
* **Namespace fixes** - Added explicit `tibble::` prefix for `column_to_rownames()`

### Testing
* **test-loadNULISAseq.R** - New tests for list input support, AQ project equivalence, and coverage of all output fields
* **test-reverse-curve-detectability.R** - New test suite for reverse curve and `noDetectability` target handling
* **test-writeNULISAseq.R** - Expanded to test entire sheets for both RQ and AQ output; added case-insensitive formatting test; improved robustness to platform differences in string handling and floating point arithmetic

### Infrastructure
* Removed vignettes folder from build
* Added Neuro220 XML files for testing
* Updated CI workflows for hybrid branch pattern

---

# NULISAseqR 1.4.2 (2026-02-15)

## Changes

### Enhancements
* **render_QC_report()** - Improved function parameter ordering and defaults:
  - `xml_files` parameter moved to first position for more intuitive usage
  - Added default values for `output_filename` ("NULISAseq_QC_Report.html") and `output_dir` (current working directory)
  - Added default value for `dataDir` (current working directory)
  - Simplified `Rmd_input_file` path construction using `system.file()`
* **lod()** - Enhanced documentation and parameter handling:
  - Improved parameter ordering (moved `data_matrix` before `blanks`)
  - Enhanced roxygen documentation with clearer return value descriptions
  - Added filtering to ensure `targetNoOutlierDetection` only includes targets present in `data_matrix`

### Bug Fixes
* **quantifiability()** - Fixed sample subsetting issue that could cause errors when sample lists don't match between AQ data and sample information:
  - Now uses `intersect()` to find common samples between `Data_AQ_aM` and `SampleNames`
  - Correctly calculates sample counts for overall and subgroup quantifiability 
  - Prevents errors when processing data with mismatched sample lists
* **loadNULISAseq()** - Added calculation of `LOD_pgmL` (limit of detection in pg/mL units) from XML data for AQ assays
* **targetBoxplot()** - Fixed parameter naming in `lod()` function call to use `data_matrix=` explicitly

### Testing
* **New comprehensive test suites** added to ensure code quality and reliability:
  - `test-importNULISAseq.R` - Tests for `importNULISAseq()` function with and without NULISAseqAQ package, including fallback mode validation and AQ data consistency checks
  - `test-reverse-curve.R` - Tests for reverse curve target handling, including correlation validation, data transformation verification, and NPQ value consistency between `loadNULISAseq()` and `importNULISAseq()`
  - `test-writeNULISAseq.R` - Tests for Excel output generation with both RQ-only and AQ data, including validation of sheet structure, column names, and specific data values
* **Test infrastructure improvements**:
  - Moved test fixtures from `inst/rmarkdown/templates/nulisaseq/skeleton/` to `tests/testthat/fixtures/` for better organization
  - Removed unnecessary `.gitignore` file from skeleton template directory

---

# NULISAseqR 1.4.1 (2026-01-16)

## Changes

### Enhancements
* **importNULISAseq()** - Improved robustness and flexibility for handling NULISAseq data files:
  - Added validation of internal `AUTO_PLATE` IDs with duplicate detection before processing
  - Enhanced parameter mapping for `excludeSamples`, `excludeTargets`, and control parameters (`IC`, `IPC`, `SC`, `NC`, `Bridge`, `Calibrator`) using prioritized keys (user-provided plate names, internal IDs, or fallback names)
  - Improved error handling with clear messages when duplicate plate IDs are detected with named exclusions
* **get_internal_plate_id()** - New utility function to extract `AUTO_PLATE` ID from NULISAseq XML file headers

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
