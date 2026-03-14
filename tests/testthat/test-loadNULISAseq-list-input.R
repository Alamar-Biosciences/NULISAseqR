# Tests for loadNULISAseq list input support (PR #573)
#
# These tests verify that loadNULISAseq can accept a pre-built list structure
# instead of a file path, enabling NAS to reprocess data with sample exclusions
# without re-parsing XML files.
#
# Related: GitHub Issue #1424 - Allow IPC/CAL exclusion and reprocessing

# =============================================================================
# Test 1: Basic list input acceptance
# =============================================================================

test_that("loadNULISAseq accepts a pre-built list structure from readNULISAseq output", {
  # Load from XML first to get canonical structure
  input_file <- test_path("fixtures", "detectability_P1_Tr03.xml")
  raw_from_file <- readNULISAseq(input_file, IPC = NULL, IC = 'mCherry', SC = NULL)

  # Add required fields that readNULISAseq includes
  raw_from_file$xmlFile <- basename(input_file)

  # Now pass the raw structure directly to loadNULISAseq
  result <- loadNULISAseq(raw_from_file, IC = 'mCherry')


  # Verify processing completed successfully
  expect_true(!is.null(result$qcPlate), "qcPlate should be present")
  expect_true(!is.null(result$qcSample), "qcSample should be present")
  expect_true(!is.null(result$NPQ), "NPQ should be present")
  expect_true(nrow(result$qcPlate) >= 5, "qcPlate should have at least 5 rows")
})


# =============================================================================
# Test 2: List input produces equivalent results to file input
# =============================================================================

test_that("loadNULISAseq with list input produces equivalent results to file input", {
  input_file <- test_path("fixtures", "detectability_P1_Tr03.xml")

  # Process via file path (original method)
  result_from_file <- loadNULISAseq(input_file, IPC = NULL, IC = 'mCherry', SC = NULL)

  # Process via list structure (new method)
  raw_structure <- readNULISAseq(input_file, IPC = NULL, IC = 'mCherry', SC = NULL)
  raw_structure$xmlFile <- basename(input_file)
  result_from_list <- loadNULISAseq(raw_structure, IC = 'mCherry')

  # --- Core data matrices (must match exactly) ---
  expect_equal(result_from_file$Data, result_from_list$Data,
               label = "Raw Data matrix should match")
  expect_equal(result_from_file$NPQ, result_from_list$NPQ,
               tolerance = 1e-10, label = "NPQ matrix should match")

  # --- Normalized data matrices ---
  expect_equal(result_from_file$IC_normed$normData, result_from_list$IC_normed$normData,
               tolerance = 1e-10, label = "IC_normed$normData should match")
  expect_equal(result_from_file$normed$interNormData, result_from_list$normed$interNormData,
               tolerance = 1e-10, label = "normed$interNormData should match")
  expect_equal(result_from_file$normed$log2_interNormData, result_from_list$normed$log2_interNormData,
               tolerance = 1e-10, label = "normed$log2_interNormData should match")

  # --- Sample/target metadata ---
  expect_equal(result_from_file$samples, result_from_list$samples,
               label = "samples dataframe should match")
  expect_equal(result_from_file$targets, result_from_list$targets,
               label = "targets dataframe should match")

  # --- Control sample vectors ---
  expect_equal(result_from_file$IC, result_from_list$IC,
               label = "IC should match")
  expect_equal(result_from_file$IPC, result_from_list$IPC,
               label = "IPC should match")
  expect_equal(result_from_file$NC, result_from_list$NC,
               label = "NC should match")
  expect_equal(result_from_file$SC, result_from_list$SC,
               label = "SC should match")
  expect_equal(result_from_file$SampleNames, result_from_list$SampleNames,
               label = "SampleNames should match")

  # --- QC metrics (critical for NAS) ---
  expect_equal(result_from_file$qcPlate, result_from_list$qcPlate,
               tolerance = 1e-10, label = "qcPlate should match")
  expect_equal(result_from_file$qcSample, result_from_list$qcSample,
               tolerance = 1e-10, label = "qcSample should match")

  # --- LOD and detectability ---
  expect_equal(result_from_file$lod$LOD, result_from_list$lod$LOD,
               tolerance = 1e-10, label = "LOD values should match")
  expect_equal(result_from_file$detectability, result_from_list$detectability,
               tolerance = 1e-10, label = "detectability should match")

  # --- Match matrix (sample-target matching) ---
  expect_equal(result_from_file$match_matrix, result_from_list$match_matrix,
               label = "match_matrix should match")
})


# =============================================================================
# Test 3: Validation - missing required fields
# =============================================================================

test_that("loadNULISAseq throws error when list input is missing required fields", {
  # Create minimal incomplete structure
  incomplete_structure <- list(
    Data = matrix(1:10, nrow = 2),
    samples = data.frame(sampleName = c("S1", "S2")),
    # Missing: targets, IC, IPC, NC, SampleNames
    plateID = "test"
  )

  expect_error(
    loadNULISAseq(incomplete_structure, IC = 'mCherry'),
    regexp = "Pre-built input missing required fields.*targets.*IC.*IPC.*NC.*SampleNames",
    label = "Should error on missing required fields"
  )
})


test_that("loadNULISAseq lists specific missing fields in error message", {
  # Missing only NC and SampleNames
  partial_structure <- list(
    Data = matrix(1:10, nrow = 2),
    samples = data.frame(sampleName = c("S1", "S2")),
    targets = data.frame(targetName = c("T1", "T2")),
    IC = "mCherry",
    IPC = c("IPC_1")
    # Missing: NC, SampleNames
  )

  expect_error(
    loadNULISAseq(partial_structure, IC = 'mCherry'),
    regexp = "NC.*SampleNames|SampleNames.*NC",
    label = "Error should list specific missing fields"
  )
})


# =============================================================================
# Test 4: Validation - invalid input types
# =============================================================================

test_that("loadNULISAseq throws error for invalid input types", {
  expect_error(
    loadNULISAseq(123, IC = 'mCherry'),
    regexp = "file must be either a file path.*or a pre-built raw structure",
    label = "Should reject numeric input"
  )

  expect_error(
    loadNULISAseq(NULL, IC = 'mCherry'),
    regexp = "file must be either a file path.*or a pre-built raw structure",
    label = "Should reject NULL input"
  )
})


# =============================================================================
# Test 5: numericCovariates recreation
# =============================================================================

test_that("loadNULISAseq recreates numericCovariates when missing from list input", {
  input_file <- test_path("fixtures", "detectability_P1_Tr03.xml")
  raw_structure <- readNULISAseq(input_file, IPC = NULL, IC = 'mCherry', SC = NULL)

  # Remove numericCovariates to test recreation
  raw_structure$numericCovariates <- NULL
  raw_structure$xmlFile <- basename(input_file)

  result <- loadNULISAseq(raw_structure, IC = 'mCherry')

  # Processing should complete without error
  expect_true(!is.null(result$qcPlate), "Processing should complete")

  # numericCovariates should be recreated
  expect_true(!is.null(result$numericCovariates),
              info = "numericCovariates should be present after recreation")
  expect_true(is.logical(result$numericCovariates),
              info = "numericCovariates should be a logical vector")
  expect_true(!is.null(names(result$numericCovariates)),
              info = "numericCovariates should be a named vector")

  # Protected columns should be marked as non-numeric in numericCovariates
  protected_cols <- c("sampleName", "sampleType", "plateID")
  for (col in protected_cols) {
    if (col %in% names(result$samples) && col %in% names(result$numericCovariates)) {
      expect_false(result$numericCovariates[[col]],
                   label = paste("Protected column", col,
                                 "should be marked non-numeric in numericCovariates"))
    }
  }
})


# =============================================================================
# Test 6: Sample exclusion with list input
# =============================================================================

test_that("loadNULISAseq with list input supports excludeSamples parameter", {
  input_file <- test_path("fixtures", "detectability_P1_Tr03.xml")
  raw_structure <- readNULISAseq(input_file, IPC = NULL, IC = 'mCherry', SC = NULL)
  raw_structure$xmlFile <- basename(input_file)

  # Get IPC samples to exclude one
  ipc_samples <- raw_structure$IPC
  skip_if(length(ipc_samples) < 2, "Need at least 2 IPC samples to test exclusion")

  sample_to_exclude <- ipc_samples[1]

  # Process with exclusion
  result <- loadNULISAseq(raw_structure, IC = 'mCherry',
                          excludeSamples = sample_to_exclude)

  # Excluded sample should not be in processed data
  expect_false(sample_to_exclude %in% colnames(result$Data),
               label = "Excluded sample should not be in Data columns")
  expect_false(sample_to_exclude %in% result$IPC,
               label = "Excluded sample should not be in IPC list")
})


# =============================================================================
# Test 7: AQ project handling with list input
# =============================================================================

test_that("loadNULISAseq validates ExecutionDetails$Abs for AQ projects", {
  # Create structure that looks like an AQ project
  input_file <- test_path("fixtures", "detectability_P1_Tr03.xml")
  raw_structure <- readNULISAseq(input_file, IPC = NULL, IC = 'mCherry', SC = NULL)
  raw_structure$xmlFile <- basename(input_file)

  # Add ExecutionDetails$Abs as non-dataframe to trigger validation error
  raw_structure$ExecutionDetails <- list(Abs = "not_a_dataframe")

  expect_error(
    loadNULISAseq(raw_structure, IC = 'mCherry'),
    regexp = "Pre-built input for AQ project must include ExecutionDetails\\$Abs as a data frame",
    label = "Should validate ExecutionDetails$Abs type for AQ projects"
  )
})


test_that("loadNULISAseq processes AQ project with valid ExecutionDetails$Abs", {
  # Skip if no AQ fixture available
  aq_file <- test_path("fixtures", "XML_v1.3.0_with_AQ.xml")
  skip_if_not(file.exists(aq_file), "AQ fixture file not available")

  raw_structure <- readNULISAseq(aq_file)
  raw_structure$xmlFile <- basename(aq_file)

  # Only test if this is actually an AQ project
  has_aq <- !is.null(raw_structure$ExecutionDetails$Abs)
  skip_if_not(has_aq, "Fixture is not an AQ project")

  # Should process without error
  result <- loadNULISAseq(raw_structure)
  expect_true(!is.null(result$qcPlate), "AQ project should process successfully")
})


test_that("loadNULISAseq with list input produces equivalent AQ results to file input", {
  # Use the dedicated AQ fixture
  aq_file <- test_path("fixtures", "Analysis_INF250_Lot4_AQ_LC_R3_20241229.xml")
  skip_if_not(file.exists(aq_file), "AQ fixture file not available")

  # Process via file path (original method)
  result_from_file <- loadNULISAseq(aq_file)

  # Verify this is actually an AQ project
  skip_if(is.null(result_from_file$AQ), "Fixture did not produce AQ output")

  # Process via list structure (new method)
  raw_structure <- readNULISAseq(aq_file)
  raw_structure$xmlFile <- basename(aq_file)
  result_from_list <- loadNULISAseq(raw_structure)

  # --- Core data matrices ---
  expect_equal(result_from_file$Data, result_from_list$Data,
               label = "AQ: Raw Data matrix should match")
  expect_equal(result_from_file$NPQ, result_from_list$NPQ,
               tolerance = 1e-10, label = "AQ: NPQ matrix should match")

  # --- AQ-specific outputs (critical for AQ projects) ---
  expect_equal(result_from_file$AQ, result_from_list$AQ,
               tolerance = 1e-10, label = "AQ: AQ matrix should match")
  expect_equal(result_from_file$AbsAssay, result_from_list$AbsAssay,
               label = "AQ: AbsAssay flag should match")

  # --- Normalized data ---
  expect_equal(result_from_file$IC_normed$normData, result_from_list$IC_normed$normData,
               tolerance = 1e-10, label = "AQ: IC_normed$normData should match")
  expect_equal(result_from_file$normed$interNormData, result_from_list$normed$interNormData,
               tolerance = 1e-10, label = "AQ: normed$interNormData should match")

  # --- Sample/target metadata ---
  expect_equal(result_from_file$samples, result_from_list$samples,
               label = "AQ: samples dataframe should match")
  expect_equal(result_from_file$targets, result_from_list$targets,
               label = "AQ: targets dataframe should match")

  # --- QC metrics ---
  expect_equal(result_from_file$qcPlate, result_from_list$qcPlate,
               tolerance = 1e-10, label = "AQ: qcPlate should match")
  expect_equal(result_from_file$qcSample, result_from_list$qcSample,
               tolerance = 1e-10, label = "AQ: qcSample should match")

  # --- LOD and detectability ---
  expect_equal(result_from_file$lod, result_from_list$lod,
               tolerance = 1e-10, label = "AQ: LOD should match")
  expect_equal(result_from_file$detectability, result_from_list$detectability,
               tolerance = 1e-10, label = "AQ: detectability should match")

  # --- Control vectors ---
  expect_equal(result_from_file$IPC, result_from_list$IPC,
               label = "AQ: IPC should match")
  expect_equal(result_from_file$NC, result_from_list$NC,
               label = "AQ: NC should match")
  expect_equal(result_from_file$Calibrator, result_from_list$Calibrator,
               label = "AQ: Calibrator should match")
})


# =============================================================================
# Test 8: xmlFile default naming
# =============================================================================

test_that("loadNULISAseq sets default xmlFile when not provided in list input", {
  input_file <- test_path("fixtures", "detectability_P1_Tr03.xml")
  raw_structure <- readNULISAseq(input_file, IPC = NULL, IC = 'mCherry', SC = NULL)

  # Ensure xmlFile is not set
  raw_structure$xmlFile <- NULL
  raw_structure$plateID <- "TestPlate_01"

  result <- loadNULISAseq(raw_structure, IC = 'mCherry')

  # xmlFile should default to plateID
  expect_equal(result$xmlFile, "TestPlate_01",
               label = "xmlFile should default to plateID when not provided")
})


test_that("loadNULISAseq uses 'reprocessed' as fallback xmlFile", {
  input_file <- test_path("fixtures", "detectability_P1_Tr03.xml")
  raw_structure <- readNULISAseq(input_file, IPC = NULL, IC = 'mCherry', SC = NULL)

  # Remove both xmlFile and plateID
  raw_structure$xmlFile <- NULL
  raw_structure$plateID <- NULL

  result <- loadNULISAseq(raw_structure, IC = 'mCherry')

  # xmlFile should fall back to "reprocessed"
  expect_equal(result$xmlFile, "reprocessed",
               label = "xmlFile should fall back to 'reprocessed'")
})


# =============================================================================
# Test 9: Full output structure validation
# =============================================================================

test_that("loadNULISAseq with list input returns identical output structure to file input", {
  input_file <- test_path("fixtures", "detectability_P1_Tr03.xml")

  # Process via both methods
  result_from_file <- loadNULISAseq(input_file, IPC = NULL, IC = 'mCherry', SC = NULL)

  raw_structure <- readNULISAseq(input_file, IPC = NULL, IC = 'mCherry', SC = NULL)
  raw_structure$xmlFile <- basename(input_file)
  result_from_list <- loadNULISAseq(raw_structure, IC = 'mCherry')

  # All field names should be identical
  expect_equal(
    sort(names(result_from_file)),
    sort(names(result_from_list)),
    label = "Output should have identical field names"
  )

  # Compare all numeric matrices and vectors field by field
  numeric_fields <- c("Data", "NPQ", "match_matrix", "detectability")
  for (field in numeric_fields) {
    if (!is.null(result_from_file[[field]]) && !is.null(result_from_list[[field]])) {
      expect_equal(
        result_from_file[[field]],
        result_from_list[[field]],
        tolerance = 1e-10,
        label = sprintf("Field '%s' should match", field)
      )
    }
  }

  # Compare all dataframes
  df_fields <- c("samples", "targets", "qcPlate", "qcSample")
  for (field in df_fields) {
    if (!is.null(result_from_file[[field]]) && !is.null(result_from_list[[field]])) {
      expect_equal(
        result_from_file[[field]],
        result_from_list[[field]],
        tolerance = 1e-10,
        label = sprintf("Dataframe '%s' should match", field)
      )
    }
  }

  # Compare nested list structures (normed, IC_normed, lod)
  if (!is.null(result_from_file$normed) && !is.null(result_from_list$normed)) {
    expect_equal(names(result_from_file$normed), names(result_from_list$normed),
                 label = "normed structure should have same fields")
    for (subfield in names(result_from_file$normed)) {
      expect_equal(
        result_from_file$normed[[subfield]],
        result_from_list$normed[[subfield]],
        tolerance = 1e-10,
        label = sprintf("normed$%s should match", subfield)
      )
    }
  }

  if (!is.null(result_from_file$lod) && !is.null(result_from_list$lod)) {
    expect_equal(
      result_from_file$lod,
      result_from_list$lod,
      tolerance = 1e-10,
      label = "lod structure should match"
    )
  }
})
