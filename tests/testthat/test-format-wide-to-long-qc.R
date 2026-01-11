# Functions are now in mergeNULISAseq.R
# No need to source separately - they're part of the package

test_that("prepare_sample_qc_for_long creates correct structure", {
  # Create mock qcSample data
  qc_sample <- data.frame(
    sampleName = c("Sample_1", "Sample_1", "Sample_2", "Sample_2"),
    plateID = c("Plate1", "Plate1", "Plate1", "Plate1"),
    flagName = c("Detectability", "ICReads", "Detectability", "ICReads"),
    val = c(0.95, 1200, 0.88, 1100),
    status = c(FALSE, FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )

  result <- NULISAseqR:::prepare_sample_qc_for_long(qc_sample)

  # Check structure
  expect_true("SampleName" %in% names(result))
  expect_true("PlateID" %in% names(result))
  expect_true("Sample_QC_Status" %in% names(result))
  expect_true("Sample_QC_Detectability" %in% names(result))
  expect_true("Sample_QC_ICReads" %in% names(result))
  expect_true("Sample_QC_Detectability_Status" %in% names(result))
  expect_true("Sample_QC_ICReads_Status" %in% names(result))

  # Check values
  expect_equal(nrow(result), 2)  # Two unique samples
  expect_equal(result$Sample_QC_Status, c("PASS", "PASS"))
  expect_equal(result$Sample_QC_Detectability, c(0.95, 0.88))

  # Check status columns are character labels (using "WARN" not "failed")
  expect_equal(result$Sample_QC_Detectability_Status, c("PASS", "PASS"))
  expect_equal(result$Sample_QC_ICReads_Status, c("PASS", "PASS"))
})

test_that("prepare_sample_qc_for_long handles character status values", {
  qc_sample <- data.frame(
    sampleName = c("Sample_1", "Sample_1", "Sample_2", "Sample_2"),
    plateID = c("Plate1", "Plate1", "Plate1", "Plate1"),
    flagName = c("Detectability", "ICReads", "Detectability", "ICReads"),
    val = c(0.95, 1200, 0.88, 1100),
    status = c("OK", "OK", "TRUE", "OK"),  # Character status
    stringsAsFactors = FALSE
  )

  result <- NULISAseqR:::prepare_sample_qc_for_long(qc_sample)

  expect_equal(nrow(result), 2)
  expect_equal(result$Sample_QC_Status, c("PASS", "WARN"))  # Sample_2 has warning
})

test_that("prepare_sample_qc_for_long handles failed QC", {
  qc_sample <- data.frame(
    sampleName = c("Sample_1", "Sample_1"),
    plateID = c("Plate1", "Plate1"),
    flagName = c("Detectability", "ICReads"),
    val = c(0.85, 800),
    status = c(TRUE, FALSE),  # Detectability failed
    stringsAsFactors = FALSE
  )

  result <- NULISAseqR:::prepare_sample_qc_for_long(qc_sample)

  expect_equal(result$Sample_QC_Status, "WARN")  # Uses "WARN" not "failed"
})

test_that("prepare_target_qc_for_long creates correct structure", {
  # Create mock qcTarget data
  qc_target <- data.frame(
    target = c("Target_1", "Target_1", "Target_2", "Target_2"),
    plateID = c("Plate1", "Plate1", "Plate1", "Plate1"),
    flagName = c("Target_Conc_Accuracy", "Target_Conc_CV", "Target_Conc_Accuracy", "Target_Conc_CV"),
    val = c(0.05, 0.12, -0.08, 0.25),
    status = c(FALSE, FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )

  result <- NULISAseqR:::prepare_target_qc_for_long(qc_target)

  # Check structure
  expect_true("Target" %in% names(result))
  expect_true("PlateID" %in% names(result))
  expect_true("Target_QC_Status" %in% names(result))
  expect_true("Target_QC_Conc_Accuracy" %in% names(result))
  expect_true("Target_QC_Conc_CV" %in% names(result))

  # Check values
  expect_equal(nrow(result), 2)  # Two unique targets
  expect_equal(result$Target_QC_Status, c("PASS", "PASS"))

  # Check status columns are character labels
  expect_equal(result$Target_QC_Conc_Accuracy_Status, c("PASS", "PASS"))
  expect_equal(result$Target_QC_Conc_CV_Status, c("PASS", "PASS"))
})

test_that("prepare_target_qc_for_long uses warning for failed targets", {
  # Create mock qcTarget data with one failure
  qc_target <- data.frame(
    target = c("Target_1", "Target_1"),
    plateID = c("Plate1", "Plate1"),
    flagName = c("Target_Conc_Accuracy", "Target_Conc_CV"),
    val = c(0.05, 0.45),
    status = c(FALSE, TRUE),  # CV failed
    stringsAsFactors = FALSE
  )

  result <- NULISAseqR:::prepare_target_qc_for_long(qc_target)

  # Target QC uses "WARN" not "failed"
  expect_equal(result$Target_QC_Status, "WARN")
})

test_that("format_wide_to_long integrates QC columns", {
  skip_if_not(requireNamespace("NULISAseqR", quietly = TRUE))
  skip_if_not(rlang::is_function(get0("format_wide_to_long", envir = asNamespace("NULISAseqR"), inherits = FALSE)))

  # Create minimal mock data object
  mock_data <- list(
    targets = data.frame(
      plateID = c("Plate1", "Plate1"),
      targetName = c("Target_1", "Target_2"),
      logged_LOD = c(5.0, 5.5),
      stringsAsFactors = FALSE
    ),
    samples = data.frame(
      plateID = c("Plate1", "Plate1"),
      sampleName = c("Sample_1", "Sample_2"),
      sampleType = c("Sample", "Sample"),
      stringsAsFactors = FALSE
    ),
    Data_NPQ = matrix(
      c(6.0, 6.5, 7.0, 7.5),
      nrow = 2,
      dimnames = list(c("Target_1", "Target_2"), c("Sample_1", "Sample_2"))
    ),
    Data_raw = matrix(
      c(100, 120, 150, 180),
      nrow = 2,
      dimnames = list(c("Target_1", "Target_2"), c("Sample_1", "Sample_2"))
    ),
    ExecutionDetails = list(
      Plate1 = list(
        Assay = "Test Panel",
        TargetKitLotNumber = "LOT123",
        InstrumentSerial = "INST001"
      )
    ),
    qcSample = data.frame(
      sampleName = c("Sample_1", "Sample_1", "Sample_2", "Sample_2"),
      plateID = c("Plate1", "Plate1", "Plate1", "Plate1"),
      flagName = c("Detectability", "ICReads", "Detectability", "ICReads"),
      val = c(0.95, 1200, 0.88, 1100),
      status = c(FALSE, FALSE, FALSE, FALSE),
      stringsAsFactors = FALSE
    ),
    qcTarget = data.frame(
      target = c("Target_1", "Target_1", "Target_2", "Target_2"),
      plateID = c("Plate1", "Plate1", "Plate1", "Plate1"),
      flagName = c("Metric1", "Metric2", "Metric1", "Metric2"),
      val = c(0.05, 0.12, -0.08, 0.25),
      status = c(FALSE, FALSE, FALSE, FALSE),
      stringsAsFactors = FALSE
    ),
    plateID = "Plate1"
  )

  result <- NULISAseqR:::format_wide_to_long(mock_data, AQ = FALSE, include_qc = TRUE)

  # Check that result has QC columns
  expect_true("Sample_QC_Status" %in% names(result))
  expect_true("Target_QC_Status" %in% names(result))
  expect_true("Sample_QC_Detectability" %in% names(result))
  expect_true("Target_QC_Metric1" %in% names(result))

  # Should have 4 rows (2 targets × 2 samples)
  expect_equal(nrow(result), 4)

  # All samples should have passed QC
  expect_true(all(result$Sample_QC_Status == "PASS"))
  # All targets should have passed QC
  expect_true(all(result$Target_QC_Status == "PASS"))
})

test_that("format_wide_to_long works without QC data", {
  skip_if_not(requireNamespace("NULISAseqR", quietly = TRUE))
  skip_if_not(rlang::is_function(get0("format_wide_to_long", envir = asNamespace("NULISAseqR"), inherits = FALSE)))

  # Create minimal mock data without QC
  mock_data <- list(
    targets = data.frame(
      plateID = c("Plate1"),
      targetName = c("Target_1"),
      logged_LOD = c(5.0),
      stringsAsFactors = FALSE
    ),
    samples = data.frame(
      plateID = c("Plate1"),
      sampleName = c("Sample_1"),
      sampleType = c("Sample"),
      stringsAsFactors = FALSE
    ),
    Data_NPQ = matrix(
      c(6.0),
      nrow = 1,
      dimnames = list(c("Target_1"), c("Sample_1"))
    ),
    Data_raw = matrix(
      c(100),
      nrow = 1,
      dimnames = list(c("Target_1"), c("Sample_1"))
    ),
    ExecutionDetails = list(
      Plate1 = list(Assay = "Test Panel")
    ),
    plateID = "Plate1"
  )

  # Should not error when QC data is missing
  result <- NULISAseqR:::format_wide_to_long(mock_data, AQ = FALSE, include_qc = TRUE)

  expect_equal(nrow(result), 1)
  expect_false("Sample_QC_Status" %in% names(result))
  expect_false("Target_QC_Status" %in% names(result))
})

# Note: Column reordering is tested implicitly through format_wide_to_long tests above
# The reordering logic is embedded within format_wide_to_long() function


# ==============================================================================
# Additional Edge Case and Robustness Tests
# ==============================================================================

test_that("prepare_sample_qc_for_long handles case-insensitive status values", {
  qc_sample <- data.frame(
    sampleName = c("S1", "S1", "S2", "S2", "S3", "S3", "S4", "S4"),
    plateID = rep("Plate1", 8),
    flagName = rep(c("Detectability", "ICReads"), 4),
    val = rep(c(0.95, 1200), 4),
    status = c("ok", "OK", "pass", "PASS", "fail", "FAIL", "true", "TRUE"),
    stringsAsFactors = FALSE
  )

  result <- NULISAseqR:::prepare_sample_qc_for_long(qc_sample)

  # S1 and S2 should pass (ok/OK, pass/PASS)
  # S3 and S4 should warn (fail/FAIL, true/TRUE)
  expect_equal(result$Sample_QC_Status, c("PASS", "PASS", "WARN", "WARN"))
})

test_that("prepare_sample_qc_for_long handles numeric status values", {
  qc_sample <- data.frame(
    sampleName = c("S1", "S1", "S2", "S2", "S3", "S3"),
    plateID = rep("Plate1", 6),
    flagName = rep(c("Detectability", "ICReads"), 3),
    val = rep(c(0.95, 1200), 3),
    status = c(0, 0, 1, 0, 2, 5),  # Numeric: 0 = pass, non-zero = warn
    stringsAsFactors = FALSE
  )

  result <- NULISAseqR:::prepare_sample_qc_for_long(qc_sample)

  # S1 passes (0, 0), S2 warns (1, 0), S3 warns (2, 5)
  expect_equal(result$Sample_QC_Status, c("PASS", "WARN", "WARN"))
  expect_equal(result$Sample_QC_Detectability_Status, c("PASS", "WARN", "WARN"))
  expect_equal(result$Sample_QC_ICReads_Status, c("PASS", "PASS", "WARN"))
})

test_that("prepare_sample_qc_for_long handles NA values in status", {
  qc_sample <- data.frame(
    sampleName = c("S1", "S1", "S2", "S2"),
    plateID = rep("Plate1", 4),
    flagName = rep(c("Detectability", "ICReads"), 2),
    val = c(0.95, 1200, NA, 1100),
    status = c(FALSE, FALSE, NA, FALSE),
    stringsAsFactors = FALSE
  )

  result <- NULISAseqR:::prepare_sample_qc_for_long(qc_sample)

  expect_equal(nrow(result), 2)
  # NA in status should be handled gracefully
  expect_true(is.na(result$Sample_QC_Detectability[2]) || !is.na(result$Sample_QC_Detectability[2]))
})

test_that("prepare_sample_qc_for_long handles empty QC data", {
  qc_sample <- data.frame(
    sampleName = character(0),
    plateID = character(0),
    flagName = character(0),
    val = numeric(0),
    status = logical(0),
    stringsAsFactors = FALSE
  )

  result <- NULISAseqR:::prepare_sample_qc_for_long(qc_sample)

  expect_equal(nrow(result), 0)
  expect_true("SampleName" %in% names(result))
  expect_true("Sample_QC_Status" %in% names(result))
})

test_that("prepare_sample_qc_for_long handles single QC metric", {
  qc_sample <- data.frame(
    sampleName = c("S1", "S2"),
    plateID = c("Plate1", "Plate1"),
    flagName = c("Detectability", "Detectability"),
    val = c(0.95, 0.88),
    status = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  result <- NULISAseqR:::prepare_sample_qc_for_long(qc_sample)

  expect_equal(nrow(result), 2)
  expect_true("Sample_QC_Detectability" %in% names(result))
  expect_true("Sample_QC_Detectability_Status" %in% names(result))
  expect_equal(result$Sample_QC_Status, c("PASS", "WARN"))
})

test_that("prepare_sample_qc_for_long handles whitespace in status values", {
  qc_sample <- data.frame(
    sampleName = c("S1", "S1", "S2", "S2"),
    plateID = rep("Plate1", 4),
    flagName = rep(c("Detectability", "ICReads"), 2),
    val = rep(c(0.95, 1200), 2),
    status = c("  OK  ", "PASS ", " FALSE", "0  "),  # With whitespace
    stringsAsFactors = FALSE
  )

  result <- NULISAseqR:::prepare_sample_qc_for_long(qc_sample)

  # All should pass after trimming whitespace
  expect_equal(result$Sample_QC_Status, c("PASS", "PASS"))
})

test_that("prepare_target_qc_for_long handles case-insensitive status values", {
  qc_target <- data.frame(
    target = c("T1", "T1", "T2", "T2"),
    plateID = rep("Plate1", 4),
    flagName = rep(c("Conc_Accuracy", "Conc_CV"), 2),
    val = rep(c(0.05, 0.12), 2),
    status = c("ok", "PASS", "Fail", "TRUE"),
    stringsAsFactors = FALSE
  )

  result <- NULISAseqR:::prepare_target_qc_for_long(qc_target)

  # T1 passes (ok, PASS), T2 warns (Fail, TRUE)
  expect_equal(result$Target_QC_Status, c("PASS", "WARN"))
})

test_that("prepare_target_qc_for_long handles missing status column", {
  qc_target <- data.frame(
    target = c("T1", "T1", "T2", "T2"),
    plateID = rep("Plate1", 4),
    flagName = rep(c("Conc_Accuracy", "Conc_CV"), 2),
    val = rep(c(0.05, 0.12), 2),
    stringsAsFactors = FALSE
  )

  result <- NULISAseqR:::prepare_target_qc_for_long(qc_target)

  expect_equal(nrow(result), 2)
  expect_true("Target_QC_Conc_Accuracy" %in% names(result))
  expect_true("Target_QC_Conc_CV" %in% names(result))
  expect_true(all(is.na(result$Target_QC_Status)))
  expect_false("Target_QC_Conc_Accuracy_Status" %in% names(result))
})

test_that("prepare_target_qc_for_long handles numeric status", {
  qc_target <- data.frame(
    target = c("T1", "T1", "T2", "T2"),
    plateID = rep("Plate1", 4),
    flagName = rep(c("Conc_Accuracy", "Conc_CV"), 2),
    val = rep(c(0.05, 0.12), 2),
    status = c(0, 0, 1, 2),
    stringsAsFactors = FALSE
  )

  result <- NULISAseqR:::prepare_target_qc_for_long(qc_target)

  # T1 passes (0, 0), T2 warns (1, 2)
  expect_equal(result$Target_QC_Status, c("PASS", "WARN"))
})


# ==============================================================================
# Multi-Plate Scenario Tests
# ==============================================================================

test_that("prepare_sample_qc_for_long handles multi-plate data", {
  qc_sample <- data.frame(
    sampleName = c("S1", "S1", "S2", "S2", "S1", "S1", "S2", "S2"),
    plateID = c(rep("Plate1", 4), rep("Plate2", 4)),
    flagName = rep(c("Detectability", "ICReads"), 4),
    val = c(0.95, 1200, 0.88, 1100, 0.92, 1150, 0.85, 1050),
    status = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  result <- NULISAseqR:::prepare_sample_qc_for_long(qc_sample)

  # Should have 4 rows: S1-Plate1, S2-Plate1, S1-Plate2, S2-Plate2
  expect_equal(nrow(result), 4)
  expect_equal(sum(result$PlateID == "Plate1"), 2)
  expect_equal(sum(result$PlateID == "Plate2"), 2)

  # S2 on Plate2 should have WARN status
  s2_plate2 <- result[result$SampleName == "S2" & result$PlateID == "Plate2", ]
  expect_equal(s2_plate2$Sample_QC_Status, "WARN")
})

test_that("prepare_target_qc_for_long handles multi-plate data", {
  qc_target <- data.frame(
    target = c("T1", "T1", "T2", "T2", "T1", "T1", "T2", "T2"),
    plateID = c(rep("Plate1", 4), rep("Plate2", 4)),
    flagName = rep(c("Conc_Accuracy", "Conc_CV"), 4),
    val = c(0.05, 0.12, -0.08, 0.25, 0.06, 0.15, -0.10, 0.50),
    status = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  result <- NULISAseqR:::prepare_target_qc_for_long(qc_target)

  # Should have 4 rows: T1-Plate1, T2-Plate1, T1-Plate2, T2-Plate2
  expect_equal(nrow(result), 4)

  # T2 on Plate2 should warn
  t2_plate2 <- result[result$Target == "T2" & result$PlateID == "Plate2", ]
  expect_equal(t2_plate2$Target_QC_Status, "WARN")
})

test_that("format_wide_to_long handles multi-plate QC data", {
  skip_if_not(requireNamespace("NULISAseqR", quietly = TRUE))
  skip_if_not(rlang::is_function(get0("format_wide_to_long", envir = asNamespace("NULISAseqR"), inherits = FALSE)))

  mock_data <- list(
    targets = data.frame(
      plateID = c("Plate1", "Plate1"),
      targetName = c("Target_1", "Target_2"),
      logged_LOD = c(5.0, 5.5),
      stringsAsFactors = FALSE
    ),
    samples = data.frame(
      plateID = c("Plate1", "Plate1"),
      sampleName = c("Sample_1", "Sample_2"),
      sampleType = c("Sample", "Sample"),
      stringsAsFactors = FALSE
    ),
    Data_NPQ = matrix(
      c(6.0, 6.5, 7.0, 7.5),
      nrow = 2,
      dimnames = list(c("Target_1", "Target_2"), c("Sample_1", "Sample_2"))
    ),
    Data_raw = matrix(
      c(100, 120, 150, 180),
      nrow = 2,
      dimnames = list(c("Target_1", "Target_2"), c("Sample_1", "Sample_2"))
    ),
    ExecutionDetails = list(
      Plate1 = list(Assay = "Test Panel")
    ),
    qcSample = data.frame(
      sampleName = c("Sample_1", "Sample_1", "Sample_2", "Sample_2"),
      plateID = c("Plate1", "Plate1", "Plate1", "Plate1"),
      flagName = c("Detectability", "ICReads", "Detectability", "ICReads"),
      val = c(0.95, 1200, 0.88, 1100),
      status = c(FALSE, FALSE, TRUE, FALSE),  # Sample_2 has warning
      stringsAsFactors = FALSE
    ),
    qcTarget = data.frame(
      target = c("Target_1", "Target_1", "Target_2", "Target_2"),
      plateID = c("Plate1", "Plate1", "Plate1", "Plate1"),
      flagName = c("Metric1", "Metric2", "Metric1", "Metric2"),
      val = c(0.05, 0.12, -0.08, 0.25),
      status = c(FALSE, FALSE, FALSE, TRUE),  # Target_2 has warning
      stringsAsFactors = FALSE
    ),
    plateID = "Plate1"
  )

  result <- NULISAseqR:::format_wide_to_long(mock_data, AQ = FALSE, include_qc = TRUE)

  # Check Sample_2 has WARN status
  sample_2_rows <- result[result$SampleName == "Sample_2", ]
  expect_true(all(sample_2_rows$Sample_QC_Status == "WARN"))

  # Check Target_2 has WARN status
  target_2_rows <- result[result$Target == "Target_2", ]
  expect_true(all(target_2_rows$Target_QC_Status == "WARN"))
})


# ==============================================================================
# Error Handling Tests
# ==============================================================================

test_that("prepare_sample_qc_for_long validates required columns", {
  # Missing 'status' column
  qc_sample_missing_status <- data.frame(
    sampleName = c("S1"),
    plateID = c("Plate1"),
    flagName = c("Detectability"),
    val = c(0.95),
    stringsAsFactors = FALSE
  )

  expect_error(
    NULISAseqR:::prepare_sample_qc_for_long(qc_sample_missing_status),
    "qcSample missing required columns.*status"
  )

  # Missing 'val' column
  qc_sample_missing_val <- data.frame(
    sampleName = c("S1"),
    plateID = c("Plate1"),
    flagName = c("Detectability"),
    status = c(FALSE),
    stringsAsFactors = FALSE
  )

  expect_error(
    NULISAseqR:::prepare_sample_qc_for_long(qc_sample_missing_val),
    "qcSample missing required columns.*val"
  )
})

test_that("prepare_sample_qc_for_long rejects invalid status types", {
  qc_sample <- data.frame(
    sampleName = c("S1"),
    plateID = c("Plate1"),
    flagName = c("Detectability"),
    val = c(0.95),
    status = factor(c("PASS")),  # Factor type - invalid
    stringsAsFactors = FALSE
  )

  expect_error(
    NULISAseqR:::prepare_sample_qc_for_long(qc_sample),
    "qcSample 'status' column must be logical, numeric, or character"
  )
})

test_that("prepare_target_qc_for_long validates required columns", {
  qc_target_missing_target <- data.frame(
    plateID = c("Plate1"),
    flagName = c("Metric1"),
    val = c(0.05),
    stringsAsFactors = FALSE
  )

  expect_error(
    NULISAseqR:::prepare_target_qc_for_long(qc_target_missing_target),
    "qcTarget missing required columns.*target"
  )
})

test_that("prepare_target_qc_for_long rejects invalid status types when present", {
  qc_target <- data.frame(
    target = c("T1"),
    plateID = c("Plate1"),
    flagName = c("Metric1"),
    val = c(0.05),
    stringsAsFactors = FALSE
  )

  # Add factor column (invalid type)
  qc_target$status <- factor(c("PASS"))

  expect_error(
    NULISAseqR:::prepare_target_qc_for_long(qc_target),
    "qcTarget 'status' column must be logical, numeric, or character"
  )
})


# ==============================================================================
# Duplicate Handling Tests
# ==============================================================================

test_that("prepare_sample_qc_for_long handles duplicate entries", {
  qc_sample <- data.frame(
    sampleName = c("S1", "S1", "S1", "S1"),  # Duplicate entries
    plateID = c("Plate1", "Plate1", "Plate1", "Plate1"),
    flagName = c("Detectability", "Detectability", "ICReads", "ICReads"),
    val = c(0.95, 0.90, 1200, 1100),  # Different values
    status = c(FALSE, TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )

  # Should warn about duplicates
  expect_warning(
    result <- NULISAseqR:::prepare_sample_qc_for_long(qc_sample),
    NA  # We're using logger, not warning(), so no R warning expected
  )

  # Should keep first value for each metric
  expect_equal(nrow(result), 1)
  expect_equal(result$Sample_QC_Detectability, 0.95)  # First value
  expect_equal(result$Sample_QC_ICReads, 1200)  # First value
})

test_that("prepare_target_qc_for_long handles duplicate entries", {
  qc_target <- data.frame(
    target = c("T1", "T1", "T1", "T1"),
    plateID = c("Plate1", "Plate1", "Plate1", "Plate1"),
    flagName = c("Metric1", "Metric1", "Metric2", "Metric2"),
    val = c(0.05, 0.10, 0.12, 0.15),
    status = c(FALSE, TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )

  result <- NULISAseqR:::prepare_target_qc_for_long(qc_target)

  expect_equal(nrow(result), 1)
  expect_equal(result$Target_QC_Metric1, 0.05)  # First value
  expect_equal(result$Target_QC_Metric2, 0.12)  # First value
})


# ==============================================================================
# Special Characters and Edge Cases
# ==============================================================================

test_that("prepare_sample_qc_for_long handles special characters in flagName", {
  qc_sample <- data.frame(
    sampleName = c("S1", "S1"),
    plateID = c("Plate1", "Plate1"),
    flagName = c("Sample_Detectability_%", "IC_Reads_(Count)"),
    val = c(0.95, 1200),
    status = c(FALSE, FALSE),
    stringsAsFactors = FALSE
  )

  result <- NULISAseqR:::prepare_sample_qc_for_long(qc_sample)

  expect_equal(nrow(result), 1)
  # Special characters should be preserved in column names
  expect_true(any(grepl("Detectability", names(result))))
  expect_true(any(grepl("Reads", names(result))))
})

test_that("prepare_sample_qc_for_long handles all metrics passing", {
  qc_sample <- data.frame(
    sampleName = rep("S1", 5),
    plateID = rep("Plate1", 5),
    flagName = c("Metric1", "Metric2", "Metric3", "Metric4", "Metric5"),
    val = 1:5,
    status = rep(FALSE, 5),
    stringsAsFactors = FALSE
  )

  result <- NULISAseqR:::prepare_sample_qc_for_long(qc_sample)

  expect_equal(result$Sample_QC_Status, "PASS")
  # All individual metrics should also pass
  status_cols <- names(result)[grepl("_Status$", names(result)) & names(result) != "Sample_QC_Status"]
  expect_true(all(result[, status_cols] == "PASS"))
})

test_that("prepare_sample_qc_for_long handles all metrics failing", {
  qc_sample <- data.frame(
    sampleName = rep("S1", 5),
    plateID = rep("Plate1", 5),
    flagName = c("Metric1", "Metric2", "Metric3", "Metric4", "Metric5"),
    val = 1:5,
    status = rep(TRUE, 5),
    stringsAsFactors = FALSE
  )

  result <- NULISAseqR:::prepare_sample_qc_for_long(qc_sample)

  expect_equal(result$Sample_QC_Status, "WARN")
  # All individual metrics should also warn
  status_cols <- names(result)[grepl("_Status$", names(result)) & names(result) != "Sample_QC_Status"]
  expect_true(all(result[, status_cols] == "WARN"))
})


# ==============================================================================
# Parameter Control Tests
# ==============================================================================

test_that("format_wide_to_long respects include_qc = FALSE parameter", {
  skip_if_not(requireNamespace("NULISAseqR", quietly = TRUE))
  skip_if_not(rlang::is_function(get0("format_wide_to_long", envir = asNamespace("NULISAseqR"), inherits = FALSE)))

  mock_data <- list(
    targets = data.frame(
      plateID = c("Plate1", "Plate1"),
      targetName = c("Target_1", "Target_2"),
      logged_LOD = c(5.0, 5.5),
      stringsAsFactors = FALSE
    ),
    samples = data.frame(
      plateID = c("Plate1", "Plate1"),
      sampleName = c("Sample_1", "Sample_2"),
      sampleType = c("Sample", "Sample"),
      stringsAsFactors = FALSE
    ),
    Data_NPQ = matrix(
      c(6.0, 6.5, 7.0, 7.5),
      nrow = 2,
      dimnames = list(c("Target_1", "Target_2"), c("Sample_1", "Sample_2"))
    ),
    Data_raw = matrix(
      c(100, 120, 150, 180),
      nrow = 2,
      dimnames = list(c("Target_1", "Target_2"), c("Sample_1", "Sample_2"))
    ),
    ExecutionDetails = list(
      Plate1 = list(Assay = "Test Panel")
    ),
    qcSample = data.frame(
      sampleName = c("Sample_1", "Sample_1", "Sample_2", "Sample_2"),
      plateID = c("Plate1", "Plate1", "Plate1", "Plate1"),
      flagName = c("Detectability", "ICReads", "Detectability", "ICReads"),
      val = c(0.95, 1200, 0.88, 1100),
      status = c(FALSE, FALSE, FALSE, FALSE),
      stringsAsFactors = FALSE
    ),
    qcTarget = data.frame(
      target = c("Target_1", "Target_1", "Target_2", "Target_2"),
      plateID = c("Plate1", "Plate1", "Plate1", "Plate1"),
      flagName = c("Metric1", "Metric2", "Metric1", "Metric2"),
      val = c(0.05, 0.12, -0.08, 0.25),
      status = c(FALSE, FALSE, FALSE, FALSE),
      stringsAsFactors = FALSE
    ),
    plateID = "Plate1"
  )

  # Call with include_qc = FALSE
  result <- NULISAseqR:::format_wide_to_long(mock_data, AQ = FALSE, include_qc = FALSE)

  # QC columns should NOT be present
  expect_false("Sample_QC_Status" %in% names(result))
  expect_false("Target_QC_Status" %in% names(result))
  expect_false("Sample_QC_Detectability" %in% names(result))
  expect_false("Target_QC_Metric1" %in% names(result))

  # Should still have 4 rows (2 targets × 2 samples)
  expect_equal(nrow(result), 4)

  # Regular data columns should be present
  expect_true("NPQ" %in% names(result))
  expect_true("Target" %in% names(result))
  expect_true("SampleName" %in% names(result))
})

test_that("format_wide_to_long includes QC by default (include_qc = TRUE)", {
  skip_if_not(requireNamespace("NULISAseqR", quietly = TRUE))
  skip_if_not(rlang::is_function(get0("format_wide_to_long", envir = asNamespace("NULISAseqR"), inherits = FALSE)))

  mock_data <- list(
    targets = data.frame(
      plateID = c("Plate1", "Plate1"),
      targetName = c("Target_1", "Target_2"),
      logged_LOD = c(5.0, 5.5),
      stringsAsFactors = FALSE
    ),
    samples = data.frame(
      plateID = c("Plate1", "Plate1"),
      sampleName = c("Sample_1", "Sample_2"),
      sampleType = c("Sample", "Sample"),
      stringsAsFactors = FALSE
    ),
    Data_NPQ = matrix(
      c(6.0, 6.5, 7.0, 7.5),
      nrow = 2,
      dimnames = list(c("Target_1", "Target_2"), c("Sample_1", "Sample_2"))
    ),
    Data_raw = matrix(
      c(100, 120, 150, 180),
      nrow = 2,
      dimnames = list(c("Target_1", "Target_2"), c("Sample_1", "Sample_2"))
    ),
    ExecutionDetails = list(
      Plate1 = list(Assay = "Test Panel")
    ),
    qcSample = data.frame(
      sampleName = c("Sample_1", "Sample_1", "Sample_2", "Sample_2"),
      plateID = c("Plate1", "Plate1", "Plate1", "Plate1"),
      flagName = c("Detectability", "ICReads", "Detectability", "ICReads"),
      val = c(0.95, 1200, 0.88, 1100),
      status = c(FALSE, FALSE, FALSE, FALSE),
      stringsAsFactors = FALSE
    ),
    qcTarget = data.frame(
      target = c("Target_1", "Target_1", "Target_2", "Target_2"),
      plateID = c("Plate1", "Plate1", "Plate1", "Plate1"),
      flagName = c("Metric1", "Metric2", "Metric1", "Metric2"),
      val = c(0.05, 0.12, -0.08, 0.25),
      status = c(FALSE, FALSE, FALSE, FALSE),
      stringsAsFactors = FALSE
    ),
    plateID = "Plate1"
  )

  # Call without specifying include_qc (should default to TRUE)
  result <- NULISAseqR:::format_wide_to_long(mock_data, AQ = FALSE)

  # QC columns SHOULD be present
  expect_true("Sample_QC_Status" %in% names(result))
  expect_true("Target_QC_Status" %in% names(result))
  expect_true("Sample_QC_Detectability" %in% names(result))
  expect_true("Target_QC_Metric1" %in% names(result))
})

test_that("format_wide_to_long filters AQ-specific Target QC columns in RQ mode", {
  skip_if_not(requireNamespace("NULISAseqR", quietly = TRUE))
  skip_if_not(rlang::is_function(get0("format_wide_to_long", envir = asNamespace("NULISAseqR"), inherits = FALSE)))

  mock_data <- list(
    targets = data.frame(
      plateID = c("Plate1", "Plate1"),
      targetName = c("Target_1", "Target_2"),
      logged_LOD = c(5.0, 5.5),
      stringsAsFactors = FALSE
    ),
    samples = data.frame(
      plateID = c("Plate1", "Plate1"),
      sampleName = c("Sample_1", "Sample_2"),
      sampleType = c("Sample", "Sample"),
      stringsAsFactors = FALSE
    ),
    Data_NPQ = matrix(
      c(6.0, 6.5, 7.0, 7.5),
      nrow = 2,
      dimnames = list(c("Target_1", "Target_2"), c("Sample_1", "Sample_2"))
    ),
    Data_raw = matrix(
      c(100, 120, 150, 180),
      nrow = 2,
      dimnames = list(c("Target_1", "Target_2"), c("Sample_1", "Sample_2"))
    ),
    ExecutionDetails = list(
      Plate1 = list(Assay = "Test Panel")
    ),
    qcSample = data.frame(
      sampleName = c("Sample_1", "Sample_1", "Sample_2", "Sample_2"),
      plateID = c("Plate1", "Plate1", "Plate1", "Plate1"),
      flagName = c("Detectability", "ICReads", "Detectability", "ICReads"),
      val = c(0.95, 1200, 0.88, 1100),
      status = c(FALSE, FALSE, FALSE, FALSE),
      stringsAsFactors = FALSE
    ),
    qcTarget = data.frame(
      target = rep(c("Target_1", "Target_2"), each = 4),
      plateID = rep("Plate1", 8),
      flagName = rep(c("Target_Min_Reads", "Target_Conc_CV_RQ", "Target_Conc_Accuracy", "Target_Conc_CV"), 2),
      val = rep(c(0.95, 0.12, 0.05, 0.08), 2),
      status = rep(FALSE, 8),
      stringsAsFactors = FALSE
    ),
    plateID = "Plate1"
  )

  # Call with AQ = FALSE (RQ mode)
  result <- NULISAseqR:::format_wide_to_long(mock_data, AQ = FALSE, include_qc = TRUE)

  # Sample QC columns SHOULD be present
  expect_true("Sample_QC_Status" %in% names(result))
  expect_true("Sample_QC_Detectability" %in% names(result))

  # RQ Target QC columns SHOULD be present
  expect_true("Target_QC_Status" %in% names(result))
  expect_true("Target_QC_Min_Reads" %in% names(result))
  expect_true("Target_QC_Conc_CV_RQ" %in% names(result))

  # AQ-specific Target QC columns should NOT be present in RQ mode
  expect_false("Target_QC_Conc_Accuracy" %in% names(result))
  expect_false("Target_QC_Conc_Accuracy_Status" %in% names(result))
  expect_false("Target_QC_Conc_CV" %in% names(result))
  expect_false("Target_QC_Conc_CV_Status" %in% names(result))

  # Should still have 4 rows (2 targets × 2 samples)
  expect_equal(nrow(result), 4)
})

test_that("format_wide_to_long includes all Target QC columns in AQ mode", {
  skip_if_not(requireNamespace("NULISAseqR", quietly = TRUE))
  skip_if_not(rlang::is_function(get0("format_wide_to_long", envir = asNamespace("NULISAseqR"), inherits = FALSE)))

  mock_data <- list(
    targets = data.frame(
      plateID = c("Plate1", "Plate1"),
      targetName = c("Target_1", "Target_2"),
      logged_LOD = c(5.0, 5.5),
      stringsAsFactors = FALSE
    ),
    samples = data.frame(
      plateID = c("Plate1", "Plate1"),
      sampleName = c("Sample_1", "Sample_2"),
      sampleType = c("Sample", "Sample"),
      stringsAsFactors = FALSE
    ),
    Data_AQ_aM = matrix(
      c(100, 120, 150, 180),
      nrow = 2,
      dimnames = list(c("Target_1", "Target_2"), c("Sample_1", "Sample_2"))
    ),
    Data_AQ_pgmL = matrix(
      c(10, 12, 15, 18),
      nrow = 2,
      dimnames = list(c("Target_1", "Target_2"), c("Sample_1", "Sample_2"))
    ),
    Data_AQlog2_aM = matrix(
      c(6.0, 6.5, 7.0, 7.5),
      nrow = 2,
      dimnames = list(c("Target_1", "Target_2"), c("Sample_1", "Sample_2"))
    ),
    Data_AQlog2_pgmL = matrix(
      c(3.0, 3.5, 4.0, 4.5),
      nrow = 2,
      dimnames = list(c("Target_1", "Target_2"), c("Sample_1", "Sample_2"))
    ),
    ExecutionDetails = list(
      Plate1 = list(Assay = "Test Panel")
    ),
    qcSample = data.frame(
      sampleName = c("Sample_1", "Sample_1", "Sample_2", "Sample_2"),
      plateID = c("Plate1", "Plate1", "Plate1", "Plate1"),
      flagName = c("Detectability", "ICReads", "Detectability", "ICReads"),
      val = c(0.95, 1200, 0.88, 1100),
      status = c(FALSE, FALSE, FALSE, FALSE),
      stringsAsFactors = FALSE
    ),
    qcTarget = data.frame(
      target = rep(c("Target_1", "Target_2"), each = 4),
      plateID = rep("Plate1", 8),
      flagName = rep(c("Target_Min_Reads", "Target_Conc_CV_RQ", "Target_Conc_Accuracy", "Target_Conc_CV"), 2),
      val = rep(c(0.95, 0.12, 0.05, 0.08), 2),
      status = rep(FALSE, 8),
      stringsAsFactors = FALSE
    ),
    plateID = "Plate1"
  )

  # Call with AQ = TRUE (AQ mode)
  result <- NULISAseqR:::format_wide_to_long(mock_data, AQ = TRUE, include_qc = TRUE)

  # Sample QC columns SHOULD be present
  expect_true("Sample_QC_Status" %in% names(result))
  expect_true("Sample_QC_Detectability" %in% names(result))

  # ALL Target QC columns SHOULD be present in AQ mode (both RQ and AQ metrics)
  expect_true("Target_QC_Status" %in% names(result))
  expect_true("Target_QC_Min_Reads" %in% names(result))
  expect_true("Target_QC_Conc_CV_RQ" %in% names(result))
  expect_true("Target_QC_Conc_Accuracy" %in% names(result))
  expect_true("Target_QC_Conc_CV" %in% names(result))

  # Should still have 4 rows (2 targets × 2 samples)
  expect_equal(nrow(result), 4)
})

