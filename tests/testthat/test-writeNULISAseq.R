test_that("writeNULISAseq creates Excel file with RQ data (no AQ)", {
  # Get fixture file path
  fixture_file <- test_path("fixtures", "XML_v1.3.0.xml")

  # Create temp directory and copy fixture there
  temp_dir <- tempdir()
  input_file <- file.path(temp_dir, basename(fixture_file))
  file.copy(fixture_file, input_file, overwrite = TRUE)

  # Create temp output file in the same temp directory
  output_file <- tempfile(tmpdir = temp_dir, fileext = ".xlsx")

  # Write the Excel file (RQ only, no AQ)
  writeNULISAseq(xml_files = basename(input_file),
                 dataDir = temp_dir,
                 output_filename = basename(output_file),
                 output_TAP_AQ = FALSE,
                 Panel = '200-Plex Inflammation Panel v2',
                 PanelLotNumber = 'panelLot022',
                 metadata = list(NULISAseqR='dlfjqp13948eo',
                                 NULISAseqAQ='e98249lkfjqke'),
                 include_unnorm_counts = TRUE)

  # Output file should be in temp_dir (no need to move)
  
  # Check file was created
  expect_true(file.exists(output_file))
  
  # Read the Excel file to check structure
  sheets <- readxl::excel_sheets(output_file)
  
  # Expected sheets for RQ-only
  expected_sheets <- c("NPQ Values", "Target Detectability", 
                       "Sample Information", "Analysis Metadata")
  
  expect_true(all(expected_sheets %in% sheets),
              info = "Expected RQ sheets should be present")
  
  # AQ sheets should NOT be present
  expect_false("AQ (pg per mL and aM)" %in% sheets,
               info = "AQ sheet should not be present when output_TAP_AQ=FALSE")
  expect_false("Target Quantifiability" %in% sheets,
               info = "Target Quantifiability sheet should not be present")
  
  # Check column structure of NPQ Values sheet
  npq_data <- readxl::read_excel(output_file, sheet = "NPQ Values")
  expected_cols <- c("Panel", "PanelLotNumber", "PlateID", "SampleName", 
                     "SampleType", "Target", "AlamarTargetID", "UniProtID", 
                     "ProteinName", "SampleQC", "LOD", "UnnormalizedCount", "NPQ")
  
  expect_true(all(expected_cols %in% colnames(npq_data)),
              info = "NPQ Values should have expected columns")
  
  # Check Target Detectability sheet structure
  detect_data <- readxl::read_excel(output_file, sheet = "Target Detectability")
  expect_true("Target" %in% colnames(detect_data),
              info = "Target Detectability should have Target column")
  expect_true("plasma (n = 86)" %in% tolower(colnames(detect_data)),
              info = "Target Detectability should have plasma (n = 86) column")
  
  # Check Sample Information sheet structure
  sample_info <- readxl::read_excel(output_file, sheet = "Sample Information")
  expect_true(all(c("PlateID", "WellPosition", "SampleName") %in% colnames(sample_info)),
              info = "Sample Information should have expected columns")
  
  # Check specific data values for a known target and sample
  test_row <- npq_data[npq_data$Target == "AGER" &
                         npq_data$SampleName == "Sample01", ]

  # Assert that the expected row exists
  expect_equal(nrow(test_row), 1,
               label = "Should find exactly one row for AGER/Sample01")

  # Check character columns
  expect_equal(test_row$Panel, "200-Plex Inflammation Panel v2")
  expect_equal(test_row$PanelLotNumber, "panelLot022")
  expect_equal(test_row$PlateID, "Plate_01")
  expect_equal(test_row$SampleName, "Sample01")
  expect_equal(test_row$SampleType, "Sample")
  expect_equal(test_row$Target, "AGER")
  expect_equal(test_row$AlamarTargetID, "t5521")
  expect_equal(test_row$UniProtID, "Q15109")
  expect_equal(test_row$ProteinName, "Advanced glycosylation end product-specific receptor")
  expect_equal(test_row$SampleQC, "PASS")

  # Check numeric columns
  expect_equal(test_row$LOD, 17.1221185470271, tolerance = 1e-4,
               label = "LOD for AGER/Sample01")
  expect_equal(test_row$UnnormalizedCount, 119,
               label = "UnnormalizedCount for AGER/Sample01")
  expect_equal(test_row$NPQ, 15.1835336023668, tolerance = 1e-4,
               label = "NPQ for AGER/Sample01")

  # Clean up
  unlink(c(input_file, output_file))
})

test_that("writeNULISAseq creates Excel file with AQ data", {
  skip_if_not(requireNamespace("NULISAseqAQ", quietly = TRUE),
              "NULISAseqAQ package not available")

  # Get fixture file path
  fixture_file <- test_path("fixtures", "XML_v1.3.0_with_AQ.xml")

  # Create temp directory and copy fixture there
  temp_dir <- tempdir()
  input_file <- file.path(temp_dir, basename(fixture_file))
  file.copy(fixture_file, input_file, overwrite = TRUE)

  # Create temp output file in the same temp directory
  output_file <- tempfile(tmpdir = temp_dir, fileext = ".xlsx")

  # Write the Excel file with AQ data
  writeNULISAseq(xml_files = basename(input_file),
                 dataDir = temp_dir,
                 output_filename = basename(output_file),
                 output_TAP_AQ = TRUE,
                 Panel = '200-Plex Inflammation Panel AQ v2',
                 PanelLotNumber = 'panelLot022',
                 metadata = list(NULISAseqR='dlfjqp13948eo',
                                 NULISAseqAQ='e98249lkfjqke'),
                 include_unnorm_counts = TRUE)

  # Output file should be in temp_dir (no need to move)
  
  # Check file was created
  expect_true(file.exists(output_file))
  
  # Read the Excel file to check structure
  sheets <- readxl::excel_sheets(output_file)
  
  # Expected sheets for RQ + AQ
  expected_sheets <- c("RQ NPQ Values", "AQ (pg per mL and aM)", 
                       "Target Detectability", "Target Quantifiability",
                       "Sample Information", "Analysis Metadata")
  
  expect_true(all(expected_sheets %in% sheets),
              info = "Expected RQ + AQ sheets should be present")
  
  # Check column structure of RQ NPQ Values sheet
  rq_data <- readxl::read_excel(output_file, sheet = "RQ NPQ Values")
  expected_rq_cols <- c("Panel", "PanelLotNumber", "PlateID", "SampleName", 
                        "SampleType", "Target", "AlamarTargetID", "UniProtID", 
                        "ProteinName", "SampleQC", "LOD", "UnnormalizedCount", "NPQ")
  
  expect_true(all(expected_rq_cols %in% colnames(rq_data)),
              info = "RQ NPQ Values should have expected columns")
  
  # Check column structure of AQ sheet (starts at row 3 due to header)
  aq_data <- readxl::read_excel(output_file, sheet = "AQ (pg per mL and aM)", skip = 2)
  expected_aq_cols <- c("Panel", "PanelLotNumber", "PlateID", "SampleName", 
                        "SampleType", "Target", "AlamarTargetID", "UniProtID", 
                        "ProteinName", "SampleQC", "TargetQC", 
                        "Conc (pg/mL)", "LOD (pg/mL)", "LLOQ (pg/mL)", "ULOQ (pg/mL)",
                        "Conc (aM)", "LOD (aM)", "LLOQ (aM)", "ULOQ (aM)")
  
  expect_true(all(expected_aq_cols %in% colnames(aq_data)),
              info = "AQ sheet should have expected columns")
  
  # Check Target Detectability sheet structure
  detect_data <- readxl::read_excel(output_file, sheet = "Target Detectability")
  expect_true("Target" %in% colnames(detect_data),
              info = "Target Detectability should have Target column")
  expect_true(any(tolower(colnames(detect_data)) == "other (n = 86)"),
              info = "Target Detectability should have other (n = 86) column")
  
  # Check detectability has expected number of rows (one per target)
  expect_equal(nrow(detect_data), 249,
               label = "Should have detectability for 249 targets")
  
  # Spot-check specific target detectability values (first 3 targets)
  test_targets <- c("AGER", "AGRP", "ANGPT1")
  expected_detect <- list(
    "AGER" = 58.1,
    "AGRP" = 44.2,
    "ANGPT1" = 25.6
  )
  
  # Find the sample type column (e.g., "other (n = 86)")
  sample_col <- grep("\\(n = ", colnames(detect_data), value = TRUE)
  
  if (length(sample_col) > 0) {
    for (target in test_targets) {
      detect_row <- detect_data[detect_data$Target == target, ]
      # Assert that the expected target exists
      expect_equal(nrow(detect_row), 1,
                   label = sprintf("Should find exactly one row for target %s", target))
      # Convert from character to numeric (readxl reads as character)
      actual_value <- as.numeric(detect_row[[sample_col[1]]])
      expect_equal(actual_value, expected_detect[[target]], tolerance = 0.001,
                   label = sprintf("%s detectability in %s", target, sample_col[1]))
    }
  }
  
  # Check that reverse curve targets are present but have NA for LOD and detectability
  reverse_targets <- c("CRP", "KNG1")  # Replace with actual reverse targets from your data

  for (target in reverse_targets) {
    # Check target exists in RQ data
    rq_row <- rq_data[rq_data$Target == target, ]
    expect_gt(nrow(rq_row), 0,
              label = sprintf("Reverse curve target %s should be in RQ data", target))

    # LOD should be NA for reverse curve targets
    expect_true(all(is.na(rq_row$LOD)),
                info = sprintf("Reverse curve target %s should have NA for LOD", target))

    # Check target exists in detectability data
    detect_row <- detect_data[detect_data$Target == target, ]
    expect_gt(nrow(detect_row), 0,
              label = sprintf("Reverse curve target %s should be in detectability table", target))

    # Detectability should be NA for reverse curve targets
    sample_col <- grep("\\(n = ", colnames(detect_data), value = TRUE)
    if (length(sample_col) > 0) {
      detect_value <- detect_row[[sample_col[1]]]
      # Check if it's character "NA" or actual NA
      expect_true(is.na(detect_value) || detect_value == "NA",
                  info = sprintf("Reverse curve target %s should have NA detectability", target))
    }
  }
  
  # Check Target Quantifiability sheet structure
  quant_data <- readxl::read_excel(output_file, sheet = "Target Quantifiability")
  expect_true("Target" %in% colnames(quant_data),
              info = "Target Quantifiability should have Target column")
  expect_true(any(tolower(colnames(quant_data)) == "other (n = 86)"),
              info = "Target Quantifiability should have other (n = 86) column")
  
  # Check quantifiability has expected number of rows
  expect_equal(nrow(quant_data), 157,
               label = "Should have quantifiability for 157 targets")
  
  # Spot-check specific target quantifiability values (first 4 targets)
  test_quant_targets <- c("AGER", "AGRP", "ANGPT1", "BST2")
  expected_quant <- list(
    "AGER" = 0,
    "AGRP" = 0,
    "ANGPT1" = 25.6,
    "BST2" = 90.7
  )
  
  # Find the sample type column (e.g., "other (n = 86)")
  sample_col <- grep("\\(n = ", colnames(quant_data), value = TRUE)
  
  if (length(sample_col) > 0) {
    for (target in test_quant_targets) {
      quant_row <- quant_data[quant_data$Target == target, ]
      # Assert that the expected target exists
      expect_equal(nrow(quant_row), 1,
                   label = sprintf("Should find exactly one row for target %s", target))
      actual_value <- quant_row[[sample_col[1]]]
      expect_equal(actual_value, expected_quant[[target]], tolerance = 0.5,
                   label = sprintf("%s quantifiability in %s", target, sample_col[1]))
    }
  }
  
  # Check Sample Information sheet structure
  sample_info <- readxl::read_excel(output_file, sheet = "Sample Information")
  expect_true(all(c("PlateID", "WellPosition", "SampleName") %in% colnames(sample_info)),
              info = "Sample Information should have expected columns")
  
  # Check specific RQ data values for a known target and sample
  test_row_rq <- rq_data[rq_data$Target == "AGRP" &
                           rq_data$SampleName == "25307004091", ]

  # Assert that the expected row exists
  expect_equal(nrow(test_row_rq), 1,
               label = "Should find exactly one row for AGRP/25307004091 in RQ data")

  # Check character columns
  expect_equal(test_row_rq$Panel, "200-Plex Inflammation Panel AQ v2")
  expect_equal(test_row_rq$PanelLotNumber, "panelLot022")
  expect_equal(test_row_rq$PlateID, "Plate_01")
  expect_equal(test_row_rq$SampleName, "25307004091")
  expect_equal(test_row_rq$SampleType, "Sample")
  expect_equal(test_row_rq$Target, "AGRP")
  expect_equal(test_row_rq$AlamarTargetID, "t5469")
  expect_equal(test_row_rq$UniProtID, "O00253")
  expect_equal(test_row_rq$ProteinName, "Agouti-related protein")
  expect_equal(test_row_rq$SampleQC, "PASS")

  # Check numeric columns
  expect_equal(test_row_rq$LOD, 6.53811867553211, tolerance = 1e-4,
               label = "LOD for AGRP/Sample01")
  expect_equal(test_row_rq$UnnormalizedCount, 288,
               label = "UnnormalizedCount for AGRP/Sample01")
  expect_equal(test_row_rq$NPQ, 7.11945270878005, tolerance = 1e-4,
               label = "NPQ for AGRP/Sample01")
  
  # Check specific AQ data values for a known target and sample
  test_row_aq <- aq_data[aq_data$Target == "AGRP" &
                           aq_data$SampleName == "25307004091", ]

  # Assert that the expected row exists
  expect_equal(nrow(test_row_aq), 1,
               label = "Should find exactly one row for AGRP/25307004091 in AQ data")

  # Check character columns
  expect_equal(test_row_aq$Panel, "200-Plex Inflammation Panel AQ v2")
  expect_equal(test_row_aq$PanelLotNumber, "panelLot022")
  expect_equal(test_row_aq$PlateID, "Plate_01")
  expect_equal(test_row_aq$SampleName, "25307004091")
  expect_equal(test_row_aq$SampleType, "Sample")
  expect_equal(test_row_aq$Target, "AGRP")
  expect_equal(test_row_aq$AlamarTargetID, "t5469")
  expect_equal(test_row_aq$UniProtID, "O00253")
  expect_equal(test_row_aq$ProteinName, "Agouti-related protein")
  expect_equal(test_row_aq$SampleQC, "PASS")
  expect_equal(test_row_aq$TargetQC, "PASS")

  # Check numeric columns (pg/mL)
  expect_equal(test_row_aq$`Conc (pg/mL)`, 0.647043583536688, tolerance = 1e-4,
               label = "Conc (pg/mL) for AGRP/25307004091")
  expect_equal(test_row_aq$`LOD (pg/mL)`, 0.328186219190365, tolerance = 1e-4,
               label = "LOD (pg/mL) for AGRP/25307004091")
  expect_equal(test_row_aq$`LLOQ (pg/mL)`, 4.61168846686682, tolerance = 1e-4,
               label = "LLOQ (pg/mL) for AGRP/25307004091")
  expect_equal(test_row_aq$`ULOQ (pg/mL)`, 258243.928426144, tolerance = 1e-4,
               label = "ULOQ (pg/mL) for AGRP/25307004091")

  # Check numeric columns (aM)
  expect_equal(test_row_aq$`Conc (aM)`, 44809.1124332887, tolerance = 1e-4,
               label = "Conc (aM) for AGRP/25307004091")
  expect_equal(test_row_aq$`LOD (aM)`, 22727.577506258, tolerance = 1e-4,
               label = "LOD (aM) for AGRP/25307004091")
  expect_equal(test_row_aq$`LLOQ (aM)`, 319369.007400749, tolerance = 1e-4,
               label = "LLOQ (aM) for AGRP/25307004091")
  expect_equal(test_row_aq$`ULOQ (aM)`, 17883928561.3673, tolerance = 1e-4,
               label = "ULOQ (aM) for AGRP/25307004091")

  # Clean up
  unlink(c(input_file, output_file))
})
