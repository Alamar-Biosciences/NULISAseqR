# Define a test for QCSampleCriteria
test_that("skeleton.Rmd can be rendered into HTML", {
  test_path <- paste0(testthat::test_path(), "./../inst/rmarkdown/templates/nulisaseq/skeleton/")
  input <- paste0(test_path, "skeleton.Rmd")
  files <- c("Analysis_INF250_Lot4_AQ_LC_R3_20241229.xml", "Analysis_LC_R3_INF250_Lot3_20240727.xml")
  output <- paste0(test_path, "skeleton.html")
  fixtures_dir <- normalizePath(testthat::test_path("fixtures"))
  withr::with_tempfile(output, {
    rmarkdown::render(input, params=list(dataDir=fixtures_dir, xmlFiles=files))  
    expect_true(file.exists(output), output) 
  })
})
test_that("skeleton.Rmd can be rendered into HTML for XMLs with differing targets", {
  test_path <- paste0(testthat::test_path(), "./../inst/rmarkdown/templates/nulisaseq/skeleton/")
  input <- paste0(test_path, "skeleton.Rmd")
  output <- paste0(test_path, "skeleton.html")
  files <- c("Analysis_INF250_Lot4_AQ_LC_R3_20241229.xml", "Analysis_LC_R3_INF250_Lot3_20240727.xml")
  fixtures_dir <- normalizePath(testthat::test_path("fixtures"))
  withr::with_tempfile(output, {
    rmarkdown::render(input, params=list(dataDir=fixtures_dir, xmlFiles=files))  
    expect_true(file.exists(output), output)
  })
})
test_that("skeleton.Rmd can be rendered into HTML with an alternative IC", {
  test_path <- paste0(testthat::test_path(), "./../inst/rmarkdown/templates/nulisaseq/skeleton/")
  input <- paste0(test_path, "skeleton.Rmd")
  output <- paste0(test_path, "skeleton.html")
  files <- c("detectability_P1_Tr03_typemCherry_CCL7_hide2.xml")
  fixtures_dir <- normalizePath(testthat::test_path("fixtures"))
  withr::with_tempfile(output, {
    rmarkdown::render(input, params=list(dataDir=fixtures_dir, xmlFiles=files))  
    expect_true(file.exists(output), output)
    orig_content <- readLines(output, warn=FALSE)
    patterns <- c("-mCherry-summary", files)
    html_content <- gsub(paste(patterns, collapse = "|"), "", orig_content)
    expect_true(any(grepl("CCL7", html_content)), info="Output does not contain 'CCL7'")
    expect_false(any(grepl("mCherry", html_content)), info="Output does contain 'mCherry'")
    expect_false(any(grepl("WNT7A", html_content)), info="Output does contain 'WNT7A'")
  })
})
test_that("skeleton.Rmd can be rendered into HTML, NAS version", { 
  test_path <- paste0(testthat::test_path(), "./../inst/rmarkdown/templates/nulisaseq/skeleton/")
  input <- paste0(test_path, "skeleton.Rmd")
  files <- c("Analysis_INF250_Lot4_AQ_LC_R3_20241229.xml", "Analysis_LC_R3_INF250_Lot3_20240727.xml")
  output <- paste0(test_path, "skeleton.html")
  fixtures_dir <- normalizePath(testthat::test_path("fixtures"))
  withr::with_tempfile(output, {
    rmarkdown::render(input, params=list(reportType="webApp", outPlateEffect=FALSE, dataDir=fixtures_dir, xmlFiles=files ))  
    expect_true(file.exists(output), output) 
  })
})
test_that("skeleton.Rmd can be rendered into HTML, NAS version, advancedQC", { 
  test_path <- paste0(testthat::test_path(), "./../inst/rmarkdown/templates/nulisaseq/skeleton/")
  input <- paste0(test_path, "skeleton.Rmd")
  files <- c("Analysis_INF250_Lot4_AQ_LC_R3_20241229.xml", "Analysis_LC_R3_INF250_Lot3_20240727.xml")
  output <- paste0(test_path, "skeleton.html")
  fixtures_dir <- normalizePath(testthat::test_path("fixtures"))
  withr::with_tempfile(output, {
    rmarkdown::render(input, params=list(reportType="webApp", outPlateEffect=TRUE, advancedQC=TRUE, dataDir=fixtures_dir, xmlFiles=files ))  
    expect_true(file.exists(output), output) 
  })
})
test_that("skeleton.Rmd can be rendered into HTML, NAS version, XML_v1.3.0", {
  test_path <- paste0(testthat::test_path(), "./../inst/rmarkdown/templates/nulisaseq/skeleton/")
  input <- paste0(test_path, "skeleton.Rmd")
  files <- c("XML_v1.3.0.xml")
  output <- paste0(test_path, "skeleton.html")
  fixtures_dir <- normalizePath(testthat::test_path("fixtures"))
  withr::with_tempfile(output, {
    rmarkdown::render(input, params=list(reportType="webApp", outPlateEffect=TRUE, advancedQC=FALSE, dataDir=fixtures_dir, xmlFiles=files ))
    expect_true(file.exists(output), output)
  })
})

test_that("skeleton.Rmd with outputPlots=TRUE generates PDF files (WITH NULISAseqAQ)", {
  skip_if_not(requireNamespace("NULISAseqAQ", quietly = TRUE),
              "NULISAseqAQ package not available")
  
  input <- system.file("rmarkdown", "templates", "nulisaseq", "skeleton", "skeleton.Rmd",
                       package = "NULISAseqR")
  files <- c("XML_v1.3.0_with_AQ.xml")
  
  fixtures_dir <- normalizePath(testthat::test_path("fixtures"))
  
  # Use temp directory - make sure it's an absolute path
  temp_output_dir <- normalizePath(tempfile(), mustWork = FALSE)
  dir.create(temp_output_dir)
  output_files_dir <- file.path(temp_output_dir, "outputFiles")
  
  # Render to temp location
  temp_html <- file.path(temp_output_dir, "skeleton.html")
  
  rmarkdown::render(
    input, 
    output_file = temp_html,
    params = list(
      dataDir = fixtures_dir, 
      xmlFiles = files, 
      outputPlots = TRUE,
      reportType = "internal",
      outfolder = temp_output_dir  # This is now an absolute path
    )
  )
  
  # Check HTML output exists
  expect_true(file.exists(temp_html), info = "HTML output should be created")
  
  # Check outputFiles directory was created
  expect_true(dir.exists(output_files_dir), info = "outputFiles directory should be created")
  
  # Check that PDF files were generated
  pdf_files <- list.files(output_files_dir, pattern = "\\.pdf$", full.names = TRUE)
  expect_gt(length(pdf_files), 0, label = "Number of PDF files generated")
  
  # Clean up
  unlink(temp_output_dir, recursive = TRUE)
})

# Helper function to mock requireNamespace for NULISAseqAQ
# Uses local() to avoid global state leakage between tests
mock_no_NULISAseqAQ_skeleton <- local({
  original_requireNamespace <- base::requireNamespace
  function(package, ...) {
    if (package == "NULISAseqAQ") return(FALSE)
    original_requireNamespace(package, ...)
  }
})

test_that("skeleton.Rmd with outputPlots=TRUE generates PDF files (WITHOUT NULISAseqAQ)", {
  input <- system.file("rmarkdown", "templates", "nulisaseq", "skeleton", "skeleton.Rmd",
                       package = "NULISAseqR")
  files <- c("XML_v1.3.0_with_AQ.xml")
  
  fixtures_dir <- normalizePath(testthat::test_path("fixtures"))
  
  # Use temp directory - make sure it's an absolute path
  temp_output_dir <- normalizePath(tempfile(), mustWork = FALSE)
  dir.create(temp_output_dir)
  output_files_dir <- file.path(temp_output_dir, "outputFiles")
  
  # Render to temp location
  temp_html <- file.path(temp_output_dir, "skeleton.html")
  
  # Mock requireNamespace to return FALSE for NULISAseqAQ
  local_mocked_bindings(
    requireNamespace = mock_no_NULISAseqAQ_skeleton,
    .package = "base"
  )
  
  rmarkdown::render(
    input, 
    output_file = temp_html,
    params = list(
      dataDir = fixtures_dir, 
      xmlFiles = files, 
      outputPlots = TRUE,
      reportType = "internal",
      outfolder = temp_output_dir  # Now an absolute path
    )
  )
  
  # Check HTML output exists
  expect_true(file.exists(temp_html), info = "HTML output should be created")
  
  # Check outputFiles directory was created
  expect_true(dir.exists(output_files_dir), info = "outputFiles directory should be created")
  
  # Check that PDF files were generated
  pdf_files <- list.files(output_files_dir, pattern = "\\.pdf$", full.names = TRUE)
  expect_gt(length(pdf_files), 0, label = "Number of PDF files generated")
  
  # Clean up
  unlink(temp_output_dir, recursive = TRUE)
})

test_that("XML_v1.3.0.xml loads correctly WITH NULISAseqAQ", {
  skip_if_not(requireNamespace("NULISAseqAQ", quietly = TRUE),
              "NULISAseqAQ package not available")

  input1 <- test_path("fixtures", "XML_v1.3.0.xml")

  data <- loadNULISAseq(input1, IPC=NULL, IC='mCherry', SC=NULL)

  # Basic structure checks
  expect_false(is.null(data$qcXML))
  expect_false(is.null(data$AQ))
  expect_true("Data_AQ_aM" %in% names(data$AQ))

  # Check AQ structure when NULISAseqAQ is available
  expect_true("targetAQ_param" %in% names(data$AQ))
  expect_true("withinDR" %in% names(data$AQ))

  # Check quantifiability was calculated
  expect_false(is.null(data$quantifiability))
})

test_that("XML_v1.3.0.xml loads correctly WITHOUT NULISAseqAQ (fallback mode)", {
  input1 <- test_path("fixtures", "XML_v1.3.0.xml")

  # Mock requireNamespace to return FALSE for NULISAseqAQ
  local_mocked_bindings(
    requireNamespace = mock_no_NULISAseqAQ_skeleton,
    .package = "base"
  )

  # Should see the fallback message
  expect_message(
    data <- loadNULISAseq(input1, IPC=NULL, IC='mCherry', SC=NULL),
    "LLOQ/ULOQ not calculated - using values from XML"
  )

  # Basic structure checks
  expect_false(is.null(data$qcXML))
  expect_false(is.null(data$AQ))

  # Check AQ structure in fallback mode
  expect_true("targetAQ_param" %in% names(data$AQ))
  expect_true("Data_AQ_aM" %in% names(data$AQ))
  expect_true("withinDR" %in% names(data$AQ))

  # Check quantifiability was calculated (uses LOD_aM populated from qcXML in fallback mode)
  expect_false(is.null(data$quantifiability))

  # Check LOD_aM and LOD_pgmL are populated from qcXML in fallback mode
  expect_false(is.null(data$lod$LOD_aM),
               info = "LOD_aM should be populated from qcXML in fallback mode")
  expect_false(is.null(data$lod$LOD_pgmL),
               info = "LOD_pgmL should be populated from qcXML in fallback mode")

  # Check LOD vectors have names (target names)
  expect_true(length(names(data$lod$LOD_aM)) > 0,
              info = "LOD_aM should be a named vector")
  expect_true(length(names(data$lod$LOD_pgmL)) > 0,
              info = "LOD_pgmL should be a named vector")
})

test_that("XML_v1.3.0.xml AQ data consistency between NULISAseqAQ and fallback modes", {
  skip_if_not(requireNamespace("NULISAseqAQ", quietly = TRUE),
              "NULISAseqAQ package not available - cannot compare modes")

  input1 <- test_path("fixtures", "XML_v1.3.0.xml")

  # Load with NULISAseqAQ
  data_with_aq <- loadNULISAseq(input1, IPC=NULL, IC='mCherry', SC=NULL)

  # Load without NULISAseqAQ (mocked)
  local_mocked_bindings(
    requireNamespace = mock_no_NULISAseqAQ_skeleton,
    .package = "base"
  )

  suppressMessages({
    data_without_aq <- loadNULISAseq(input1, IPC=NULL, IC='mCherry', SC=NULL)
  })

  # Both should have the same basic structure (same named elements)
  expect_setequal(names(data_with_aq$AQ), names(data_without_aq$AQ))

  # Both should have Data_AQ_aM with identical dimensions
  if (!is.null(data_with_aq$AQ$Data_AQ_aM) && !is.null(data_without_aq$AQ$Data_AQ_aM)) {
    expect_equal(ncol(data_with_aq$AQ$Data_AQ_aM), ncol(data_without_aq$AQ$Data_AQ_aM))
    expect_equal(nrow(data_with_aq$AQ$Data_AQ_aM), nrow(data_without_aq$AQ$Data_AQ_aM))
    # Both should have the same targets
    expect_setequal(rownames(data_with_aq$AQ$Data_AQ_aM), rownames(data_without_aq$AQ$Data_AQ_aM))
  }

  # Both should have targetAQ_param with identical row counts
  expect_false(is.null(data_with_aq$AQ$targetAQ_param))
  expect_false(is.null(data_without_aq$AQ$targetAQ_param))
  expect_equal(nrow(data_with_aq$AQ$targetAQ_param), nrow(data_without_aq$AQ$targetAQ_param))
})

# Tests for XML_v1.3.0_with_AQ.xml (has pre-computed aM values)
test_that("XML_v1.3.0_with_AQ.xml loads correctly WITH NULISAseqAQ", {
  skip_if_not(requireNamespace("NULISAseqAQ", quietly = TRUE),
              "NULISAseqAQ package not available")

  input1 <- test_path("fixtures", "XML_v1.3.0_with_AQ.xml")

  data <- loadNULISAseq(input1, IPC=NULL, IC='mCherry', SC=NULL)

  # Basic structure checks
  expect_false(is.null(data$qcXML))
  expect_false(is.null(data$AQ))
  expect_true("Data_AQ_aM" %in% names(data$AQ))

  # Check AQ structure when NULISAseqAQ is available
  expect_true("targetAQ_param" %in% names(data$AQ))
  expect_true("withinDR" %in% names(data$AQ))

  # Check quantifiability was calculated
  expect_false(is.null(data$quantifiability))

  # Should have 157 targets (matching encrypted params in ExecutionDetails$Abs)
  expect_equal(nrow(data$AQ$Data_AQ_aM), 157)
})

test_that("XML_v1.3.0_with_AQ.xml loads correctly WITHOUT NULISAseqAQ (fallback mode)", {
  input1 <- test_path("fixtures", "XML_v1.3.0_with_AQ.xml")

  # Mock requireNamespace to return FALSE for NULISAseqAQ
  local_mocked_bindings(
    requireNamespace = mock_no_NULISAseqAQ_skeleton,
    .package = "base"
  )

  # Should see the fallback message
  expect_message(
    data <- loadNULISAseq(input1, IPC=NULL, IC='mCherry', SC=NULL),
    "LLOQ/ULOQ not calculated - using values from XML"
  )

  # Basic structure checks
  expect_false(is.null(data$qcXML))
  expect_false(is.null(data$AQ))

  # Check AQ structure in fallback mode
  expect_true("targetAQ_param" %in% names(data$AQ))
  expect_true("Data_AQ_aM" %in% names(data$AQ))
  expect_true("withinDR" %in% names(data$AQ))

  # Should have 157 targets (filtered to match encrypted params)
  expect_equal(nrow(data$AQ$Data_AQ_aM), 157)

  # Check LOD_aM and LOD_pgmL are populated from qcXML in fallback mode
  expect_false(is.null(data$lod$LOD_aM),
               info = "LOD_aM should be populated from qcXML in fallback mode")
  expect_false(is.null(data$lod$LOD_pgmL),
               info = "LOD_pgmL should be populated from qcXML in fallback mode")
})

test_that("process_loadNULISAseq creates aqParams in fallback mode (without MW_kDa)", {
  input1 <- test_path("fixtures", "XML_v1.3.0_with_AQ.xml")

  # Mock requireNamespace to return FALSE for NULISAseqAQ
  local_mocked_bindings(
    requireNamespace = mock_no_NULISAseqAQ_skeleton,
    .package = "base"
  )

  suppressMessages({
    data <- loadNULISAseq(input1, IPC=NULL, IC='mCherry', SC=NULL)
  })

  # Process the data
  processed <- NULISAseqR:::process_loadNULISAseq(data)

  # aqParams should be created even without MW_kDa
  expect_false(is.null(processed$aqParams),
               info = "aqParams should be created from targetAQ_param even without MW_kDa")

  # aqParams should have targetName column
  expect_true("targetName" %in% colnames(processed$aqParams),
              info = "aqParams should have targetName column")

  # Check Data_AQ exists (aM units)
  expect_false(is.null(processed$Data_AQ),
               info = "Data_AQ should exist")

  # Check Data_AQ_pgmL is created from pre-existing Data_AQ if available
  if (!is.null(data$AQ$Data_AQ)) {
    expect_false(is.null(processed$Data_AQ_pgmL),
                 info = "Data_AQ_pgmL should be created from pre-existing Data_AQ")
    expect_false(is.null(processed$Data_AQlog2_pgmL),
                 info = "Data_AQlog2_pgmL should be created")
  }
})

test_that("XML_v1.3.0_with_AQ.xml AQ data consistency between NULISAseqAQ and fallback modes", {
  skip_if_not(requireNamespace("NULISAseqAQ", quietly = TRUE),
              "NULISAseqAQ package not available - cannot compare modes")

  input1 <- test_path("fixtures", "XML_v1.3.0_with_AQ.xml")

  # Load with NULISAseqAQ
  data_with_aq <- loadNULISAseq(input1, IPC=NULL, IC='mCherry', SC=NULL)

  # Load without NULISAseqAQ (mocked)
  local_mocked_bindings(
    requireNamespace = mock_no_NULISAseqAQ_skeleton,
    .package = "base"
  )

  suppressMessages({
    data_without_aq <- loadNULISAseq(input1, IPC=NULL, IC='mCherry', SC=NULL)
  })

  # Both should have the same basic structure (same named elements)
  expect_setequal(names(data_with_aq$AQ), names(data_without_aq$AQ))

  # Both should have Data_AQ_aM with identical dimensions
  expect_equal(ncol(data_with_aq$AQ$Data_AQ_aM), ncol(data_without_aq$AQ$Data_AQ_aM))
  expect_equal(nrow(data_with_aq$AQ$Data_AQ_aM), nrow(data_without_aq$AQ$Data_AQ_aM))
  # Both should have the same targets
  expect_setequal(rownames(data_with_aq$AQ$Data_AQ_aM), rownames(data_without_aq$AQ$Data_AQ_aM))

  # Data values should match where fallback has values
  # Note: Fallback may have additional NAs for targets where XML has aM="nan"
  common_rows <- intersect(rownames(data_with_aq$AQ$Data_AQ_aM),
                           rownames(data_without_aq$AQ$Data_AQ_aM))

  with_aq <- data_with_aq$AQ$Data_AQ_aM[common_rows, ]
  without_aq <- data_without_aq$AQ$Data_AQ_aM[common_rows, ]

  # Where fallback has non-NA, non-zero values, they should match NULISAseqAQ (within tolerance)
  # Note: Fallback may have zeros where NULISAseqAQ calculates non-zero values because the
  # embedded XML values may be zero while NULISAseqAQ does full calculations from raw data
  # Use relative tolerance based on the magnitude of values for robustness across scales
  comparable_mask <- !is.na(without_aq) & !is.na(with_aq) & without_aq != 0 & with_aq != 0
  if (any(comparable_mask)) {
    max_diff <- max(abs(with_aq[comparable_mask] - without_aq[comparable_mask]), na.rm = TRUE)
    ref_vals <- abs(without_aq[comparable_mask])
    scale <- suppressWarnings(max(ref_vals, na.rm = TRUE))
    # Use relative tolerance (1e-6) scaled by max value
    rel_tol <- 1e-6
    tolerance <- if (!is.finite(scale) || scale == 0) rel_tol else rel_tol * scale
    expect_true(max_diff <= tolerance,
                info = paste("Max diff:", max_diff, "Tolerance:", tolerance))
  }

  # Both should have targetAQ_param with identical row counts
  expect_false(is.null(data_with_aq$AQ$targetAQ_param))
  expect_false(is.null(data_without_aq$AQ$targetAQ_param))
  expect_equal(nrow(data_with_aq$AQ$targetAQ_param), nrow(data_without_aq$AQ$targetAQ_param))
})

# =============================================================================
# Tests for SC_conc (sample control concentration) parsing
# =============================================================================

test_that("SC_conc values are parsed from XML TargetQC section (WITH NULISAseqAQ)", {
  skip_if_not(requireNamespace("NULISAseqAQ", quietly = TRUE),
              "NULISAseqAQ package not available")

  input1 <- test_path("fixtures", "XML_v1.3.0_with_AQ.xml")

  data <- loadNULISAseq(input1, IPC=NULL, IC='mCherry', SC=NULL)

  # Check that SC_conc was parsed into qcXML
  expect_false(is.null(data$qcXML$target$SC_conc),
               info = "qcXML$target$SC_conc should be populated from XML")

  # Check SC_conc has expected columns (renamed by readQCXMLNode)
  if (!is.null(data$qcXML$target$SC_conc)) {
    expect_true("SC_conc_aM" %in% colnames(data$qcXML$target$SC_conc) ||
                "target" %in% colnames(data$qcXML$target$SC_conc),
                info = "SC_conc should have SC_conc_aM or target column")
  }

  # Check that SC_conc values are merged into targetAQ_param
  expect_false(is.null(data$AQ$targetAQ_param),
               info = "targetAQ_param should exist")

  if (!is.null(data$AQ$targetAQ_param)) {
    # SC_conc_aM should be present if XML contains SC_conc data
    sc_cols <- grep("SC_conc", colnames(data$AQ$targetAQ_param), value = TRUE)
    expect_true(length(sc_cols) > 0,
                info = "targetAQ_param should contain SC_conc columns")
  }
})

test_that("SC_conc values are parsed from XML TargetQC section (WITHOUT NULISAseqAQ)", {
  input1 <- test_path("fixtures", "XML_v1.3.0_with_AQ.xml")

  # Mock requireNamespace to return FALSE for NULISAseqAQ
  local_mocked_bindings(
    requireNamespace = mock_no_NULISAseqAQ_skeleton,
    .package = "base"
  )

  suppressMessages({
    data <- loadNULISAseq(input1, IPC=NULL, IC='mCherry', SC=NULL)
  })

  # Check that SC_conc was parsed into qcXML (fallback mode)
  expect_false(is.null(data$qcXML$target$SC_conc),
               info = "qcXML$target$SC_conc should be populated from XML in fallback mode")

  # Check that SC_conc values are merged into targetAQ_param
  expect_false(is.null(data$AQ$targetAQ_param),
               info = "targetAQ_param should exist in fallback mode")

  if (!is.null(data$AQ$targetAQ_param)) {
    # SC_conc_aM should be present if XML contains SC_conc data
    sc_cols <- grep("SC_conc", colnames(data$AQ$targetAQ_param), value = TRUE)
    expect_true(length(sc_cols) > 0,
                info = "targetAQ_param should contain SC_conc columns in fallback mode")
  }
})

test_that("SC_conc values are numeric in targetAQ_param", {
  skip_if_not(requireNamespace("NULISAseqAQ", quietly = TRUE),
              "NULISAseqAQ package not available")

  input1 <- test_path("fixtures", "XML_v1.3.0_with_AQ.xml")

  data <- loadNULISAseq(input1, IPC=NULL, IC='mCherry', SC=NULL)

  expect_false(is.null(data$AQ$targetAQ_param),
               info = "targetAQ_param should exist")

  # Find all SC_conc-related columns
  sc_cols <- grep("SC_conc", colnames(data$AQ$targetAQ_param), value = TRUE)
  expect_true(length(sc_cols) > 0,
              info = "targetAQ_param should have SC_conc columns")

  # Check each SC_conc column is numeric
  for (col in sc_cols) {
    expect_true(is.numeric(data$AQ$targetAQ_param[[col]]),
                info = paste(col, "should be numeric"))
  }
})

test_that("SC_conc columns present in both NULISAseqAQ and fallback modes", {
  skip_if_not(requireNamespace("NULISAseqAQ", quietly = TRUE),
              "NULISAseqAQ package not available - cannot compare modes")

  input1 <- test_path("fixtures", "XML_v1.3.0_with_AQ.xml")

  # Load with NULISAseqAQ
  data_with_aq <- loadNULISAseq(input1, IPC=NULL, IC='mCherry', SC=NULL)

  # Load without NULISAseqAQ (mocked)
  local_mocked_bindings(
    requireNamespace = mock_no_NULISAseqAQ_skeleton,
    .package = "base"
  )

  suppressMessages({
    data_without_aq <- loadNULISAseq(input1, IPC=NULL, IC='mCherry', SC=NULL)
  })

  # Both modes should have at least one SC_conc-related column in targetAQ_param
  with_sc_cols <- grep("SC_conc", colnames(data_with_aq$AQ$targetAQ_param), value = TRUE)
  without_sc_cols <- grep("SC_conc", colnames(data_without_aq$AQ$targetAQ_param), value = TRUE)

  expect_true(length(with_sc_cols) > 0,
              info = "WITH NULISAseqAQ should have SC_conc columns")
  expect_true(length(without_sc_cols) > 0,
              info = "WITHOUT NULISAseqAQ should have SC_conc columns")

  # Both should have qcXML SC_conc parsed from XML
  expect_false(is.null(data_with_aq$qcXML$target$SC_conc),
               info = "WITH mode should have qcXML$target$SC_conc")
  expect_false(is.null(data_without_aq$qcXML$target$SC_conc),
               info = "WITHOUT mode should have qcXML$target$SC_conc")

  # qcXML SC_conc should be identical (same XML source)
  expect_equal(nrow(data_with_aq$qcXML$target$SC_conc),
               nrow(data_without_aq$qcXML$target$SC_conc),
               info = "qcXML SC_conc row count should match between modes")
})
