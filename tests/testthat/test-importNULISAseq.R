# Helper function to mock requireNamespace for NULISAseqAQ
# Uses local() to avoid global state leakage between tests
mock_no_NULISAseqAQ_skeleton <- local({
  original_requireNamespace <- base::requireNamespace
  function(package, ...) {
    if (package == "NULISAseqAQ") return(FALSE)
    original_requireNamespace(package, ...)
  }
})

test_that("importNULISAseq reads XML with AQ data WITH NULISAseqAQ", {
  skip_if_not(requireNamespace("NULISAseqAQ", quietly = TRUE),
              "NULISAseqAQ package not available")
  
  input_file <- test_path("fixtures", "XML_v1.3.0_with_AQ.xml")
  
  data <- importNULISAseq(files = input_file)
  
  # Basic structure checks
  expect_false(is.null(data$merged$Data_NPQ_long))
  expect_false(is.null(data$merged$Data_AQ_long))
  
  # Check AQ structure when NULISAseqAQ is available
  expect_true("Data_AQ_aM" %in% names(data$merged))
  expect_true("Data_AQlog2_aM" %in% names(data$merged))
  expect_true("Data_AQ_pgmL" %in% names(data$merged))
  expect_true("Data_AQlog2_pgmL" %in% names(data$merged))
  
  # Check quantifiability was calculated
  expect_false(is.null(data$merged$quantifiability))
  
  # Should have 157 targets
  expect_equal(nrow(data$merged$Data_AQ_long), 15072)
})

test_that("importNULISAseq reads XML with AQ data WITHOUT NULISAseqAQ (fallback mode)", {
  input_file <- test_path("fixtures", "XML_v1.3.0_with_AQ.xml")
  
  # Mock requireNamespace to return FALSE for NULISAseqAQ
  local_mocked_bindings(
    requireNamespace = mock_no_NULISAseqAQ_skeleton,
    .package = "base"
  )
  
  expect_no_error(
    data <- importNULISAseq(files = input_file)
  )
  
  # Basic structure checks
  expect_false(is.null(data$merged$Data_NPQ_long))
  expect_false(is.null(data$merged$Data_AQ_long))
  
  # Check AQ structure in fallback mode
  expect_true("Data_AQ_aM" %in% names(data$merged))
  expect_true("Data_AQlog2_aM" %in% names(data$merged))
  expect_true("Data_AQ_pgmL" %in% names(data$merged))
  expect_true("Data_AQlog2_pgmL" %in% names(data$merged))
  
  # Should have 15072 rows (157 targets × 96 samples)
  expect_equal(nrow(data$merged$Data_AQ_long), 15072)
  
  # Check quantifiability was calculated using fallback LOD_aM
  expect_false(is.null(data$merged$quantifiability))
})

test_that("importNULISAseq AQ data consistency between NULISAseqAQ and fallback modes", {
  skip_if_not(requireNamespace("NULISAseqAQ", quietly = TRUE),
              "NULISAseqAQ package not available - cannot compare modes")
  
  input_file <- test_path("fixtures", "XML_v1.3.0_with_AQ.xml")
  
  # Load with NULISAseqAQ
  data_with_aq <- importNULISAseq(files = input_file)
  
  # Load without NULISAseqAQ (mocked)
  local_mocked_bindings(
    requireNamespace = mock_no_NULISAseqAQ_skeleton,
    .package = "base"
  )
  
  suppressMessages({
    data_without_aq <- importNULISAseq(files = input_file)
  })
  
  # Both should have the same basic structure
  expect_true(all(c("Data_AQ_aM", "Data_AQlog2_aM", "Data_AQ_pgmL", "Data_AQlog2_pgmL") %in% 
                    names(data_with_aq$merged)))
  expect_true(all(c("Data_AQ_aM", "Data_AQlog2_aM", "Data_AQ_pgmL", "Data_AQlog2_pgmL") %in% 
                    names(data_without_aq$merged)))
  
  # Both should have Data_AQ_long with identical dimensions
  expect_equal(nrow(data_with_aq$merged$Data_AQ_long), 
               nrow(data_without_aq$merged$Data_AQ_long))
  expect_equal(ncol(data_with_aq$merged$Data_AQ_long), 
               ncol(data_without_aq$merged$Data_AQ_long))
  
  # Both should have the same targets
  expect_setequal(unique(data_with_aq$merged$Data_AQ_long$Target), 
                  unique(data_without_aq$merged$Data_AQ_long$Target))
})

# Tests for XML_v1.3.0.xml (without pre-computed AQ data)
test_that("importNULISAseq reads XML_v1.3.0.xml WITH NULISAseqAQ", {
  skip_if_not(requireNamespace("NULISAseqAQ", quietly = TRUE),
              "NULISAseqAQ package not available")
  
  input_file <- test_path("fixtures", "XML_v1.3.0.xml")
  
  data <- importNULISAseq(files = input_file)
  
  # Basic structure checks
  expect_false(is.null(data$merged$Data_NPQ_long))
  expect_false(is.null(data$merged$Data_AQ_long))
  
  # Check AQ structure when NULISAseqAQ is available
  expect_true("Data_AQ_aM" %in% names(data$merged))
  expect_true("Data_AQlog2_aM" %in% names(data$merged))
  expect_true("Data_AQ_pgmL" %in% names(data$merged))
  expect_true("Data_AQlog2_pgmL" %in% names(data$merged))
  
  # Check quantifiability was calculated
  expect_false(is.null(data$merged$quantifiability))
})

test_that("importNULISAseq reads XML_v1.3.0.xml WITHOUT NULISAseqAQ (fallback mode)", {
  input_file <- test_path("fixtures", "XML_v1.3.0.xml")
  
  # Mock requireNamespace to return FALSE for NULISAseqAQ
  local_mocked_bindings(
    requireNamespace = mock_no_NULISAseqAQ_skeleton,
    .package = "base"
  )
  
  # Should see the fallback message
  expect_message(
    data <- importNULISAseq(files = input_file),
    "LLOQ/ULOQ not calculated - using values from XML"
  )
  
  # Basic structure checks
  expect_false(is.null(data$merged$Data_NPQ_long))
  expect_false(is.null(data$merged$Data_AQ_long))
  
  # Check AQ structure in fallback mode
  expect_true("Data_AQ_aM" %in% names(data$merged))
  expect_true("Data_AQlog2_aM" %in% names(data$merged))
  expect_true("Data_AQ_pgmL" %in% names(data$merged))
  expect_true("Data_AQlog2_pgmL" %in% names(data$merged))
  
  # Check quantifiability was calculated using fallback LOD_aM
  expect_false(is.null(data$merged$quantifiability))
})

test_that("importNULISAseq XML_v1.3.0_with_AQ.xml AQ data consistency between modes", {
  skip_if_not(requireNamespace("NULISAseqAQ", quietly = TRUE),
              "NULISAseqAQ package not available - cannot compare modes")
  
  input_file <- test_path("fixtures", "XML_v1.3.0_with_AQ.xml")
  
  # Load with NULISAseqAQ
  data_with_aq <- importNULISAseq(files = input_file)
  
  # Load without NULISAseqAQ (mocked)
  local_mocked_bindings(
    requireNamespace = mock_no_NULISAseqAQ_skeleton,
    .package = "base"
  )
  
  suppressMessages({
    data_without_aq <- importNULISAseq(files = input_file)
  })
  
  # Both should have the same basic structure
  expect_true(all(c("Data_AQ_aM", "Data_AQlog2_aM", "Data_AQ_pgmL", "Data_AQlog2_pgmL") %in% 
                    names(data_with_aq$merged)))
  expect_true(all(c("Data_AQ_aM", "Data_AQlog2_aM", "Data_AQ_pgmL", "Data_AQlog2_pgmL") %in% 
                    names(data_without_aq$merged)))
  
  # Both should have Data_AQ_long with identical dimensions
  expect_equal(nrow(data_with_aq$merged$Data_AQ_long), 
               nrow(data_without_aq$merged$Data_AQ_long))
  expect_equal(ncol(data_with_aq$merged$Data_AQ_long), 
               ncol(data_without_aq$merged$Data_AQ_long))
  
  # Both should have the same targets
  expect_setequal(unique(data_with_aq$merged$Data_AQ_long$Target), 
                  unique(data_without_aq$merged$Data_AQ_long$Target))
})