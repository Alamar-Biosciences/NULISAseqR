# Helper function to lazily load test data
# Only loads when called within a test, allowing skip_if_not to work properly
local_test_data <- function() {
  xml_file <- test_path("fixtures", "Analysis_INF250_Lot4_AQ_LC_R3_20241229.xml")
  skip_if_not(file.exists(xml_file), "Test fixture not available")

  list(
    xml_file = xml_file,
    data_load = loadNULISAseq(file = xml_file),
    data_import = importNULISAseq(files = xml_file)
  )
}

test_that("Reverse curve targets show expected correlations in loadNULISAseq", {
  test_data <- local_test_data()

  reverse_targets <- test_data$data_load$targets$targetName[
    test_data$data_load$targets$Curve_Quant == 'R'
  ]
  skip_if(length(reverse_targets) == 0, "No reverse curve targets in test data")
  
  # Test first reverse target
  target <- reverse_targets[1]
  npq <- test_data$data_load$NPQ[target, ]
  
  # Pre-transformation: negative correlations (should be < -0.5 for strong negative)
  cor_data <- cor(log2(test_data$data_load$Data[target, ] + 1), npq, use = "complete.obs")
  expect_true(cor_data < 0, info = "Raw data should be negatively correlated with NPQ")
  expect_lt(cor_data, -0.5, label = "Correlation should be strongly negative")
  
  cor_ic <- cor(log2(test_data$data_load$IC_normed$normData[target, ] + 1), npq, use = "complete.obs")
  expect_true(cor_ic < 0, info = "IC normalized data should be negatively correlated with NPQ")
  
  # Post-transformation: positive correlations (should be > 0.5 for strong positive)
  cor_normed <- cor(
    log2(test_data$data_load$normed$interNormData[[1]][target, ] + 1), 
    npq, 
    use = "complete.obs"
  )
  expect_true(cor_normed > 0, info = "Transformed data should be positively correlated with NPQ")
  expect_gt(cor_normed, 0.5, label = "Correlation should be strongly positive")
  
  cor_log2 <- cor(
    test_data$data_load$normed$log2_interNormData[[1]][target, ], 
    npq, 
    use = "complete.obs"
  )
  expect_true(cor_log2 > 0, info = "Log2 transformed data should be positively correlated with NPQ")
})

test_that("Reverse curve transformation inverts the data correctly", {
  test_data <- local_test_data()

  reverse_targets <- test_data$data_load$targets$targetName[
    test_data$data_load$targets$Curve_Quant == 'R'
  ]
  skip_if(length(reverse_targets) == 0, "No reverse curve targets in test data")
  
  target <- reverse_targets[1]
  
  # Get untransformed and transformed data
  untransformed <- test_data$data_load$normed_untransformedReverse$interNormData[[1]][target, ]
  transformed <- test_data$data_load$normed$interNormData[[1]][target, ]
  
  # Check that transformation inverts: for reverse curves, lower raw counts should become higher transformed values
  # Find indices of min and max untransformed values
  min_idx <- which.min(untransformed)
  max_idx <- which.max(untransformed)
  
  # After transformation, the min should become relatively high and max should become relatively low
  expect_true(
    transformed[min_idx] > transformed[max_idx],
    info = "Lowest untransformed value should become higher than highest after reverse transformation"
  )
})

test_that("Reverse curve NPQ values are within expected range", {
  test_data <- local_test_data()

  reverse_targets <- test_data$data_load$targets$targetName[
    test_data$data_load$targets$Curve_Quant == 'R'
  ]
  skip_if(length(reverse_targets) == 0, "No reverse curve targets in test data")
  
  target <- reverse_targets[1]
  npq <- test_data$data_load$NPQ[target, ]
  
  # NPQ is on log2 scale, should be non-negative
  expect_true(all(npq >= 0, na.rm = TRUE),
              info = "NPQ values should be non-negative (log2 scale)")
  
  # Should have some variation (not all the same)
  expect_gt(sd(npq, na.rm = TRUE), 0,
            label = "NPQ values should have variation")
  
  # Sanity check: NPQ values should typically be in a reasonable range (e.g., < 40 for log2 scale)
  # This is a soft check - if it fails, might indicate data issues but not necessarily algorithm issues
  median_npq <- median(npq, na.rm = TRUE)
  expect_lt(median_npq, 40,
            label = "Median NPQ should be in typical range (soft check for data quality)")
})

test_that("loadNULISAseq and importNULISAseq produce identical NPQ values for reverse targets", {
  test_data <- local_test_data()

  reverse_targets <- test_data$data_load$targets$targetName[
    test_data$data_load$targets$Curve_Quant == 'R'
  ]
  skip_if(length(reverse_targets) == 0, "No reverse curve targets in test data")
  
  # Test first target
  target <- reverse_targets[1]
  
  # Get NPQ from load (wide format)
  npq_load <- test_data$data_load$NPQ[target, ]
  
  # Get NPQ from import (long format)
  import_subset <- test_data$data_import$merged$Data_NPQ_long[
    test_data$data_import$merged$Data_NPQ_long$Target == target,
  ]
  
  # Match by sample name
  npq_import <- import_subset$NPQ[match(names(npq_load), import_subset$SampleName)]
  
  expect_equal(as.numeric(npq_load), npq_import, tolerance = 1e-10,
               info = "NPQ values should match between loadNULISAseq and importNULISAseq")
})

test_that("Reverse curve targets show negative correlation between NPQ and raw counts", {
  test_data <- local_test_data()

  reverse_targets <- test_data$data_import$merged$targets$targetName[
    test_data$data_import$merged$targets$Curve_Quant=='R'
    ]
  skip_if(length(reverse_targets) == 0, "No reverse curve targets in test data")
  
  # Test first target
  data_reverse <- test_data$data_import$merged$Data_NPQ_long[
    test_data$data_import$merged$Data_NPQ_long$Target == reverse_targets[1], 
  ]
  
  cor_val <- cor(
    data_reverse$NPQ, 
    log2(data_reverse$UnnormalizedCount + 1), 
    use = "complete.obs"
  )
  
  expect_true(cor_val < 0,
              info = "Reverse targets should show negative correlation between NPQ and raw counts")
  expect_lt(cor_val, -0.5, label = "Correlation should be strongly negative")
})

test_that("Reverse curve specific NPQ values match expected results", {
  test_data <- local_test_data()

  reverse_targets <- test_data$data_load$targets$targetName[
    test_data$data_load$targets$Curve_Quant == 'R'
  ]
  skip_if(length(reverse_targets) == 0, "No reverse curve targets in test data")
  
  target <- reverse_targets[1]
  
  # Test specific NPQ values from loadNULISAseq
  npq_load <- test_data$data_load$NPQ[target, ]
  
  expected_values_load <- c(
    "LP_Bridge_2" = 25.0840257805487, 
    "LC_HP_Std4_3" = 24.6552950934552,  
    "LC_HP_Std2_2" = 23.7384333981208
  )
  
  # Check specific samples
  for (sample_name in names(expected_values_load)) {
    if (sample_name %in% names(npq_load)) {
      expect_equal(
        npq_load[sample_name], 
        expected_values_load[sample_name],
        tolerance = 1e-4,
        label = sprintf("NPQ for %s from loadNULISAseq", sample_name)
      )
    }
  }
  
  # Test specific NPQ values from importNULISAseq
  import_subset <- test_data$data_import$merged$Data_NPQ_long[
    test_data$data_import$merged$Data_NPQ_long$Target == target,
  ]
  
  expected_values_import <- data.frame(
    SampleName = c("LP_Bridge_2", "LC_HP_Std4_3", "LC_HP_Std2_2"),
    NPQ = c(25.0840257805487, 24.6552950934552, 23.7384333981208),
    stringsAsFactors = FALSE
  )
  
  # Check specific samples from import
  for (i in 1:nrow(expected_values_import)) {
    sample_name <- expected_values_import$SampleName[i]
    expected_npq <- expected_values_import$NPQ[i]
    
    actual_npq <- import_subset$NPQ[import_subset$SampleName == sample_name]
    
    if (length(actual_npq) > 0) {
      expect_equal(
        actual_npq,
        expected_npq,
        tolerance = 1e-4,
        label = sprintf("NPQ for %s from importNULISAseq", sample_name)
      )
    }
  }
})