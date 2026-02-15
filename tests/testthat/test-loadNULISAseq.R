# Define a test processXML
test_that("loadNULISAseq reads in an XML file and adds lists of appropriate size", {
  
  input1 <- test_path("fixtures", "detectability_P1_Tr03.xml")
  data <- loadNULISAseq(input1, IPC=NULL, IC='mCherry', SC=NULL)
  expect_true(nrow(data$qcPlate) >= 5)
  expect_true(nrow(data$qcSample) > ncol(data$Data))
  expect_true(all(dim(data$IC_normed$normData) == dim(data$normed$interNormData)))
  expect_true(all(dim(data$IC_normed$normData) == dim(data$normed$log2_interNormData)))
  expect_true(all(dim(data$normed$normData) == dim(data$normed$interNormData)))
  
})

test_that("loadNULISAseq can handle v1.3.x version of the XML which incorporates QC and NPQ /AQ values",{
  input1 <- test_path("fixtures", "XML_v1.3.0.xml")
  data <- loadNULISAseq(input1, IPC=NULL, IC='mCherry', SC=NULL)
  expect_false(is.null(data$qcXML))
})

test_that("loadNULISAseq returns expected forward and reverse curve target values", {
  
  input1 <- test_path("fixtures", "Analysis_INF250_Lot4_AQ_LC_R3_20241229.xml")
  data <- loadNULISAseq(input1)
  
  # Spot-check specific reverse target NPQ values
  reverse_target <- "CRP"
  expected_reverse_NPQ <- c(25.0840257805487, 24.6552950934552, 23.7384333981208, 
                            24.3538798566515, 10.6233775377862, 9.13682978406353, 
                            9.04802404141014, 10.6394783856514, 25.0198117002621, 
                            14.232003081347)
  observed_reverse_NPQ <- data$NPQ[reverse_target, 1:10]
  
  # Check that target is marked as reverse curve
  expect_equal(data$targets$Curve_Quant[data$targets$targetName == reverse_target], 'R',
               label = sprintf("%s should be marked as reverse curve target", reverse_target))
  
  # Check reverse curve NPQ values match expected (convert to numeric first)
  expect_equal(as.numeric(observed_reverse_NPQ), expected_reverse_NPQ, tolerance = 1e-4,
               label = sprintf("NPQ values for reverse target %s", reverse_target))
  
  # Check that LOD is NA for reverse curve targets
  expect_true(is.na(data$lod$LOD[data$targets$targetName == reverse_target]),
              info = sprintf("%s should have NA for LOD", reverse_target))
  
  # Spot-check specific forward (normal) target NPQ values
  forward_target <- "AGER"  
  expected_forward_NPQ <- c(0, 2.87205556533307, 0, 
                            0, 7.85584349234945, 3.29032128894335, 
                            5.20967166909267, 0, 0, 
                            11.2948882462736) 
  observed_forward_NPQ <- data$NPQ[forward_target, 1:10]
  
  # Check that target is marked as forward curve 
  expect_true(data$targets$Curve_Quant[data$targets$targetName == forward_target] == 'F',
              label = sprintf("%s should be marked as forward curve quant target", forward_target))
  
  # Check forward curve NPQ values match expected
  expect_equal(as.numeric(observed_forward_NPQ), expected_forward_NPQ, tolerance = 1e-4,
               label = sprintf("NPQ values for forward target %s", forward_target))
  
  # Check that LOD exists for forward curve targets and is the correct value
  # use target with non-zero LOD
  forward_target <- "AGRP"  
  forward_lod <- data$lod$LOD[data$targets$targetName == forward_target]
  expect_false(is.na(forward_lod),
               info = sprintf("%s should have LOD value", forward_target))
  
  # Check the actual LOD value matches expected
  expected_forward_LOD <- 103.14203914238
  expect_equal(as.numeric(forward_lod), expected_forward_LOD, tolerance = 1e-4,
               label = sprintf("LOD for forward target %s", forward_target))
  
  # check LOD for a target that had outlier removal
  forward_target <- "CCL3"  
  forward_lod <- data$lod$LOD[data$targets$targetName == forward_target]
  expect_false(is.na(forward_lod),
               info = sprintf("%s should have LOD value", forward_target))
  
  # Check the actual LOD value matches expected
  expected_forward_LOD <- 207.48118071394
  expect_equal(as.numeric(forward_lod), expected_forward_LOD, tolerance = 1e-4,
               label = sprintf("LOD for forward target %s", forward_target))
})