# Define a test processXML
test_that("Test that loadNULISAseq reads in an XML file and adds lists of appropriate size", {
  
  test_path <- paste0(testthat::test_path(), "./../inst/rmarkdown/templates/nulisaseq/skeleton/")

  input1 <- paste0(test_path, "detectability_P1_Tr03.xml")
  data <- loadNULISAseq(input1, IPC=NULL, IC='mCherry', SC=NULL)
  expect_true(nrow(data$qcPlate) >= 5)
  expect_true(nrow(data$qcSample) > ncol(data$Data))
  expect_true(all(dim(data$IC_normed$normData) == dim(data$normed$interNormData)))
  expect_true(all(dim(data$IC_normed$normData) == dim(data$normed$log2_interNormData)))
  expect_true(all(dim(data$normed$normData) == dim(data$normed$interNormData)))

})
