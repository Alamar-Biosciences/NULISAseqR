# Define a test processXML
test_that("Test that processXML reads in an XML file and outputs an XML according to the schema", {
  
  test_path <- paste0(testthat::test_path(), "./../inst/rmarkdown/templates/nulisaseq/skeleton/")

  input1 <- paste0(test_path, "20230325_TAP_Qatar_plate02_no_nonmatching_XML.xml")
  data <- readNULISAseq(input1)
  covarFile <- paste0(test_path, "covar.txt")
  data2 <- readCovariateFile(covarFile, list(data))
  expect_true("var1" %in% names(data2[[1]]$numericCovariates))
  expect_true("var3" %in% names(data2[[1]]$numericCovariates)) 
  expect_true(unname(data2[[1]]$numericCovariates["var1"]))
  expect_false(unname(data2[[1]]$numericCovariates["var3"])) 
 
  #These shouldn't match at all because the covar.txt and the input XML do not have the same samples 
  input1 <- paste0(test_path, "20230325_TAP_Qatar_plate01_no_nonmatching.xml")
  data <- readNULISAseq(input1)
  data2 <- readCovariateFile(covarFile, list(data))
  expect_true("var1" %in% names(data2[[1]]$numericCovariates))
  expect_true("var3" %in% names(data2[[1]]$numericCovariates)) 
  expect_false(unname(data2[[1]]$numericCovariates["var1"]))
  expect_false(unname(data2[[1]]$numericCovariates["var3"]))

})
