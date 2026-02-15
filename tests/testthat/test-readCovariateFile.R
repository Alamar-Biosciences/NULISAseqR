# Define a test processXML
test_that("Test that processXML reads in an XML file and outputs an XML according to the schema", {

  input1 <- test_path("fixtures", "detectability_P1_Tr03.xml")
  data <- readNULISAseq(input1)
  covarFile <- test_path("fixtures", "covar.csv")
  data2 <- readCovariateFile(covarFile, list(data))
  expect_true("var1" %in% names(data2[[1]]$numericCovariates))
  expect_true("var3" %in% names(data2[[1]]$numericCovariates)) 
  expect_false(unname(data2[[1]]$numericCovariates["var1"]))
  expect_true(unname(data2[[1]]$numericCovariates["var3"])) 

  #These shouldn't match at all because the covar.txt and the input XML do not have the same samples
  input1 <- test_path("fixtures", "detectability_P2_Tr03.xml")
  data <- readNULISAseq(input1)
  data2 <- readCovariateFile(covarFile, list(data))
  expect_true("var1" %in% names(data2[[1]]$numericCovariates))
  expect_true("var3" %in% names(data2[[1]]$numericCovariates)) 
  expect_true(all(is.na(unname(data2[[1]]$samples$var1))))
  expect_true(all(is.na(unname(data2[[1]]$samples$var3))))

})
