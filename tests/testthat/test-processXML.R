# Define a test for QCSampleCriteria
test_that("", {
  
  test_path <- paste0(testthat::test_path(), "./../inst/rmarkdown/templates/nulisaseq/skeleton/")
  schema_file <- paste0(testthat::test_path(), "./../NGSNormalize.xsd")

  input1 <- paste0(test_path, "20230325_TAP_Qatar_plate01_no_nonmatching.xml")
  input2 <- paste0(test_path, "20230325_TAP_Qatar_plate02_no_nonmatching_XML.xml")
  val <- withr::with_tempfile("file1", {
    out <- processXML(in_xml=input1, out_XML=file1)
    schema <- xml2::read_xml(schema_file, package="xml2")
    input <- xml2::read_xml(file1, package="xml2")
    xml_validate(input, schema)
  })
  expect_equal(val[1], TRUE)
})
