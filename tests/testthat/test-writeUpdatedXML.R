# Define a test processXML
test_that("Test that we can turn a pre-1.3.0 XML into a 1.3.0 XML that includes QC, etc.", {
  
  test_path <- paste0(testthat::test_path(), "./../inst/rmarkdown/templates/nulisaseq/skeleton/")

  input1 <- paste0(test_path, "detectability_P1_Tr03.xml")
  W <- writeUpdatedXML(input1)
  output <- "output.xml"
  withr::with_tempfile(output, {
    write(W, output) 
    data2<-loadNULISAseq(output)
    expect_true(file.exists(output),output)
  })
})
test_that("Test that we can turn a pre-1.3.0 AQ XML into a 1.3.0 AQ XML that includes QC, etc.", {
  skip_if_not(requireNamespace("NULISAseqAQ", quietly = TRUE),
              "NULISAseqAQ package not available")

  test_path <- paste0(testthat::test_path(), "./../inst/rmarkdown/templates/nulisaseq/skeleton/")

  input1 <- paste0(test_path, "Analysis_INF250_Lot4_AQ_LC_R3_20241229.xml")
  W <- writeUpdatedXML(input1)
  output <- "output.xml"
  withr::with_tempfile(output, {
    write(W, output)
    data2<-loadNULISAseq(output)
    expect_true(file.exists(output),output)
  })
})

test_that("Pre-1.3.0 AQ XML warns and loads as RQ when NULISAseqAQ is not available", {
  # Mock NULISAseqAQ as unavailable
  local_mocked_bindings(
    requireNamespace = function(pkg, ...) {
      if (pkg == "NULISAseqAQ") return(FALSE)
      base::requireNamespace(pkg, ...)
    },
    .package = "base"
  )

  test_path <- paste0(testthat::test_path(), "./../inst/rmarkdown/templates/nulisaseq/skeleton/")
  input1 <- paste0(test_path, "Analysis_INF250_Lot4_AQ_LC_R3_20241229.xml")

  # Should warn but load successfully as RQ
  expect_warning(data <- loadNULISAseq(input1), "AQ metadata")

  # Verify it loaded successfully as RQ
  expect_false(is.null(data$Data))
  expect_true(is.null(data$AQ) || is.null(data$AQ$Data_AQ_aM))
})
