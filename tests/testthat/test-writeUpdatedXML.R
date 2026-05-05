# Define a test processXML
test_that("Test that we can turn a pre-1.3.0 XML into a 1.3.0 XML that includes QC, etc.", {

  input1 <- test_path("fixtures", "detectability_P1_Tr03.xml")
  W <- writeUpdatedXML(input1, data = suppressWarnings(loadNULISAseq(input1, allowMissingCurveQuant = TRUE)))
  output <- "output.xml"
  withr::with_tempfile(output, {
    write(W, output)
    data2<-suppressWarnings(loadNULISAseq(output, allowMissingCurveQuant=TRUE))
    expect_true(file.exists(output),output)
  })
})
test_that("Test that we can turn a pre-1.3.0 AQ XML into a 1.3.0 AQ XML that includes QC, etc.", {
  skip_if_not(requireNamespace("NULISAseqAQ", quietly = TRUE),
              "NULISAseqAQ package not available")

  input1 <- test_path("fixtures", "Analysis_INF250_Lot4_AQ_LC_R3_20241229.xml")
  W <- writeUpdatedXML(input1)
  output <- "output.xml"
  withr::with_tempfile(output, {
    write(W, output)
    data2<-loadNULISAseq(output)
    expect_true(file.exists(output),output)
  })
})

test_that("CAL/IPC samples do not get AQ, aM, or dr attributes on ReadCount nodes", {
  skip_if_not(requireNamespace("NULISAseqAQ", quietly = TRUE),
              "NULISAseqAQ package not available")

  input1 <- test_path("fixtures", "XML_v1.3.0_with_AQ.xml")
  data <- loadNULISAseq(input1)
  W <- writeUpdatedXML(input1, data = data)
  doc <- xml2::read_xml(W)

  # Identify CAL/IPC sample barcodes
  cal_ipc_barcodes <- data$samples$sampleBarcode[data$samples$sampleType %in% c("Calibrator", "IPC")]
  sample_barcodes <- data$samples$sampleBarcode[data$samples$sampleType == "Sample"]

  # CAL/IPC ReadCount nodes should NOT have AQ, aM, or dr
  for (bc in cal_ipc_barcodes) {
    xpath <- sprintf(".//Data//Sample[@barcode='%s']/ReadCount", bc)
    readcounts <- xml2::xml_find_all(doc, xpath)
    if (length(readcounts) > 0) {
      expect_true(all(is.na(xml2::xml_attr(readcounts, "AQ"))),
                  info = paste("CAL/IPC barcode", bc, "should not have AQ"))
      expect_true(all(is.na(xml2::xml_attr(readcounts, "aM"))),
                  info = paste("CAL/IPC barcode", bc, "should not have aM"))
      expect_true(all(is.na(xml2::xml_attr(readcounts, "dr"))),
                  info = paste("CAL/IPC barcode", bc, "should not have dr"))
    }
  }

  # Regular sample ReadCount nodes SHOULD have AQ, aM, and dr
  has_aq <- FALSE
  for (bc in sample_barcodes[1:min(3, length(sample_barcodes))]) {
    xpath <- sprintf(".//Data//Sample[@barcode='%s']/ReadCount", bc)
    readcounts <- xml2::xml_find_all(doc, xpath)
    if (length(readcounts) > 0 && any(!is.na(xml2::xml_attr(readcounts, "AQ")))) {
      has_aq <- TRUE
      break
    }
  }
  expect_true(has_aq, info = "At least one regular sample should have AQ attribute")
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

  input1 <- test_path("fixtures", "Analysis_INF250_Lot4_AQ_LC_R3_20241229.xml")

  # Should warn but load successfully as RQ
  expect_warning(data <- loadNULISAseq(input1), "AQ metadata")

  # Verify it loaded successfully as RQ
  expect_false(is.null(data$Data))
  expect_true(is.null(data$AQ) || is.null(data$AQ$Data_AQ_aM))
})
