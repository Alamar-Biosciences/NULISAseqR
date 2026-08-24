# Tests for bug fixes in insertCovariatesXML.

# clean_covariate_names() emits logger::log_info() lines that suppressMessages()
# does not catch (logger uses its own appender, not base R's message condition).
# Lift the threshold to WARN for the test scope and restore it on exit.
local_silence_logger <- function(.local_envir = parent.frame()) {
  old <- logger::log_threshold()
  withr::defer(logger::log_threshold(old), envir = .local_envir)
  logger::log_threshold(logger::WARN)
  invisible()
}

test_that("insertCovariatesXML writes covariate when sampleName is plate-suffixed (issue #3199)", {
  # Regression: when NAS exports an updated XML, it passes a covariates df
  # whose sampleName has a plate-filename suffix to disambiguate cross-plate
  # duplicates (e.g. QC samples that appear on every plate). The XML Barcode
  # element's text is the bare sample name — a strict prefix of sampleName.
  # The function previously silently skipped these rows, leaving the attribute
  # unwritten and re-imports producing NA for QC samples.
  xml_string <- '<?xml version="1.0"?>
  <NULISAseq>
    <BarcodeB>
      <Barcode AUTO_SAMPLENAME="LowQC_Plasma_begin" AUTO_WELLPOSITION="A01">LowQC_Plasma_begin</Barcode>
      <Barcode AUTO_SAMPLENAME="Patient_001" AUTO_WELLPOSITION="B01">Patient_001</Barcode>
    </BarcodeB>
  </NULISAseq>'

  tmp_xml <- withr::local_tempfile(fileext = ".xml")
  writeLines(xml_string, tmp_xml)

  covariates <- data.frame(
    sampleName = c("LowQC_Plasma_begin_20260330-1430_Bay1_Panel.xml",
                   "Patient_001"),
    diagnosis  = c("exemplarydiagnosis", "exemplarydiagnosis"),
    stringsAsFactors = FALSE
  )

  result <- insertCovariatesXML(tmp_xml, covariates)

  barcode_nodes <- xml2::xml_find_all(result, "//BarcodeB/Barcode")
  expect_equal(xml2::xml_attr(barcode_nodes[[1]], "diagnosis"), "exemplarydiagnosis")
  expect_equal(xml2::xml_attr(barcode_nodes[[2]], "diagnosis"), "exemplarydiagnosis")
})

test_that("insertCovariatesXML writes nothing when sampleName is absent", {
  xml_string <- '<?xml version="1.0"?>
  <NULISAseq>
    <BarcodeB>
      <Barcode AUTO_SAMPLENAME="Sample1" AUTO_WELLPOSITION="A01">Sample1</Barcode>
    </BarcodeB>
  </NULISAseq>'

  tmp_xml <- withr::local_tempfile(fileext = ".xml")
  writeLines(xml_string, tmp_xml)

  covariates <- data.frame(
    sampleName  = "NonExistent",
    MyCovariate = "ValueA",
    stringsAsFactors = FALSE
  )

  result <- insertCovariatesXML(tmp_xml, covariates)
  barcode_nodes <- xml2::xml_find_all(result, "//BarcodeB/Barcode")
  expect_true(is.na(xml2::xml_attr(barcode_nodes[[1]], "MyCovariate")))
})

test_that("insertCovariatesXML preserves exact-match happy path", {
  # The plate-suffix fix must not regress callers that pass an exact-match
  # sampleName (the historical contract before the suffix workflow).
  xml_string <- '<?xml version="1.0"?>
  <NULISAseq>
    <BarcodeB>
      <Barcode AUTO_SAMPLENAME="S1" AUTO_WELLPOSITION="A01">S1</Barcode>
    </BarcodeB>
  </NULISAseq>'

  tmp_xml <- withr::local_tempfile(fileext = ".xml")
  writeLines(xml_string, tmp_xml)

  covariates <- data.frame(
    sampleName = "S1",
    cohort     = "A",
    stringsAsFactors = FALSE
  )

  result <- insertCovariatesXML(tmp_xml, covariates)
  expect_equal(xml2::xml_attr(xml2::xml_find_first(result, "//Barcode"), "cohort"), "A")
})

test_that("insertCovariatesXML does not over-match a different sample with a prefix-shared name", {
  # Pathological substring-match: element text "S1" is a substring of a single
  # covariates row "S100". An unanchored grep would treat this as a unique hit
  # and write S100's covariates onto S1's Barcode. The prefix-with-"_" match
  # rejects this — there is no row equal to "S1" nor starting with "S1_".
  xml_string <- '<?xml version="1.0"?>
  <NULISAseq>
    <BarcodeB>
      <Barcode AUTO_SAMPLENAME="S1" AUTO_WELLPOSITION="A01">S1</Barcode>
    </BarcodeB>
  </NULISAseq>'

  tmp_xml <- withr::local_tempfile(fileext = ".xml")
  writeLines(xml_string, tmp_xml)

  covariates <- data.frame(
    sampleName = "S100",
    cohort     = "WRONG",
    stringsAsFactors = FALSE
  )

  result <- insertCovariatesXML(tmp_xml, covariates)
  expect_true(is.na(xml2::xml_attr(xml2::xml_find_first(result, "//Barcode"), "cohort")))
})

test_that("insertCovariatesXML skips Barcode elements with empty or whitespace-only text", {
  # Empty / blank barcode text would make the prefix check match every row
  # whose sampleName starts with "_". Guard so no attributes are written.
  xml_string <- '<?xml version="1.0"?>
  <NULISAseq>
    <BarcodeB>
      <Barcode AUTO_SAMPLENAME="" AUTO_WELLPOSITION="A01"></Barcode>
      <Barcode AUTO_SAMPLENAME="" AUTO_WELLPOSITION="B01">   </Barcode>
    </BarcodeB>
  </NULISAseq>'

  tmp_xml <- withr::local_tempfile(fileext = ".xml")
  writeLines(xml_string, tmp_xml)

  covariates <- data.frame(
    sampleName = "_prefix_match",
    cohort     = "WRONG",
    stringsAsFactors = FALSE
  )

  result <- insertCovariatesXML(tmp_xml, covariates)
  barcode_nodes <- xml2::xml_find_all(result, "//BarcodeB/Barcode")
  expect_true(all(is.na(xml2::xml_attr(barcode_nodes, "cohort"))))
})

test_that("insertCovariatesXML preserves verbatim attribute names by default", {
  # The historical contract is write-verbatim. Verifies sanitize_names defaults
  # to FALSE and an unsanitized column name reaches the XML unchanged.
  xml_string <- '<?xml version="1.0"?>
  <NULISAseq>
    <BarcodeB>
      <Barcode AUTO_WELLPOSITION="A01">S1</Barcode>
    </BarcodeB>
  </NULISAseq>'

  tmp_xml <- withr::local_tempfile(fileext = ".xml")
  writeLines(xml_string, tmp_xml)

  covariates <- data.frame(
    sampleName = "S1",
    Disease    = "AD",
    stringsAsFactors = FALSE
  )

  result <- insertCovariatesXML(tmp_xml, covariates)
  node <- xml2::xml_find_first(result, "//Barcode")
  expect_equal(xml2::xml_attr(node, "Disease"), "AD")
  expect_true(is.na(xml2::xml_attr(node, "disease")))
})

test_that("insertCovariatesXML sanitizes attribute names when sanitize_names = TRUE", {
  # When opted in, attribute names are sanitized via clean_covariate_names
  # (case = "lower_camel") so they match NAS's covariate-upload sanitation.
  # This keeps the attribute name stable across NAS export -> NULISAseqR
  # import -> NAS re-import (issue #3199 follow-up).
  xml_string <- '<?xml version="1.0"?>
  <NULISAseq>
    <BarcodeB>
      <Barcode AUTO_WELLPOSITION="A01">S1</Barcode>
    </BarcodeB>
  </NULISAseq>'

  tmp_xml <- withr::local_tempfile(fileext = ".xml")
  writeLines(xml_string, tmp_xml)

  # The multi-word column requires check.names = FALSE and exercises the
  # lower_camel transform (whitespace -> camelCase), not just single-word
  # lowercasing.
  covariates <- data.frame(
    sampleName       = "S1",
    Disease          = "AD",
    `Sample Cohort`  = "PD",
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )

  local_silence_logger()
  result <- insertCovariatesXML(tmp_xml, covariates, sanitize_names = TRUE)
  node <- xml2::xml_find_first(result, "//Barcode")
  expect_equal(xml2::xml_attr(node, "disease"), "AD")
  expect_true(is.na(xml2::xml_attr(node, "Disease")))
  expect_equal(xml2::xml_attr(node, "sampleCohort"), "PD")
  expect_true(is.na(xml2::xml_attr(node, "Sample Cohort")))
})

test_that("insertCovariatesXML handles multi-match candidates that share no well position", {
  # multi_check > 1 routes into the well-position disambiguation branch. If
  # the well filter yields zero matching rows, the old `!is.null && !is.na &&
  # length > 0` guard relied on short-circuit order to avoid an R >= 4.3
  # length-mismatch error in `&&`. Verify the function returns cleanly and
  # writes no covariate attribute on the unmatched element.
  xml_string <- '<?xml version="1.0"?>
  <NULISAseq>
    <BarcodeB>
      <Barcode AUTO_WELLPOSITION="A01" AUTO_WELLROW="A" AUTO_WELLCOL="1">LowQC</Barcode>
    </BarcodeB>
  </NULISAseq>'

  tmp_xml <- withr::local_tempfile(fileext = ".xml")
  writeLines(xml_string, tmp_xml)

  # Two prefix-matching rows, neither at well A01 — well filter -> 0 rows.
  covariates <- data.frame(
    sampleName = c("LowQC_plateA.xml", "LowQC_plateB.xml"),
    wellRow    = c("Z", "Z"),
    wellCol    = c(99, 99),
    diagnosis  = c("a", "b"),
    stringsAsFactors = FALSE
  )

  result <- insertCovariatesXML(tmp_xml, covariates)
  node <- xml2::xml_find_first(result, "//Barcode")
  expect_true(is.na(xml2::xml_attr(node, "diagnosis")))
})

test_that("insertCovariatesXML errors when sanitize_names produces colliding attribute names", {
  # `sampleName` is protected (kept verbatim) and `"Sample Name"` sanitizes
  # to `sampleName` — both would target the same XML attribute. Surface as
  # an error rather than silently letting one overwrite the other.
  xml_string <- '<?xml version="1.0"?>
  <NULISAseq>
    <BarcodeB>
      <Barcode AUTO_WELLPOSITION="A01">S1</Barcode>
    </BarcodeB>
  </NULISAseq>'

  tmp_xml <- withr::local_tempfile(fileext = ".xml")
  writeLines(xml_string, tmp_xml)

  covariates <- data.frame(
    sampleName       = "S1",
    `Sample Name`    = "duplicate",
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )

  local_silence_logger()
  expect_error(
    insertCovariatesXML(tmp_xml, covariates, sanitize_names = TRUE),
    "colliding XML attribute names"
  )
})

test_that("insertCovariatesXML sanitization leaves sampleName (protected) untouched", {
  # clean_covariate_names protects identifier columns like sampleName. Confirm
  # the sanitation map does not mangle them when sanitize_names = TRUE.
  xml_string <- '<?xml version="1.0"?>
  <NULISAseq>
    <BarcodeB>
      <Barcode AUTO_WELLPOSITION="A01">S1</Barcode>
    </BarcodeB>
  </NULISAseq>'

  tmp_xml <- withr::local_tempfile(fileext = ".xml")
  writeLines(xml_string, tmp_xml)

  covariates <- data.frame(
    sampleName = "S1",
    cohort     = "A",
    stringsAsFactors = FALSE
  )

  local_silence_logger()
  result <- insertCovariatesXML(tmp_xml, covariates, sanitize_names = TRUE)
  node <- xml2::xml_find_first(result, "//Barcode")
  expect_equal(xml2::xml_attr(node, "sampleName"), "S1")
  expect_equal(xml2::xml_attr(node, "cohort"), "A")
})
