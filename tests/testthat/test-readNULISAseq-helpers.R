# Tests for helper functions in readNULISAseq.R

# =============================================================================
# Tests for calcSampleTargetNAs()
# =============================================================================

test_that("calcSampleTargetNAs returns empty matrix when no matches", {
  curve_quant <- c("F", "F", "F")
  sample_matrices <- c("PLASMA", "SERUM", "CSF")

  result <- calcSampleTargetNAs(curve_quant, sample_matrices)

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 2)
  expect_equal(colnames(result), c("row", "col"))
  expect_equal(nrow(result), 0)
})

test_that("calcSampleTargetNAs identifies explicit NA patterns", {
  curve_quant <- c("F;CSF-NA", "F", "F;PLASMA-NA")
  sample_matrices <- c("CSF", "PLASMA", "SERUM")

  result <- calcSampleTargetNAs(curve_quant, sample_matrices)

  expect_true(is.matrix(result))
  # Should find CSF-NA in row 1 (col 1) and PLASMA-NA in row 3 (col 2)
  expect_equal(nrow(result), 2)
  expect_true(any(result[, "row"] == 1 & result[, "col"] == 1))  # CSF-NA
  expect_true(any(result[, "row"] == 3 & result[, "col"] == 2))  # PLASMA-NA
})

test_that("calcSampleTargetNAs handles reverse curve targets correctly", {
  # Reverse curve targets (R) should be NA for non-PLASMA/SERUM/CONTROL matrices
  curve_quant <- c("R", "F", "R")
  sample_matrices <- c("PLASMA", "CSF", "SERUM")

  result <- calcSampleTargetNAs(curve_quant, sample_matrices)

  # R targets should be marked NA for CSF (col 2) but not PLASMA or SERUM
  # Row 1 (R) + col 2 (CSF) = match
  # Row 3 (R) + col 2 (CSF) = match
  expect_true(any(result[, "row"] == 1 & result[, "col"] == 2))
  expect_true(any(result[, "row"] == 3 & result[, "col"] == 2))

  # PLASMA and SERUM should NOT be marked as NA for R targets
  expect_false(any(result[, "row"] == 1 & result[, "col"] == 1))  # R + PLASMA
  expect_false(any(result[, "row"] == 1 & result[, "col"] == 3))  # R + SERUM
})

test_that("calcSampleTargetNAs handles CONTROL matrix with reverse curves", {
  curve_quant <- c("R", "R")
  sample_matrices <- c("CONTROL", "URINE")

  result <- calcSampleTargetNAs(curve_quant, sample_matrices)

  # CONTROL should NOT be marked as NA for R targets
  # URINE should be marked as NA for R targets
  expect_false(any(result[, "col"] == 1))  # No CONTROL matches
  expect_true(any(result[, "row"] == 1 & result[, "col"] == 2))   # R + URINE
  expect_true(any(result[, "row"] == 2 & result[, "col"] == 2))   # R + URINE
})

test_that("calcSampleTargetNAs is case insensitive for NA patterns", {
  curve_quant <- c("r;csf-na", "F")
  sample_matrices <- c("CSF", "PLASMA")

  result <- calcSampleTargetNAs(curve_quant, sample_matrices)

  # Should find csf-na match (case insensitive - curve_quant is lowercased)
  expect_true(any(result[, "row"] == 1 & result[, "col"] == 1))
  # R target should NOT mark PLASMA as NA (PLASMA is in the allowed list)
  expect_false(any(result[, "row"] == 1 & result[, "col"] == 2))
})

test_that("calcSampleTargetNAs exception list is case sensitive", {
  # Note: The exception list (PLASMA, SERUM, CONTROL) is case-sensitive
  # lowercase sample_matrices will NOT match the exception list
  curve_quant <- c("R", "R")
  sample_matrices <- c("plasma", "PLASMA")

  result <- calcSampleTargetNAs(curve_quant, sample_matrices)

  # lowercase "plasma" should be marked as NA (doesn't match exception list)
  expect_true(any(result[, "row"] == 1 & result[, "col"] == 1))
  expect_true(any(result[, "row"] == 2 & result[, "col"] == 1))
  # uppercase "PLASMA" should NOT be marked as NA (matches exception list)
  expect_false(any(result[, "col"] == 2))
})

test_that("calcSampleTargetNAs handles empty inputs", {
  # Empty curve_quant
  result1 <- calcSampleTargetNAs(character(0), c("PLASMA", "CSF"))
  expect_equal(nrow(result1), 0)


  # Empty sample_matrices
  result2 <- calcSampleTargetNAs(c("R", "F"), character(0))
  expect_equal(nrow(result2), 0)
})

# =============================================================================
# Tests for readQCXMLNode()
# =============================================================================

test_that("readQCXMLNode returns NULL for NULL or empty input", {
  expect_null(readQCXMLNode(NULL))
  expect_null(readQCXMLNode(list()))
})

test_that("readQCXMLNode parses plate-level QCFlag nodes", {
  xml_text <- '
  <PlateQC>
    <QCFlag name="TotalReads" set="pass" explanation="Reads above threshold">1000000</QCFlag>
    <QCFlag name="MappedReads" set="warn" explanation="Mapped reads low">500000</QCFlag>
  </PlateQC>'

  xml_doc <- xml2::read_xml(xml_text)
  nodes <- xml2::xml_find_all(xml_doc, ".//QCFlag")

  result <- readQCXMLNode(nodes, rename = TRUE)

  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2)
  expect_true("flagName" %in% colnames(result))  # renamed from "name"
  expect_true("status" %in% colnames(result))    # renamed from "set"
  expect_true("val" %in% colnames(result))       # renamed from "value"
  expect_true("explanations" %in% colnames(result))  # renamed from "explanation"

  expect_true("TotalReads" %in% result$flagName)
  expect_true("MappedReads" %in% result$flagName)
})

test_that("readQCXMLNode parses sample-level nested QCFlag nodes", {
  xml_text <- '
  <SampleQC>
    <Sample name="Sample1">
      <QCFlag name="IC_Median" set="pass">100</QCFlag>
      <QCFlag name="SampleReads" set="pass">50000</QCFlag>
    </Sample>
    <Sample name="Sample2">
      <QCFlag name="IC_Median" set="fail">10</QCFlag>
    </Sample>
  </SampleQC>'

  xml_doc <- xml2::read_xml(xml_text)
  nodes <- xml2::xml_find_all(xml_doc, ".//Sample")

  result <- readQCXMLNode(nodes, rename = TRUE)

  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 3)  # 2 flags for Sample1 + 1 for Sample2
  expect_true("sample" %in% colnames(result))
  expect_true("Sample1" %in% result$sample)
  expect_true("Sample2" %in% result$sample)
})

test_that("readQCXMLNode parses target-level nodes with custom tag", {
  xml_text <- '
  <TargetQC>
    <Target name="TARGET1">
      <SC_conc aq_aM="500.0" aq_pgmL="250.0"/>
    </Target>
    <Target name="TARGET2">
      <SC_conc aq_aM="1000.0" aq_pgmL="500.0"/>
    </Target>
  </TargetQC>'

  xml_doc <- xml2::read_xml(xml_text)
  nodes <- xml2::xml_find_all(xml_doc, ".//Target")

  result <- readQCXMLNode(nodes, rename = TRUE, tag = "SC_conc")

  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2)
  expect_true("target" %in% colnames(result))
  expect_true("SC_conc_aM" %in% colnames(result))      # renamed from aq_aM
  expect_true("SC_conc_pg_ml" %in% colnames(result))   # renamed from aq_pgmL

  expect_true("TARGET1" %in% result$target)
  expect_true("TARGET2" %in% result$target)
})

test_that("readQCXMLNode respects rename=FALSE", {
  xml_text <- '
  <PlateQC>
    <QCFlag name="TotalReads" set="pass">1000000</QCFlag>
  </PlateQC>'

  xml_doc <- xml2::read_xml(xml_text)
  nodes <- xml2::xml_find_all(xml_doc, ".//QCFlag")

  result <- readQCXMLNode(nodes, rename = FALSE)

  expect_true("name" %in% colnames(result))   # NOT renamed
  expect_true("set" %in% colnames(result))    # NOT renamed
  expect_true("value" %in% colnames(result))  # NOT renamed
})

# =============================================================================
# Tests for readQCThresholdXMLNode()
# =============================================================================

test_that("readQCThresholdXMLNode returns NULL for NULL or empty input", {
  expect_null(readQCThresholdXMLNode(NULL))
  expect_null(readQCThresholdXMLNode(list()))
})

test_that("readQCThresholdXMLNode parses threshold attributes", {
  xml_text <- '
  <Thresholds>
    <Threshold name="IC_Median" operator="gte" properName="IC Median" explanation="Internal control median">50</Threshold>
    <Threshold name="TotalReads" operator="gte" properName="Total Reads" explanation="Total read count">100000</Threshold>
  </Thresholds>'

  xml_doc <- xml2::read_xml(xml_text)
  nodes <- xml2::xml_find_all(xml_doc, ".//Threshold")

  result <- readQCThresholdXMLNode(nodes, rename = TRUE)

  expect_true(is.list(result))
  expect_true("thresholds" %in% names(result))    # renamed from "value"
  expect_true("operators" %in% names(result))     # renamed from "operator"
  expect_true("properNames" %in% names(result))   # renamed from "properName"
  expect_true("explanations" %in% names(result))  # renamed from "explanation"

  # Check named vectors
  expect_equal(result$thresholds["IC_Median"], c(IC_Median = "50"))
  expect_equal(result$thresholds["TotalReads"], c(TotalReads = "100000"))
  expect_equal(result$operators["IC_Median"], c(IC_Median = "gte"))
})

test_that("readQCThresholdXMLNode uses text value over attribute", {
  xml_text <- '
  <Thresholds>
    <Threshold name="Test" value="wrong">correct</Threshold>
  </Thresholds>'

  xml_doc <- xml2::read_xml(xml_text)
  nodes <- xml2::xml_find_all(xml_doc, ".//Threshold")

  result <- readQCThresholdXMLNode(nodes, rename = TRUE)

  # Text content should override value attribute
  expect_equal(result$thresholds["Test"], c(Test = "correct"))
})

test_that("readQCThresholdXMLNode respects rename=FALSE", {
  xml_text <- '
  <Thresholds>
    <Threshold name="Test" operator="gte">100</Threshold>
  </Thresholds>'

  xml_doc <- xml2::read_xml(xml_text)
  nodes <- xml2::xml_find_all(xml_doc, ".//Threshold")

  result <- readQCThresholdXMLNode(nodes, rename = FALSE)

  expect_true("value" %in% names(result))     # NOT renamed
  expect_true("operator" %in% names(result))  # NOT renamed
})
