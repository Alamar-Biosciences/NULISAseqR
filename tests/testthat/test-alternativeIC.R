# Define a test processXML

# These fixtures have no Curve_Quant values, so readNULISAseq() always emits
# "No Curve_Quant values found in XML..." alongside whatever else it warns
# about. Muffle just that one so expect_warning() below can match on the
# IC-removal warning without depending on warning order/count.
muffle_curve_quant_warning <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("Curve_Quant", conditionMessage(w))) invokeRestart("muffleWarning")
  })
}

test_that("Test that loadNULISAseq reads in an XML file and can use an alternate IC", {

  # Test older XMLs that do not have target types, but manually allow setting the IC
  input1 <- test_path("fixtures", "detectability_P1_Tr03.xml")
  data1 <- suppressWarnings(loadNULISAseq(input1, IPC=NULL, IC='mCherry', SC=NULL, allowMissingCurveQuant=TRUE))
  data2 <- suppressWarnings(loadNULISAseq(input1, IPC=NULL, IC='CCL7', SC=NULL, allowMissingCurveQuant=TRUE))
  expect_true(data1$IC_normed$normData[1,1] != data2$IC_normed$normData[1,1])
  expect_true(data1$IC == "mCherry")
  expect_true(data2$IC == "CCL7")

  # Test a newer XML that has targetType but no modifiers
  input3 <- test_path("fixtures", "detectability_P1_Tr03_typemCherry.xml")
  data3 <- suppressWarnings(loadNULISAseq(input3, IPC=NULL, IC=NULL, SC=NULL, allowMissingCurveQuant=TRUE))
  expect_true(data3$IC == "mCherry")
  expect_true(all(data3$IC_normed$normData == data1$IC_normed$normData))

  input4 <- test_path("fixtures", "detectability_P1_Tr03_typeCCL7.xml")
  data4 <- suppressWarnings(loadNULISAseq(input4, IPC=NULL, IC=NULL, SC=NULL, allowMissingCurveQuant=TRUE))
  expect_true(data4$IC == "CCL7")
  expect_true(all(data4$IC_normed$normData == data2$IC_normed$normData))

  # Test a newer XML that has targetType and has multiple ICs (geometric mean)
  input5 <- test_path("fixtures", "detectability_P1_Tr03_typemCherry_CCL7.xml")
  data5 <- suppressWarnings(loadNULISAseq(input5, IPC=NULL, IC=NULL, SC=NULL, allowMissingCurveQuant=TRUE))
  expect_true("mCherry" %in% data5$IC)
  expect_true("CCL7" %in% data5$IC)
  expect_false(all(data5$IC_normed$normData == data4$IC_normed$normData))
  expect_false(all(data5$IC_normed$normData == data3$IC_normed$normData))

  # Test that when multiple ICs exist, that only one can be selected
  data5 <- suppressWarnings(loadNULISAseq(input5, IPC=NULL, IC="mCherry", SC=NULL, allowMissingCurveQuant=TRUE))
  expect_true("mCherry" %in% data5$IC) 
  expect_false("CCL7" %in% data5$IC)
  expect_true(all(data5$IC_normed$normData == data3$IC_normed$normData))
  
  data5 <- suppressWarnings(loadNULISAseq(input5, IPC=NULL, IC="CCL7", SC=NULL, allowMissingCurveQuant=TRUE))
  expect_false("mCherry" %in% data5$IC) 
  expect_true("CCL7" %in% data5$IC)
  expect_true(all(data5$IC_normed$normData == data4$IC_normed$normData))

  # Test a newer XML that has targetType and uses hide on a control
  # Hidden targets are completely removed from data
  input6 <- test_path("fixtures", "detectability_P1_Tr03_typemCherry_CCL7_hide.xml")
  data6 <- suppressWarnings(loadNULISAseq(input6, IPC=NULL, IC=NULL, SC=NULL, allowMissingCurveQuant=TRUE))
  expect_true("CCL7" %in% data6$IC)
  expect_false("mCherry" %in% data6$IC)
  expect_false("mCherry" %in% data6$targets$targetName, info = "Hidden target mCherry should be removed from targets")

  # Test a newer XML that has targetType and uses noDetectability
  input7 <- test_path("fixtures", "detectability_P1_Tr03_typemCherry_noDetectability_AGER.xml")
  data7 <- suppressWarnings(loadNULISAseq(input7, IPC=NULL, IC=NULL, SC=NULL, allowMissingCurveQuant=TRUE))
  expect_false("AGER" %in% data7$detectability$all$detectable)
  expect_false("AGER" %in% data7$detectability$all$detectablity)

  # Test a newer XML that has targetType and uses hide on a control and a target
  # Hidden targets are completely removed from data
  input6 <- test_path("fixtures", "detectability_P1_Tr03_typemCherry_CCL7_hide2.xml")
  data6 <- suppressWarnings(loadNULISAseq(input6, IPC=NULL, IC=NULL, SC=NULL, allowMissingCurveQuant=TRUE))
  expect_true("CCL7" %in% data6$IC)
  expect_false("mCherry" %in% data6$IC)
  expect_false("mCherry" %in% data6$targets$targetName, info = "Hidden target mCherry should be removed from targets")
  expect_false("WNT7A" %in% data6$targets$targetName, info = "Hidden target WNT7A should be removed from targets")
  expect_false("mCherry" %in% rownames(data6$Data), info = "Hidden target mCherry should be removed from Data matrix")
  expect_false("WNT7A" %in% rownames(data6$Data), info = "Hidden target WNT7A should be removed from Data matrix")
})

test_that("excludeTargets removing the only explicitly-passed IC errors clearly, naming the target, instead of a bare subscript-out-of-bounds", {
  # Regression: excludeTargets dropped the IC target's row from Data without
  # updating the (explicitly-passed) IC value, so intraPlateNorm() crashed on
  # data_matrix[IC,] with a bare "subscript out of bounds" instead of a
  # message a caller can act on.
  input1 <- test_path("fixtures", "detectability_P1_Tr03.xml")
  expect_error(
    suppressWarnings(loadNULISAseq(input1, IPC=NULL, IC='mCherry', SC=NULL,
                                    excludeTargets='mCherry', allowMissingCurveQuant=TRUE)),
    "IC target.*not present.*mCherry"
  )

  # excludeTargets on an unrelated target still normalizes fine using IC.
  data8 <- suppressWarnings(loadNULISAseq(input1, IPC=NULL, IC='mCherry', SC=NULL,
                                           excludeTargets='CCL7', allowMissingCurveQuant=TRUE))
  expect_identical(data8$IC, "mCherry")
  expect_false("CCL7" %in% rownames(data8$Data))
  expect_true("mCherry" %in% rownames(data8$Data))
})

test_that("excludeTargets removing one of several explicitly-passed ICs warns and normalizes on the survivor, rather than silently changing the normalization basis", {
  input5 <- test_path("fixtures", "detectability_P1_Tr03_typemCherry_CCL7.xml")
  expect_warning(
    data9 <- muffle_curve_quant_warning(
      loadNULISAseq(input5, IPC=NULL, IC=c('mCherry', 'CCL7'), SC=NULL,
                    excludeTargets='CCL7', allowMissingCurveQuant=TRUE)
    ),
    "CCL7"
  )
  expect_identical(data9$IC, "mCherry")
  expect_false("CCL7" %in% rownames(data9$Data))

  # Must match single-IC normalization on the same (CCL7-excluded) data, not the
  # geometric-mean result computed when both ICs are present.
  dataSingle <- suppressWarnings(loadNULISAseq(input5, IPC=NULL, IC='mCherry', SC=NULL,
                                                excludeTargets='CCL7', allowMissingCurveQuant=TRUE))
  expect_equal(data9$IC_normed$normData, dataSingle$IC_normed$normData)
})

test_that("hide=TRUE removing an explicitly-passed IC (not just an auto-detected one) errors clearly", {
  # The XML-driven auto-detect path already filters hidden controls out of
  # raw$IC inside readNULISAseq() itself (see the earlier hide tests above);
  # this exercises the gap where a caller passes IC= explicitly and the
  # matching target happens to carry hide=TRUE.
  input6 <- test_path("fixtures", "detectability_P1_Tr03_typemCherry_CCL7_hide.xml")
  expect_error(
    suppressWarnings(loadNULISAseq(input6, IPC=NULL, IC='mCherry', SC=NULL,
                                    allowMissingCurveQuant=TRUE)),
    "IC target.*not present.*mCherry"
  )

  # Passing both ICs explicitly recovers by normalizing on the one that survived hiding.
  expect_warning(
    data10 <- muffle_curve_quant_warning(
      loadNULISAseq(input6, IPC=NULL, IC=c('mCherry', 'CCL7'), SC=NULL,
                    allowMissingCurveQuant=TRUE)
    ),
    "mCherry"
  )
  expect_identical(data10$IC, "CCL7")
})
