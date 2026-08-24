# Tests for XML-driven QC threshold overrides (issue #670)
# mergeQCCriteria overlays XML <QCThresholds> onto hardcoded QC*Criteria defaults.

test_that("NULL / empty xmlThresh returns defaults unchanged", {
  defaults <- NULISAseqR:::QCPlateCriteria(AQ = FALSE)
  expect_identical(NULISAseqR:::mergeQCCriteria(defaults, NULL), defaults)
  expect_identical(NULISAseqR:::mergeQCCriteria(defaults, list()), defaults)
})

test_that("NULL defaults (e.g. plain-RQ QCTargetCriteria) short-circuit without warning", {
  # QCTargetCriteria(advancedQC=FALSE, AQ=FALSE) returns NULL on a plain RQ
  # run; with defaults NULL, every xmlThresh key would otherwise be flagged
  # as "no matching computation" even though there's nothing to merge into.
  xmlThresh <- list(thresholds = c(Target_Min_Reads = "0.80"),
                     operators  = c(Target_Min_Reads = "<"))
  expect_no_warning(merged <- NULISAseqR:::mergeQCCriteria(NULL, xmlThresh))
  expect_null(merged)
})

test_that("forceDefaults ignores XML overrides", {
  defaults <- NULISAseqR:::QCPlateCriteria(AQ = FALSE)
  xmlThresh <- list(thresholds = c(ICRead_CV = "0.99"), operators = c(ICRead_CV = ">"))
  merged <- NULISAseqR:::mergeQCCriteria(defaults, xmlThresh, forceDefaults = TRUE)
  expect_identical(merged, defaults)
})

test_that("XML overrides existing plate thresholds but keeps flag set / thresholdNames", {
  defaults <- NULISAseqR:::QCPlateCriteria(AQ = FALSE)
  xmlThresh <- list(
    thresholds = c(ICRead_CV = "0.50", MinReads = "2e8"),
    operators  = c(ICRead_CV = ">"),
    format     = c(ICRead_CV = "percentage")
  )
  merged <- NULISAseqR:::mergeQCCriteria(defaults, xmlThresh)

  expect_equal(as.numeric(merged$thresholds[["ICRead_CV"]]), 0.50)
  expect_equal(as.numeric(merged$thresholds[["MinReads"]]), 2e8)
  # Untouched flag keeps default
  expect_equal(as.numeric(merged$thresholds[["IPCRead_CV"]]),
               as.numeric(defaults$thresholds[["IPCRead_CV"]]))
  # Flag identity/order and code-linked thresholdNames are preserved
  expect_identical(names(merged$thresholds), names(defaults$thresholds))
  expect_identical(merged$thresholdNames, defaults$thresholdNames)
})

test_that("XML flags with no matching computation are ignored, with a warning", {
  defaults <- NULISAseqR:::QCPlateCriteria(AQ = FALSE)
  xmlThresh <- list(thresholds = c(SampleCount = "96"), operators = c(SampleCount = "none"))
  expect_warning(
    merged <- NULISAseqR:::mergeQCCriteria(defaults, xmlThresh),
    "SampleCount"
  )
  expect_false("SampleCount" %in% names(merged$thresholds))
})

test_that("Detectability is permanently informational-only: XML cannot re-enable flagging (sample)", {
  # Real panel XML files commonly carry a Detectability threshold that mirrors
  # the *old* hardcoded default (e.g. plasma=0.9, operator="<"); if XML could
  # override Detectability, every such panel would silently re-enable
  # flagging the moment the hardcoded default became informational-only.
  defaults <- NULISAseqR:::QCSampleCriteria(TAP = TRUE)
  xmlThresh <- list(
    thresholds = c(Detectability_PLASMA = "0.900000"),
    operators  = c(Detectability_PLASMA = "<")
  )
  expect_warning(
    merged <- NULISAseqR:::mergeQCCriteria(defaults, xmlThresh, TAP = TRUE),
    "informational-only"
  )
  expect_identical(unname(merged$operators[["Detectability"]]), "none")
  expect_equal(as.numeric(merged$thresholds[["Detectability.plasma"]]), 0)
})

test_that("Detectability lockout also blocks adding a new per-matrix key from XML", {
  defaults <- NULISAseqR:::QCSampleCriteria(TAP = FALSE)  # lacks saliva
  expect_false("Detectability.saliva" %in% names(defaults$thresholds))
  xmlThresh <- list(
    thresholds = c(Detectability_SALIVA = "0.33"),
    operators  = c(Detectability_SALIVA = "<")
  )
  merged <- suppressWarnings(NULISAseqR:::mergeQCCriteria(defaults, xmlThresh, TAP = FALSE))
  expect_false("Detectability.saliva" %in% names(merged$thresholds))
})

test_that("Detectability is permanently informational-only: XML cannot re-enable flagging (plate)", {
  defaults <- NULISAseqR:::QCPlateCriteria(AQ = TRUE)
  xmlThresh <- list(thresholds = c(Detectability = "0.900000"), operators = c(Detectability = "<"))
  merged <- suppressWarnings(NULISAseqR:::mergeQCCriteria(defaults, xmlThresh))
  expect_identical(unname(merged$operators[["Detectability"]]), "none")
  expect_equal(as.numeric(merged$thresholds[["Detectability"]]), 0)
})

test_that("Detectability is permanently informational-only: XML cannot re-enable flagging (target)", {
  defaults <- NULISAseqR:::QCTargetCriteria(AQ = TRUE, advancedQC = TRUE)
  xmlThresh <- list(thresholds = c(Target_Detectability = "0.500000"), operators = c(Target_Detectability = "<"))
  merged <- suppressWarnings(NULISAseqR:::mergeQCCriteria(defaults, xmlThresh))
  expect_identical(unname(merged$operators[["Target_Detectability"]]), "none")
  expect_equal(as.numeric(merged$thresholds[["Target_Detectability"]]), 0)
})

test_that("appliesWhen='TAP' rows are dropped when TAP is FALSE", {
  defaults <- NULISAseqR:::QCSampleCriteria(TAP = TRUE)
  base_ic <- defaults$thresholds[["ICReads"]]
  xmlThresh <- list(
    thresholds  = c(ICReads = "7777"),
    operators   = c(ICReads = "<"),
    appliesWhen = c(ICReads = "TAP")
  )
  # TAP on -> applied
  merged_on <- NULISAseqR:::mergeQCCriteria(defaults, xmlThresh, TAP = TRUE)
  expect_equal(as.numeric(merged_on$thresholds[["ICReads"]]), 7777)
  # TAP off -> dropped, default retained
  merged_off <- NULISAseqR:::mergeQCCriteria(defaults, xmlThresh, TAP = FALSE)
  expect_equal(merged_off$thresholds[["ICReads"]], base_ic)
})

test_that("appliesWhen='TAP' gating applies to Plate criteria too (not sample-only)", {
  defaults <- NULISAseqR:::QCPlateCriteria(AQ = FALSE)
  base_min <- defaults$thresholds[["MinReads"]]
  xmlThresh <- list(
    thresholds  = c(MinReads = "12345"),
    operators   = c(MinReads = "<"),
    appliesWhen = c(MinReads = "TAP")
  )
  expect_equal(as.numeric(NULISAseqR:::mergeQCCriteria(defaults, xmlThresh, TAP = TRUE)$thresholds[["MinReads"]]), 12345)
  expect_equal(NULISAseqR:::mergeQCCriteria(defaults, xmlThresh, TAP = FALSE)$thresholds[["MinReads"]], base_min)
})

test_that("operator-only XML override (no value) is still applied", {
  defaults <- NULISAseqR:::QCPlateCriteria(AQ = FALSE)
  xmlThresh <- list(operators = c(ICRead_CV = "<"))  # no thresholds at all
  merged <- NULISAseqR:::mergeQCCriteria(defaults, xmlThresh)
  expect_equal(unname(merged$operators[["ICRead_CV"]]), "<")
})

test_that("single-value '<,>' XML threshold is expanded to symmetric -V,V form", {
  # XML expresses +/- ranges as a single positive value; evalCriterion needs
  # one threshold per operator. (Regression: Target_Conc_Accuracy value="0.3".)
  defaults <- NULISAseqR:::QCTargetCriteria(AQ = TRUE, advancedQC = FALSE)
  xmlThresh <- list(
    thresholds = c(Target_Conc_Accuracy = "0.300000"),
    operators  = c(Target_Conc_Accuracy = "<,>")
  )
  merged <- NULISAseqR:::mergeQCCriteria(defaults, xmlThresh)
  expect_equal(merged$thresholds[["Target_Conc_Accuracy"]], "-0.3,0.3")
  # evalCriterion must not error on the merged criterion
  expect_false(NULISAseqR:::evalCriterion(
    "Target_Conc_Accuracy", 0.1,
    merged$operators[["Target_Conc_Accuracy"]],
    merged$thresholds[["Target_Conc_Accuracy"]]))
  expect_true(NULISAseqR:::evalCriterion(
    "Target_Conc_Accuracy", 0.5,
    merged$operators[["Target_Conc_Accuracy"]],
    merged$thresholds[["Target_Conc_Accuracy"]]))
})

test_that("two-value '<,>' threshold is left unchanged", {
  defaults <- NULISAseqR:::QCSampleCriteria(TAP = TRUE)
  xmlThresh <- list(
    thresholds = c(IC_Median = "-0.25,0.25"),
    operators  = c(IC_Median = "<,>")
  )
  merged <- NULISAseqR:::mergeQCCriteria(defaults, xmlThresh, TAP = TRUE)
  expect_equal(merged$thresholds[["IC_Median"]], "-0.25,0.25")
})

test_that("evalCriterion treats operator='none' as never-flagging", {
  expect_false(NULISAseqR:::evalCriterion("Info", 5, "none", NA))
  expect_false(NULISAseqR:::evalCriterion("Info", NA, "none", NA))
})

test_that("evalCriterion still flags normal operators", {
  expect_true(NULISAseqR:::evalCriterion("x", 10, ">", 5))
  expect_false(NULISAseqR:::evalCriterion("x", 3, ">", 5))
  expect_true(NULISAseqR:::evalCriterion("x", NA, "<", 5))  # NA value flags
})

# End-to-end: QCFlagPlate/QCFlagSample actually consume mergeQCCriteria's
# output when producing flags, not just when comparing threshold values in
# isolation. These confirm an XML override changes the flagged status of a
# real QC table, and forceDefaults restores the hardcoded behavior.

test_that("QCFlagPlate applies an XML-overridden MinReads threshold to its output", {
  targets <- data.frame(targetName = c("IC1", "T1", "T2"),
                         targetType = c("control", "target", "target"),
                         stringsAsFactors = FALSE)
  samples <- data.frame(sampleName = c("S1", "S2", "S3"),
                         sampleType = c("IPC", "IPC", "Sample"),
                         stringsAsFactors = FALSE)
  raw <- matrix(c(100, 200, 50, 120, 180, 60, 90, 210, 55), nrow = 3, byrow = TRUE,
                dimnames = list(c("IC1", "T1", "T2"), c("S1", "S2", "S3")))
  normed <- raw / rowMeans(raw)
  aboveLOD <- matrix(TRUE, nrow = 2, ncol = 1, dimnames = list(c("T1", "T2"), "S3"))
  nReads <- sum(raw, na.rm = TRUE)

  base <- NULISAseqR:::QCFlagPlate(raw, normed, aboveLOD, targets, samples, AQ = FALSE, TAP = FALSE)
  expect_true(as.logical(base$status[base$flagName == "MinReads"]))  # nReads << hardcoded default (1e8) -> flags

  xmlThresh <- list(thresholds = c(MinReads = as.character(nReads - 1)),
                     operators  = c(MinReads = "<"))
  overridden <- NULISAseqR:::QCFlagPlate(raw, normed, aboveLOD, targets, samples, AQ = FALSE, TAP = FALSE,
                                          xmlThresh = xmlThresh)
  expect_false(as.logical(overridden$status[overridden$flagName == "MinReads"]))  # nReads not < nReads-1 -> no flag
  expect_equal(as.numeric(overridden$QCthreshold[overridden$flagName == "MinReads"]), nReads - 1)

  forced <- NULISAseqR:::QCFlagPlate(raw, normed, aboveLOD, targets, samples, AQ = FALSE, TAP = FALSE,
                                      xmlThresh = xmlThresh, forceDefaults = TRUE)
  expect_identical(forced, base)
})

test_that("QCFlagSample applies an XML-overridden NumReads threshold to its output", {
  targets <- data.frame(targetName = c("IC1", "T1", "T2"),
                         targetType = c("control", "target", "target"),
                         noDetectability = c(FALSE, FALSE, FALSE),
                         stringsAsFactors = FALSE)
  samples <- data.frame(sampleName = c("S1", "S2", "S3"),
                         sampleType = c("IPC", "IPC", "Sample"),
                         sampleBarcode = c("B1", "B2", "B3"),
                         SAMPLE_MATRIX = c("plasma", "plasma", "plasma"),
                         stringsAsFactors = FALSE)
  raw <- matrix(c(1000, 1000, 1000, 20, 30, 25, 40, 35, 45), nrow = 3, byrow = TRUE,
                dimnames = list(c("IC1", "T1", "T2"), c("S1", "S2", "S3")))
  aboveLOD <- matrix(TRUE, nrow = 2, ncol = 3, dimnames = list(c("T1", "T2"), c("S1", "S2", "S3")))

  base <- NULISAseqR:::QCFlagSample(raw, aboveLOD, samples, targets, TAP = FALSE)
  base_numreads <- base[base$flagName == "NumReads", ]
  expect_true(all(as.logical(base_numreads$status)))

  xmlThresh <- list(thresholds = c(NumReads = "1000"), operators = c(NumReads = "<"))
  overridden <- NULISAseqR:::QCFlagSample(raw, aboveLOD, samples, targets, TAP = FALSE, xmlThresh = xmlThresh)
  overridden_numreads <- overridden[overridden$flagName == "NumReads", ]
  expect_true(all(!as.logical(overridden_numreads$status)))
  expect_true(all(overridden_numreads$QCthreshold == "1000"))

  forced <- NULISAseqR:::QCFlagSample(raw, aboveLOD, samples, targets, TAP = FALSE,
                                       xmlThresh = xmlThresh, forceDefaults = TRUE)
  expect_identical(forced, base)
})
