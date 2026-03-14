test_that("Reverse curve targets automatically get noDetectability flag", {
  input <- test_path("fixtures", "Analysis_INF250_Lot4_AQ_LC_R3_20241229.xml")
  data <- loadNULISAseq(input)

  reverse_targets <- NULISAseqR:::get_reverse_curve_targets(data$targets)
  expect_true(length(reverse_targets) > 0, info = "Fixture should have reverse curve targets")

  # All reverse curve targets should have noDetectability = TRUE
  for (target in reverse_targets) {
    expect_true(data$targets$noDetectability[data$targets$targetName == target],
                info = sprintf("%s should have noDetectability = TRUE", target))
  }

  # Forward curve targets should NOT have noDetectability set (unless from XML modifiers)
  forward_targets <- data$targets$targetName[substr(data$targets$Curve_Quant, 1, 1) == "F"]
  for (target in forward_targets) {
    expect_false(data$targets$noDetectability[data$targets$targetName == target],
                 info = sprintf("%s should have noDetectability = FALSE", target))
  }
})

test_that("Reverse curve targets are excluded from detectability output", {
  input <- test_path("fixtures", "Analysis_INF250_Lot4_AQ_LC_R3_20241229.xml")
  data <- loadNULISAseq(input)

  reverse_targets <- NULISAseqR:::get_reverse_curve_targets(data$targets)

  # Reverse curve targets should NOT appear in detectability results
  for (target in reverse_targets) {
    expect_false(target %in% names(data$detectability$all$detectability),
                 info = sprintf("%s should not be in detectability results", target))
    expect_false(target %in% names(data$detectability$all$detectable),
                 info = sprintf("%s should not be in detectable results", target))
  }
})

test_that("Reverse curve targets retain NPQ values and NA LOD after noDetectability change", {
  input <- test_path("fixtures", "Analysis_INF250_Lot4_AQ_LC_R3_20241229.xml")
  data <- loadNULISAseq(input)

  # LOD should still be NA for reverse targets
  expect_true(is.na(data$lod$LOD["CRP"]))
  expect_true(is.na(data$lod$LODNPQ["CRP"]))

  # NPQ values should still be computed (not NA) for reverse targets
  expect_false(all(is.na(data$NPQ["CRP", ])))

  # Spot-check a known NPQ value to ensure no regression
  expect_equal(as.numeric(data$NPQ["CRP", 1]), 25.0840257805487, tolerance = 1e-4)
})

test_that("mergeNULISAseq detectability table has NA for reverse curve targets", {
  input <- test_path("fixtures", "Analysis_INF250_Lot4_AQ_LC_R3_20241229.xml")
  all_data <- importNULISAseq(files = input)
  merged <- all_data$merged

  detect_table <- merged$detectability
  reverse_targets <- c("CRP", "KNG1")

  # Detectability columns should be numeric
  for (col in colnames(detect_table)[-1]) {
    expect_true(is.numeric(detect_table[[col]]),
                info = sprintf("Column '%s' should be numeric", col))
  }

  for (target in reverse_targets) {
    row <- detect_table[detect_table$Target == target, ]
    expect_equal(nrow(row), 1,
                 info = sprintf("%s should appear in detectability table", target))
    # All non-Target columns should be NA (numeric)
    for (col in colnames(detect_table)[-1]) {
      expect_true(is.na(row[[col]]),
                  info = sprintf("%s %s should be NA", target, col))
    }
  }
})

test_that("label_detectability_for_display adds High Abundance for RC targets", {
  detect_table <- data.frame(
    Target = c("A", "B", "CRP"),
    `plasma (n = 10)` = c(80.0, 60.0, NA_real_),
    `csf (n = 5)` = c(70.0, NA_real_, NA_real_),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  display <- NULISAseqR:::label_detectability_for_display(detect_table, "CRP")

  expect_equal(display[display$Target == "CRP", "plasma (n = 10)"], "High Abundance")
  expect_true(is.na(display[display$Target == "CRP", "csf (n = 5)"]))
  expect_equal(display[display$Target == "A", "plasma (n = 10)"], "80.0")
})

test_that("XML noDetectability modifier targets have detectability computed but flag preserved", {
  # This fixture has noDetectability on AGER (a forward target, rare case)
  input <- test_path("fixtures", "detectability_P1_Tr03_typemCherry_noDetectability_AGER.xml")
  data <- loadNULISAseq(input, IPC = NULL, IC = NULL, SC = NULL)

  # AGER should have noDetectability flag from XML modifier
  expect_true(data$targets$noDetectability[data$targets$targetName == "AGER"])
  # But AGER should now have individual detectability computed (not excluded)
  # Non-RC noDetectability targets get computed but excluded from summary stats only
  expect_true("AGER" %in% names(data$detectability$all$detectability))
})

test_that("NULISAseqR:::get_noDetectability_targets returns flagged targets", {
  targets <- data.frame(
    targetName = c("A", "B", "C"),
    noDetectability = c(TRUE, FALSE, TRUE)
  )
  result <- NULISAseqR:::get_noDetectability_targets(targets)
  expect_equal(result, c("A", "C"))
})

test_that("NULISAseqR:::get_noDetectability_targets returns empty when column missing", {
  targets <- data.frame(targetName = c("A", "B"))
  expect_equal(NULISAseqR:::get_noDetectability_targets(targets), character(0))
})

test_that("NULISAseqR:::get_noDetectability_targets returns empty when none flagged", {
  targets <- data.frame(
    targetName = c("A", "B"),
    noDetectability = c(FALSE, FALSE)
  )
  expect_equal(NULISAseqR:::get_noDetectability_targets(targets), character(0))
})

test_that("NULISAseqR:::get_noDetectability_targets matches inline pattern from fixture", {
  input <- test_path("fixtures", "detectability_P1_Tr03_typemCherry_noDetectability_AGER.xml")
  data <- loadNULISAseq(input, IPC = NULL, IC = NULL, SC = NULL)

  # Helper should produce same result as inline pattern it replaces
  helper_result <- NULISAseqR:::get_noDetectability_targets(data$targets)
  inline_result <- data$targets$targetName[which(data$targets$noDetectability)]
  expect_equal(helper_result, inline_result)
  expect_true("AGER" %in% helper_result)
})
