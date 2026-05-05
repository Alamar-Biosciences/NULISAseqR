# Tests for NULISAseqR:::format_detectability_stats_row() and NULISAseqR:::format_detectability_report()

# --- NULISAseqR:::format_detectability_stats_row() ---

test_that("NULISAseqR:::format_detectability_stats_row returns correct 1-row matrix for basic input", {
  vals <- c(A = 80, B = 60, C = 30)
  result <- NULISAseqR:::format_detectability_stats_row(vals, "all", 10)
  
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 8)
  expect_equal(colnames(result),
               c("type", "# samples", "mean", "sd", "median",
                 "min", "max", "# of detectable targets (%)"))
  expect_equal(unname(result[1, "type"]), "all")
  expect_equal(unname(result[1, "# samples"]), "10")
  # mean of 80, 60, 30 = 56.666... -> 56.7
  expect_equal(trimws(unname(result[1, "mean"])), "56.7")
  # 2 of 3 are > 50
  expect_match(result[1, "# of detectable targets (%)"], "2 / 3")
})

test_that("NULISAseqR:::format_detectability_stats_row handles all-NA input", {
  vals <- c(A = NA_real_, B = NA_real_)
  result <- NULISAseqR:::format_detectability_stats_row(vals, "PLASMA", 5)
  
  expect_equal(nrow(result), 1)
  expect_equal(unname(result[1, "type"]), "PLASMA")
  expect_equal(unname(result[1, "mean"]), "")
  expect_equal(unname(result[1, "sd"]), "")
  expect_equal(unname(result[1, "median"]), "")
  expect_equal(unname(result[1, "min"]), "")
  expect_equal(unname(result[1, "max"]), "")
  expect_match(result[1, "# of detectable targets (%)"], "0 / 0")
})

test_that("NULISAseqR:::format_detectability_stats_row handles empty vector", {
  vals <- numeric(0)
  result <- NULISAseqR:::format_detectability_stats_row(vals, "all", 0)
  
  expect_equal(nrow(result), 1)
  expect_equal(unname(result[1, "mean"]), "")
  expect_match(result[1, "# of detectable targets (%)"], "0 / 0")
})

test_that("NULISAseqR:::format_detectability_stats_row handles single value", {
  vals <- c(X = 75.0)
  result <- NULISAseqR:::format_detectability_stats_row(vals, "all", 1)
  
  expect_equal(nrow(result), 1)
  expect_equal(trimws(unname(result[1, "mean"])), "75.0")
  # sd of single value is NA
  expect_equal(trimws(unname(result[1, "sd"])), "NA")
  expect_match(result[1, "# of detectable targets (%)"], "1 / 1")
})

test_that("NULISAseqR:::format_detectability_stats_row preserves type_label for sample groups", {
  vals <- c(A = 90, B = 10)
  result <- NULISAseqR:::format_detectability_stats_row(vals, "SERUM", 3)
  expect_equal(unname(result[1, "type"]), "SERUM")
  expect_equal(unname(result[1, "# samples"]), "3")
})

# --- NULISAseqR:::format_detectability_report() ---

test_that("NULISAseqR:::format_detectability_report returns correct structure without sample groups", {
  detect_result <- list(
    all = list(
      sampleNumber = 5,
      detectability = c(T1 = 80, T2 = 60, T3 = 20),
      detectable = c(T1 = TRUE, T2 = TRUE, T3 = FALSE)
    )
  )
  
  result <- NULISAseqR:::format_detectability_report(detect_result)
  
  expect_true(is.list(result))
  expect_equal(names(result), c("summary", "targets"))
  
  # summary: 1 row (all), 8 columns
  expect_equal(nrow(result$summary), 1)
  expect_equal(ncol(result$summary), 8)
  expect_equal(unname(result$summary[1, "type"]), "all")
  
  # targets: 3 rows, 1 column named "all"
  expect_equal(nrow(result$targets), 3)
  expect_equal(ncol(result$targets), 1)
  expect_equal(colnames(result$targets), "all")
  expect_equal(rownames(result$targets), c("T1", "T2", "T3"))
})

test_that("NULISAseqR:::format_detectability_report handles sample groups", {
  detect_result <- list(
    all = list(
      sampleNumber = 6,
      detectability = c(T1 = 80, T2 = 60),
      detectable = c(T1 = TRUE, T2 = TRUE)
    ),
    sample_group = list(
      sampleNumber = list(PLASMA = 4, SERUM = 2),
      detectability = list(
        PLASMA = c(T1 = 90, T2 = 70),
        SERUM = c(T1 = 60, T2 = 40)
      ),
      detectable = list(
        PLASMA = c(T1 = TRUE, T2 = TRUE),
        SERUM = c(T1 = TRUE, T2 = FALSE)
      )
    )
  )
  
  result <- NULISAseqR:::format_detectability_report(detect_result)
  
  # summary: 3 rows (all, PLASMA, SERUM)
  expect_equal(nrow(result$summary), 3)
  expect_equal(unname(result$summary[1, "type"]), "all")
  expect_equal(unname(result$summary[2, "type"]), "PLASMA")
  expect_equal(unname(result$summary[3, "type"]), "SERUM")
  
  # targets: 2 rows, 3 columns
  expect_equal(ncol(result$targets), 3)
  expect_equal(colnames(result$targets), c("all", "PLASMA", "SERUM"))
})

test_that("NULISAseqR:::format_detectability_report excludes noDetectability targets from summary only", {
  detect_result <- list(
    all = list(
      sampleNumber = 5,
      detectability = c(T1 = 80, T2 = 60, AGER = 90),
      detectable = c(T1 = TRUE, T2 = TRUE, AGER = TRUE)
    )
  )
  
  result <- NULISAseqR:::format_detectability_report(
    detect_result,
    noDetectability_targets = "AGER"
  )
  
  # summary should be computed from T1 and T2 only (excludes AGER)
  # mean of 80, 60 = 70.0
  expect_equal(trimws(unname(result$summary[1, "mean"])), "70.0")
  
  # targets should include AGER
  expect_true("AGER" %in% rownames(result$targets))
  expect_equal(nrow(result$targets), 3)
})

test_that("NULISAseqR:::format_detectability_report appends High Abundance rows for RC targets", {
  detect_result <- list(
    all = list(
      sampleNumber = 5,
      detectability = c(T1 = 80, T2 = 60),
      detectable = c(T1 = TRUE, T2 = TRUE)
    ),
    sample_group = list(
      sampleNumber = list(PLASMA = 6, serum = 7),
      detectability = list(PLASMA = c(T1 = 90, T2 = 70),
                           serum = c(T1 = 93, T2 = 23)),
      detectable = list(PLASMA = c(T1 = TRUE, T2 = TRUE),
                        serum = c(T1 = TRUE, T2 = FALSE))
    )
  )
  
  result <- NULISAseqR:::format_detectability_report(
    detect_result,
    reverse_curve_targets = c("CRP", "KNG1")
  )
  
  # targets should have 4 rows: T1, T2, CRP, KNG1
  expect_equal(nrow(result$targets), 4)
  expect_equal(ncol(result$targets), 3)
  expect_true("CRP" %in% rownames(result$targets))
  expect_true("KNG1" %in% rownames(result$targets))
  expect_equal(unname(result$targets["CRP", "all"]), NA_character_)
  expect_equal(unname(result$targets["KNG1", "all"]), NA_character_)
  expect_equal(unname(result$targets["CRP", "PLASMA"]), "High Abundance")
  expect_equal(unname(result$targets["KNG1", "PLASMA"]), "High Abundance")
  expect_equal(unname(result$targets["CRP", "serum"]), "High Abundance")
  expect_equal(unname(result$targets["KNG1", "serum"]), "High Abundance")
})

test_that("NULISAseqR:::format_detectability_report target values are formatted correctly", {
  detect_result <- list(
    all = list(
      sampleNumber = 5,
      detectability = c(T1 = 66.66667),
      detectable = c(T1 = TRUE)
    )
  )
  
  result <- NULISAseqR:::format_detectability_report(detect_result)
  
  # Should be formatted with format(round(x, 1), nsmall=1) — no width/justify
  expect_equal(trimws(unname(result$targets["T1", "all"])), "66.7")
})

test_that("NULISAseqR:::format_detectability_report handles NA values from cross-plate aggregation", {
  combined_detect <- list(
    all = list(
      sampleNumber = 10,
      detectability = c(T1 = 80, T2 = NA),
      detectable = c(T1 = TRUE, T2 = NA)
    )
  )
  
  result <- NULISAseqR:::format_detectability_report(combined_detect)
  
  expect_equal(unname(result$targets["T2", "all"]), "")
  expect_equal(trimws(unname(result$targets["T1", "all"])), "80.0")
})

test_that("NULISAseqR:::format_detectability_report sorts targets case-insensitively", {
  detect_result <- list(
    all = list(
      sampleNumber = 5,
      detectability = c(pTau = 80, AGER = 60, mHTT = 40),
      detectable = c(pTau = TRUE, AGER = TRUE, mHTT = FALSE)
    )
  )
  
  result <- NULISAseqR:::format_detectability_report(detect_result)
  expect_equal(rownames(result$targets), c("AGER", "mHTT", "pTau"))
})

# --- plasma_serum_groups parameter tests ---

test_that("format_detectability_report uses plasma_serum_groups for High Abundance labeling", {
  detect_result <- list(
    all = list(
      sampleNumber = 10,
      detectability = c(T1 = 80, T2 = 60),
      detectable = c(T1 = TRUE, T2 = TRUE)
    ),
    sample_group = list(
      sampleNumber = list(HUMAN_PLASMA = 5, MOUSE_PLASMA = 3, MOUSE_CSF = 2),
      detectability = list(
        HUMAN_PLASMA = c(T1 = 90, T2 = 70),
        MOUSE_PLASMA = c(T1 = 85, T2 = 65),
        MOUSE_CSF = c(T1 = 50, T2 = 30)
      ),
      detectable = list(
        HUMAN_PLASMA = c(T1 = TRUE, T2 = TRUE),
        MOUSE_PLASMA = c(T1 = TRUE, T2 = TRUE),
        MOUSE_CSF = c(T1 = FALSE, T2 = FALSE)
      )
    )
  )

  # With plasma_serum_groups, only HUMAN_PLASMA and MOUSE_PLASMA get High Abundance
  result <- NULISAseqR:::format_detectability_report(
    detect_result,
    reverse_curve_targets = c("CRP"),
    plasma_serum_groups = c("HUMAN_PLASMA", "MOUSE_PLASMA")
  )

  expect_equal(unname(result$targets["CRP", "all"]), NA_character_)
  expect_equal(unname(result$targets["CRP", "HUMAN_PLASMA"]), "High Abundance")
  expect_equal(unname(result$targets["CRP", "MOUSE_PLASMA"]), "High Abundance")
  expect_equal(unname(result$targets["CRP", "MOUSE_CSF"]), NA_character_)
})

test_that("format_detectability_report without plasma_serum_groups uses strict regex", {
  detect_result <- list(
    all = list(
      sampleNumber = 10,
      detectability = c(T1 = 80),
      detectable = c(T1 = TRUE)
    ),
    sample_group = list(
      sampleNumber = list(HUMAN_PLASMA = 5, MOUSE_CSF = 5),
      detectability = list(
        HUMAN_PLASMA = c(T1 = 90),
        MOUSE_CSF = c(T1 = 70)
      ),
      detectable = list(
        HUMAN_PLASMA = c(T1 = TRUE),
        MOUSE_CSF = c(T1 = TRUE)
      )
    )
  )

  # Without plasma_serum_groups, strict regex ^(PLASMA|SERUM)$ should NOT match HUMAN_PLASMA
  result <- NULISAseqR:::format_detectability_report(
    detect_result,
    reverse_curve_targets = c("CRP")
  )

  expect_equal(unname(result$targets["CRP", "HUMAN_PLASMA"]), NA_character_)
  expect_equal(unname(result$targets["CRP", "MOUSE_CSF"]), NA_character_)
})

test_that("format_detectability_report strict regex matches exact PLASMA and SERUM", {
  detect_result <- list(
    all = list(
      sampleNumber = 10,
      detectability = c(T1 = 80),
      detectable = c(T1 = TRUE)
    ),
    sample_group = list(
      sampleNumber = list(PLASMA = 6, SERUM = 4),
      detectability = list(
        PLASMA = c(T1 = 90),
        SERUM = c(T1 = 70)
      ),
      detectable = list(
        PLASMA = c(T1 = TRUE),
        SERUM = c(T1 = TRUE)
      )
    )
  )

  # Without plasma_serum_groups, strict regex should match exact PLASMA and SERUM
  result <- NULISAseqR:::format_detectability_report(
    detect_result,
    reverse_curve_targets = c("CRP")
  )

  expect_equal(unname(result$targets["CRP", "PLASMA"]), "High Abundance")
  expect_equal(unname(result$targets["CRP", "SERUM"]), "High Abundance")
})

# --- label_detectability_for_display plasma_serum_groups tests ---

test_that("label_detectability_for_display uses plasma_serum_groups when provided", {
  detect_table <- data.frame(
    Target = c("A", "CRP"),
    `human_plasma (n = 10)` = c(80.0, NA_real_),
    `mouse_csf (n = 5)` = c(70.0, NA_real_),
    check.names = FALSE, stringsAsFactors = FALSE
  )

  display <- NULISAseqR:::label_detectability_for_display(
    detect_table, "CRP",
    plasma_serum_groups = c("HUMAN_PLASMA")
  )

  expect_equal(display[display$Target == "CRP", "human_plasma (n = 10)"], "High Abundance")
  expect_true(is.na(display[display$Target == "CRP", "mouse_csf (n = 5)"]))
})

test_that("label_detectability_for_display strict regex fallback does not match substrings", {
  detect_table <- data.frame(
    Target = c("A", "CRP"),
    `human_plasma (n = 10)` = c(80.0, NA_real_),
    `csf (n = 5)` = c(70.0, NA_real_),
    check.names = FALSE, stringsAsFactors = FALSE
  )

  # Without plasma_serum_groups, strict regex ^(plasma|serum)\s*\( should NOT match "human_plasma (n = 10)"
  display <- NULISAseqR:::label_detectability_for_display(detect_table, "CRP")

  expect_true(is.na(display[display$Target == "CRP", "human_plasma (n = 10)"]))
  expect_true(is.na(display[display$Target == "CRP", "csf (n = 5)"]))
})
