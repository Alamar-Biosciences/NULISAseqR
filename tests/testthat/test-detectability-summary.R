# Tests for detectability_summary() formatted output

# Helper to build a minimal mock run with targets and detectability
mock_run <- function(detect_all, detect_groups = NULL, sample_numbers = NULL,
                     rc_targets = character(0), noDetect_targets = character(0)) {
  # Include all target names: those in detectability vectors + RC/noDetect targets
  all_target_names <- unique(c(names(detect_all), rc_targets, noDetect_targets))
  targets_df <- data.frame(
    targetName = all_target_names,
    Curve_Quant = ifelse(all_target_names %in% rc_targets, "R1", "S1"),
    noDetectability = all_target_names %in% noDetect_targets,
    stringsAsFactors = FALSE
  )

  detectability_obj <- list(
    all = list(
      sampleNumber = if (is.null(sample_numbers)) 10L else sample_numbers[["all"]],
      detectability = detect_all,
      detectable = detect_all > 50
    )
  )

  if (!is.null(detect_groups)) {
    detectability_obj$sample_group <- list(
      sampleNumber = sample_numbers[names(sample_numbers) != "all"],
      detectability = detect_groups,
      detectable = lapply(detect_groups, function(x) x > 50)
    )
  }

  list(
    RunSummary = NULL,  # prevent single-run wrapping when in a list
    targets = targets_df,
    detectability = detectability_obj
  )
}

# --- Basic structure tests ---

test_that("detectability_summary returns run_summary and run_targets when format=TRUE", {
  run1 <- mock_run(c(T1 = 80, T2 = 60, T3 = 20))

  result <- detectability_summary(list("Plate 01" = run1), format = TRUE,
                                  exclude_noDetect_targets = FALSE)

  expect_true("run_summary" %in% names(result))
  expect_true("run_targets" %in% names(result))
  expect_true("Plate 01" %in% names(result$run_summary))
  expect_true("Plate 01" %in% names(result$run_targets))

  # summary: 1 row (all), 8 columns
  expect_equal(nrow(result$run_summary[["Plate 01"]]), 1)
  expect_equal(ncol(result$run_summary[["Plate 01"]]), 8)

  # targets: 3 rows, 1 column (all)
  expect_equal(nrow(result$run_targets[["Plate 01"]]), 3)
  expect_equal(ncol(result$run_targets[["Plate 01"]]), 1)
})

test_that("detectability_summary omits run_summary/run_targets when format=FALSE", {
  run1 <- mock_run(c(T1 = 80, T2 = 60))

  result <- detectability_summary(list("Plate 01" = run1), format = FALSE,
                                  exclude_noDetect_targets = FALSE)

  expect_false("run_summary" %in% names(result))
  expect_false("run_targets" %in% names(result))
  # original fields still present
  expect_true("all" %in% names(result))
  expect_true("run_detectability" %in% names(result))
})

# --- Reverse curve targets ---

test_that("RC targets appear as High Abundance in target tables", {
  run1 <- mock_run(
    detect_all = c(T1 = 80, T2 = 60),
    detect_groups = list(PLASMA = c(T1 = 80, T2 = 60)),
    sample_numbers = list(all = 10L, PLASMA = 10L),
    rc_targets = "CRP",
    noDetect_targets = "CRP"
  )
  result <- detectability_summary(list("Plate 01" = run1), format = TRUE,
                                  exclude_noDetect_targets = TRUE)
  targets_mat <- result$run_targets[["Plate 01"]]
  expect_true("CRP" %in% rownames(targets_mat))
  expect_equal(unname(targets_mat["CRP", "all"]), NA_character_)
  expect_equal(unname(targets_mat["CRP", "PLASMA"]), "High Abundance")
})

test_that("reverse_curve_targets field is populated", {
  run1 <- mock_run(
    detect_all = c(T1 = 80),
    rc_targets = "CRP",
    noDetect_targets = "CRP"
  )

  expect_message(
    detectability_summary(list("P1" = run1), format = FALSE,
                          exclude_noDetect_targets = TRUE),
    "Reverse curve targets excluded"
  )

  result <- suppressMessages(
    detectability_summary(list("P1" = run1), format = FALSE,
                          exclude_noDetect_targets = TRUE)
  )
  expect_equal(result$reverse_curve_targets, "CRP")
})

# --- Non-RC noDetectability targets ---

test_that("non-RC noDetectability targets excluded from summary stats but present in target tables", {
  # AGER is noDetectability but not RC
  run1 <- mock_run(
    detect_all = c(T1 = 80, T2 = 60, AGER = 90),
    noDetect_targets = "AGER"
  )

  result <- detectability_summary(list("Plate 01" = run1), format = TRUE,
                                  exclude_noDetect_targets = TRUE)

  # AGER should be in target table
  expect_true("AGER" %in% rownames(result$run_targets[["Plate 01"]]))

  # summary stats should be based on T1 and T2 only (mean = 70.0)
  expect_equal(trimws(unname(result$run_summary[["Plate 01"]][1, "mean"])), "70.0")
})

test_that("noDetect_nonRC_targets field is populated", {
  run1 <- mock_run(
    detect_all = c(T1 = 80, AGER = 90),
    noDetect_targets = "AGER"
  )

  expect_message(
    detectability_summary(list("P1" = run1), format = FALSE,
                          exclude_noDetect_targets = TRUE),
    "Rare case targets excluded"
  )

  result <- suppressMessages(
    detectability_summary(list("P1" = run1), format = FALSE,
                          exclude_noDetect_targets = TRUE)
  )
  expect_equal(result$noDetect_nonRC_targets, "AGER")
})

# --- Message output ---

test_that("messages list RC and rare case targets", {
  run1 <- mock_run(
    detect_all = c(T1 = 80, AGER = 90),
    rc_targets = "CRP",
    noDetect_targets = c("CRP", "AGER")
  )

  expect_message(
    detectability_summary(list("P1" = run1), format = FALSE,
                          exclude_noDetect_targets = TRUE),
    "Reverse curve targets excluded from detectability: CRP"
  )

  expect_message(
    detectability_summary(list("P1" = run1), format = FALSE,
                          exclude_noDetect_targets = TRUE),
    "Rare case targets excluded from summary statistics: AGER"
  )
})

test_that("no messages when exclude_noDetect_targets = FALSE", {
  run1 <- mock_run(
    detect_all = c(T1 = 80),
    rc_targets = "CRP",
    noDetect_targets = "CRP"
  )

  expect_no_message(
    detectability_summary(list("P1" = run1), format = FALSE,
                          exclude_noDetect_targets = FALSE)
  )
})

# --- Multi-plate Overall tables ---

test_that("multi-plate produces Overall entry in run_summary and run_targets", {
  run1 <- mock_run(c(T1 = 80, T2 = 60))
  run2 <- mock_run(c(T1 = 70, T2 = 50))

  result <- detectability_summary(list("P1" = run1, "P2" = run2),
                                  format = TRUE, exclude_noDetect_targets = FALSE)

  expect_true("Overall" %in% names(result$run_summary))
  expect_true("Overall" %in% names(result$run_targets))
  expect_true("P1" %in% names(result$run_summary))
  expect_true("P2" %in% names(result$run_summary))
})

test_that("single plate does not produce Overall entry", {
  run1 <- mock_run(c(T1 = 80, T2 = 60))

  result <- detectability_summary(list("P1" = run1),
                                  format = TRUE, exclude_noDetect_targets = FALSE)

  expect_false("Overall" %in% names(result$run_summary))
  expect_false("Overall" %in% names(result$run_targets))
})

# --- Backward compatibility ---

test_that("original fields (all, sample_group, run_detectability) are unchanged", {
  run1 <- mock_run(
    detect_all = c(T1 = 80, T2 = 60),
    detect_groups = list(PLASMA = c(T1 = 90, T2 = 70), SERUM = c(T1 = 60, T2 = 40)),
    sample_numbers = list(all = 10L, PLASMA = 6L, SERUM = 4L)
  )

  result <- detectability_summary(list("P1" = run1),
                                  format = TRUE, exclude_noDetect_targets = FALSE)

  expect_true("all" %in% names(result))
  expect_true("run_detectability" %in% names(result))
  expect_true("sample_group" %in% names(result))
  expect_equal(result$all$sampleNumber, 10L)
  expect_equal(names(result$sample_group$detectability), c("PLASMA", "SERUM"))
})

test_that("sample groups appear in formatted tables", {
  run1 <- mock_run(
    detect_all = c(T1 = 80, T2 = 60),
    detect_groups = list(PLASMA = c(T1 = 90, T2 = 70), SERUM = c(T1 = 60, T2 = 40)),
    sample_numbers = list(all = 10L, PLASMA = 6L, SERUM = 4L)
  )

  result <- detectability_summary(list("P1" = run1),
                                  format = TRUE, exclude_noDetect_targets = FALSE)

  # summary should have 3 rows: all, PLASMA, SERUM
  expect_equal(nrow(result$run_summary[["P1"]]), 3)
  expect_equal(unname(result$run_summary[["P1"]][2, "type"]), "PLASMA")

  # targets should have 3 columns: all, PLASMA, SERUM
  expect_equal(ncol(result$run_targets[["P1"]]), 3)
  expect_equal(colnames(result$run_targets[["P1"]]), c("all", "PLASMA", "SERUM"))
})

# --- exclude_targets parameter ---

test_that("exclude_targets as character vector removes targets from all plates", {
  run1 <- mock_run(c(T1 = 80, T2 = 60, T3 = 20))
  run2 <- mock_run(c(T1 = 70, T2 = 50, T3 = 30))

  result <- detectability_summary(list("P1" = run1, "P2" = run2),
                                  exclude_targets = "T3",
                                  format = TRUE, exclude_noDetect_targets = FALSE)

  # T3 should be absent from aggregated output
  expect_false("T3" %in% names(result$all$detectability))
  # T1 and T2 should remain
  expect_true("T1" %in% names(result$all$detectability))
  expect_true("T2" %in% names(result$all$detectability))
  # per-plate formatted targets should also exclude T3
  expect_false("T3" %in% rownames(result$run_targets[["P1"]]))
  expect_false("T3" %in% rownames(result$run_targets[["P2"]]))
})

test_that("exclude_targets as list applies per-plate exclusion", {
  run1 <- mock_run(c(T1 = 80, T2 = 60, T3 = 20))
  run2 <- mock_run(c(T1 = 70, T2 = 50, T3 = 30))

  # exclude T3 from plate 1 only, T2 from plate 2 only
  result <- detectability_summary(list("P1" = run1, "P2" = run2),
                                  exclude_targets = list(c("T3"), c("T2")),
                                  format = TRUE, exclude_noDetect_targets = FALSE)

  # T3 excluded from P1, T2 excluded from P2
  expect_false("T3" %in% rownames(result$run_targets[["P1"]]))
  expect_false("T2" %in% rownames(result$run_targets[["P2"]]))
  # T3 still present on P2, T2 still present on P1
  expect_true("T3" %in% rownames(result$run_targets[["P2"]]))
  expect_true("T2" %in% rownames(result$run_targets[["P1"]]))
  # Overall should still have all three (each target present on at least one plate)
  expect_true("T1" %in% names(result$all$detectability))
  expect_true("T2" %in% names(result$all$detectability))
  expect_true("T3" %in% names(result$all$detectability))
})

test_that("exclude_targets with nonexistent target is a no-op", {
  run1 <- mock_run(c(T1 = 80, T2 = 60))

  result <- detectability_summary(list("P1" = run1),
                                  exclude_targets = "NONEXISTENT",
                                  format = FALSE, exclude_noDetect_targets = FALSE)

  expect_equal(length(result$all$detectability), 2)
  expect_true("T1" %in% names(result$all$detectability))
  expect_true("T2" %in% names(result$all$detectability))
})

test_that("exclude_targets removing all targets produces empty output", {
  run1 <- mock_run(c(T1 = 80, T2 = 60))

  result <- detectability_summary(list("P1" = run1),
                                  exclude_targets = c("T1", "T2"),
                                  format = FALSE, exclude_noDetect_targets = FALSE)

  expect_equal(length(result$all$detectability), 0)
  expect_equal(length(result$all$detectable), 0)
})
