# Regression tests for sampleBoxplot() when the internal-control target is
# absent from the NPQ / normalized matrix, plus the filter_run_data() single-
# column drop=FALSE guard.
#
# Some panels drop the IC target row during NPQ normalization (e.g. mCherry),
# and importNULISAseq-merged data produces an NPQ matrix without the IC row.
# sampleBoxplot() draws an IC reference line via `data[IC_name, ]`; when the IC
# row is missing this was an out-of-bounds subscript that crashed the whole
# plot ("subscript out of bounds"). See sampleBoxplot.R internal sampleboxplot().

# Build a minimal, valid run object. `drop_ic_from_normed = TRUE` reproduces the
# bug condition: raw Data keeps the IC target row, but the normalized matrix
# does not. `ic` may be a single name or several (multi-IC panels).
make_run <- function(n_samples = 8, ic = "mCherry", drop_ic_from_normed = TRUE) {
  targets_data <- paste0("T", seq_len(10))
  all_targets  <- c(targets_data, ic)
  sample_names <- paste0("S", formatC(seq_len(n_samples), width = 2, flag = "0"))

  set.seed(1)
  raw <- matrix(abs(rnorm(length(all_targets) * n_samples, 50, 10)),
                nrow = length(all_targets), dimnames = list(all_targets, sample_names))
  normed_targets <- if (drop_ic_from_normed) targets_data else all_targets
  normed <- matrix(abs(rnorm(length(normed_targets) * n_samples, 5, 1)),
                   nrow = length(normed_targets), dimnames = list(normed_targets, sample_names))

  samples <- data.frame(
    sampleName    = sample_names,
    sampleType    = "Sample",
    AUTO_WELLROW  = "A",
    AUTO_WELLCOL  = seq_len(n_samples),
    stringsAsFactors = FALSE
  )
  list(
    samples          = samples,
    # Real target tables carry several columns; keep >1 so row-subsetting in
    # filter_run_data() doesn't drop the data.frame to a vector.
    targets          = data.frame(
      targetName = all_targets,
      targetType = ifelse(all_targets %in% ic, "control", "analyte"),
      stringsAsFactors = FALSE
    ),
    Data             = raw,
    normed           = list(interNormData = list(normed)),
    plateID          = "Plate_01",
    IC               = ic,
    ExecutionDetails = list(Panel = "TestPanel")
  )
}

test_that("sampleBoxplot(metrics='normed') does not crash when IC row is absent from normed", {
  run <- make_run(drop_ic_from_normed = TRUE)
  expect_false(run$IC %in% rownames(run$normed$interNormData[[1]]))  # bug precondition

  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off())
  expect_no_error(
    sampleBoxplot(runs = list(run), metrics = "normed", output_dir = NULL)
  )
})

test_that("sampleBoxplot(metrics='all') does not crash when IC row is absent from normed", {
  run <- make_run(drop_ic_from_normed = TRUE)
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off())
  expect_no_error(
    sampleBoxplot(runs = list(run), metrics = "all", output_dir = NULL)
  )
})

test_that("sampleBoxplot still works normally when the IC row IS present", {
  run <- make_run(drop_ic_from_normed = FALSE)
  expect_true(run$IC %in% rownames(run$normed$interNormData[[1]]))

  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off())
  expect_no_error(
    sampleBoxplot(runs = list(run), metrics = "normed", output_dir = NULL)
  )
})

test_that("sampleBoxplot handles MULTIPLE internal-control rows (colMeans branch)", {
  # Two IC targets both present in the normed matrix -> the length>1 colMeans
  # path in sampleboxplot() must run without error.
  run <- make_run(ic = c("mCherry", "mCherry2"), drop_ic_from_normed = FALSE)
  expect_true(all(run$IC %in% rownames(run$normed$interNormData[[1]])))
  expect_length(run$IC, 2L)

  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off())
  expect_no_error(
    sampleBoxplot(runs = list(run), metrics = "normed", output_dir = NULL)
  )
})

test_that("filter_run_data keeps single-column samples/targets as data.frames (drop=FALSE)", {
  # Guards the drop=FALSE fix: a 1-column table row-subset without drop=FALSE
  # collapses to a vector, so nrow() returns NULL and the downstream
  # nrow(...)==0 check throws "missing value where TRUE/FALSE needed".
  run <- list(
    samples = data.frame(sampleName = c("S1", "S2", "S3"), stringsAsFactors = FALSE),
    targets = data.frame(targetName = c("T1", "T2"), stringsAsFactors = FALSE),
    IPC = NULL, NC = NULL, SC = NULL
  )
  res <- NULISAseqR:::filter_run_data(run)  # internal, non-exported
  expect_true(is.data.frame(res$samples))
  expect_true(is.data.frame(res$targets))
  expect_equal(nrow(res$samples), 3L)
  expect_equal(nrow(res$targets), 2L)
})
