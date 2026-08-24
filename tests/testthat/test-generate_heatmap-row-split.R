# Regression tests for the default row_split gating (issue #669).
#
# generate_heatmap() defaults row_split <- 2 when targetInfo is supplied. A
# numeric row_split needs a row dendrogram to cut; with cluster_rows = FALSE
# there is none, so ComplexHeatmap reinterprets 2 as a per-row grouping vector
# and aborts at draw time with:
#   Length or nrow of `row_split` should be same as nrow of `matrix`.
# The gate must set the default only when clustering is not explicitly disabled.

make_test_inputs <- function(n_targets = 10, n_samples = 8) {
  set.seed(669)
  mat <- matrix(
    stats::rnorm(n_targets * n_samples),
    nrow = n_targets,
    dimnames = list(
      paste0("Target", seq_len(n_targets)),
      paste0("Sample", seq_len(n_samples))
    )
  )
  targetInfo <- data.frame(
    Target = rownames(mat),
    stringsAsFactors = FALSE
  )
  sampleInfo <- data.frame(
    SampleName = colnames(mat),
    stringsAsFactors = FALSE
  )
  list(data = mat, targetInfo = targetInfo, sampleInfo = sampleInfo)
}

# Draw to a throwaway device so row_split validation (which fires at draw time)
# actually runs without leaving artifacts behind.
draw_offscreen <- function(h) {
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  }, add = TRUE)
  # draw() returns the initialized heatmap, needed before querying row_order().
  ComplexHeatmap::draw(h)
}

test_that("cluster_rows = FALSE with targetInfo builds and draws without error (#669)", {
  inp <- make_test_inputs()
  h <- generate_heatmap(
    data = inp$data,
    sampleInfo = inp$sampleInfo,
    sampleName_var = "SampleName",
    targetInfo = inp$targetInfo,
    cluster_rows = FALSE
  )
  expect_no_error(draw_offscreen(h$heatmap))
})

test_that("cluster_rows = TRUE with targetInfo still applies the default 2-way split", {
  inp <- make_test_inputs()
  h <- generate_heatmap(
    data = inp$data,
    sampleInfo = inp$sampleInfo,
    sampleName_var = "SampleName",
    targetInfo = inp$targetInfo,
    cluster_rows = TRUE
  )
  drawn <- draw_offscreen(h$heatmap)
  expect_equal(length(ComplexHeatmap::row_order(drawn)), 2)
})
