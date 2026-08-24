test_that("interCV returns a 0-column matrix when no replicate sets are given", {
  # Two plates, two targets. `samples` encodes NO replicate sets (all NA), so
  # unique_samples is empty. Regression for the `1:length(unique_samples)`
  # (1:0) bug that iterated into a 0-column matrix and threw
  # "subscript out of bounds".
  p1 <- matrix(c(10, 11, 12, 13), nrow = 2, dimnames = list(c("T1", "T2"), c("A", "B")))
  p2 <- matrix(c(20, 21, 22, 23), nrow = 2, dimnames = list(c("T1", "T2"), c("C", "D")))
  data_list <- list(p1, p2)
  lod <- lapply(data_list, function(m) {
    matrix(TRUE, nrow = nrow(m), ncol = ncol(m), dimnames = dimnames(m))
  })
  samples_none <- list(c(NA, NA), c(NA, NA))

  res_count <- expect_no_error(
    interCV(data_list, samples_none, aboveLOD = lod, useMean = FALSE, method = "count")
  )
  expect_equal(ncol(res_count), 0)
  expect_equal(rownames(res_count), c("T1", "T2"))

  # Same guarantee for the log2 branch (the second fixed loop).
  res_log2 <- expect_no_error(
    interCV(data_list, samples_none, aboveLOD = lod, useMean = FALSE, method = "log2")
  )
  expect_equal(ncol(res_log2), 0)
})

test_that("interCV still computes CV for a valid replicate set", {
  # Guards against the seq_along change regressing the happy path: one shared
  # replicate label ("S1") spanning both plates -> a single-column result.
  p1 <- matrix(c(10, 11, 12, 13), nrow = 2, dimnames = list(c("T1", "T2"), c("A", "B")))
  p2 <- matrix(c(20, 21, 22, 23), nrow = 2, dimnames = list(c("T1", "T2"), c("C", "D")))
  data_list <- list(p1, p2)
  lod <- lapply(data_list, function(m) {
    matrix(TRUE, nrow = nrow(m), ncol = ncol(m), dimnames = dimnames(m))
  })
  # Both columns of each plate are replicates of sample "S1" (two replicates
  # per plate), so useMean = TRUE collapses each plate to a mean column rather
  # than dropping a single-column slice to a bare vector.
  samples_one <- list(c("S1", "S1"), c("S1", "S1"))

  res <- interCV(data_list, samples_one, aboveLOD = lod, useMean = FALSE, method = "count")
  expect_equal(ncol(res), 1)
  expect_equal(colnames(res), "S1")
  expect_true(all(is.finite(res[, 1])))

  # Default useMean = TRUE path (collapses each plate to its replicate mean).
  res_mean <- interCV(data_list, samples_one, aboveLOD = lod, useMean = TRUE, method = "count")
  expect_equal(ncol(res_mean), 1)
  expect_true(all(is.finite(res_mean[, 1])))

  # log2 branch (the second loop touched by the seq_along fix) on a valid set.
  res_log2 <- interCV(data_list, samples_one, aboveLOD = lod, useMean = FALSE, method = "log2")
  expect_equal(ncol(res_log2), 1)
  expect_equal(colnames(res_log2), "S1")
  expect_true(all(is.finite(res_log2[, 1])))
})
