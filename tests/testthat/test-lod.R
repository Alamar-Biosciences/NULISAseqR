test_that("lod computes normally with 2+ blank wells", {
  D <- matrix(c(1, 2, 3, 4,
                1, 2, 3, 4,
                10, 20, 30, 40,
                5, 6, 7, 8),
              nrow = 4, dimnames = list(paste0("t", 1:4), c("NC_1", "NC_2", "s1", "s2")))

  res <- lod(D, blanks = c("NC_1", "NC_2"))
  expect_equal(dim(res$aboveLOD), dim(D))
})

test_that("lod errors clearly with fewer than 2 blank wells, instead of crashing on drop=TRUE (regression for #697)", {
  # A stale/excluded NC well can leave exactly one (or zero) blank wells behind
  # on a plate. matrix[, single_name] drops to a vector, and apply(vector, 1, ...)
  # then fails with the opaque "dim(X) must have a positive length" -- and an
  # LOD from a single blank is meaningless (mean +/- MAD of one value) even if
  # it didn't crash.
  D <- matrix(c(1, 2, 3, 4,
                10, 20, 30, 40,
                5, 6, 7, 8),
              nrow = 4, dimnames = list(paste0("t", 1:4), c("NC_1", "s1", "s2")))

  expect_error(lod(D, blanks = "NC_1"), "At least 2 NC/blank wells")
  expect_error(lod(D, blanks = character(0)), "At least 2 NC/blank wells")
})
