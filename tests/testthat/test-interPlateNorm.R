test_that("interPlateNorm IPC-normalizes correctly with a normal multi-well IPC list", {
  D <- matrix(c(100, 200, 300, 400,
                 10,  20,  30,  40,
                  5,   6,   7,   8,
                 50,  60,  70,  80),
              nrow = 4, byrow = FALSE,
              dimnames = list(paste0("t", 1:4), c("IPC_1", "IPC_2", "s1", "s2")))

  res <- interPlateNorm(list(D), IPC_wells = list(c("IPC_1", "IPC_2")))
  expect_equal(dim(res$interNormData[[1]]), dim(D))
})

test_that("interPlateNorm still works when only one IPC well remains (regression for #697)", {
  # A stale/excluded well can leave exactly one IPC well behind on a plate;
  # matrix[, single_name] drops to a vector and apply(x, 1, median) then fails
  # with "dim(X) must have a positive length" unless the subset uses drop=FALSE.
  D <- matrix(c(100, 200, 300, 400,
                 10,  20,  30,  40,
                  5,   6,   7,   8,
                 50,  60,  70,  80),
              nrow = 4, byrow = FALSE,
              dimnames = list(paste0("t", 1:4), c("IPC_1", "IPC_2", "s1", "s2")))

  res <- expect_no_error(interPlateNorm(list(D), IPC_wells = list("IPC_1")))
  expect_equal(res$interNormData[[1]]["t1", "IPC_1"], 1e4)
})

test_that("interPlateNorm errors clearly when no IPC wells remain, instead of silently skipping normalization", {
  # If every IPC well on a plate has been excluded, IPC_wells[[i]] is
  # character(0). Without a guard, data_list[[i]][, character(0)] returns a
  # valid 0-column matrix, so IPC_factors_i comes back all-NA, gets replaced
  # with 1, and IPC normalization silently becomes a no-op instead of failing.
  D <- matrix(c(100, 200, 300, 400,
                 10,  20,  30,  40),
              nrow = 4, byrow = FALSE,
              dimnames = list(paste0("t", 1:4), c("s1", "s2")))

  expect_error(
    interPlateNorm(list(D), IPC_wells = list(character(0))),
    "No IPC wells remain"
  )
})

test_that("interPlateNorm errors clearly when no samples remain for intensity normalization (Bridge case)", {
  # IN_samples is used directly for Bridge-based normalization (skeleton.Rmd's
  # IC_Bridge call). If every Bridge well on a plate has been excluded,
  # IN_samples[[i]] is character(0); data_list[[i]][, character(0)] is a
  # valid 0-column matrix, so the plate's medians (and therefore its whole
  # normalized matrix) would silently come back all-NA instead of erroring.
  D <- matrix(c(100, 200, 300, 400,
                 10,  20,  30,  40),
              nrow = 4, byrow = FALSE,
              dimnames = list(paste0("t", 1:4), c("s1", "s2")))

  expect_error(
    interPlateNorm(list(D), IPC = FALSE, IN = TRUE, IN_samples = list(character(0))),
    "No samples remain"
  )
})

test_that("interPlateNorm errors clearly when IPC_wells names a well not present in the data", {
  # A stale IPC well name that survives some other exclusion path (or a plain
  # caller typo) previously reached data_list[[i]][, IPC_wells[[i]]] and
  # crashed with the generic "subscript out of bounds". Validate explicitly.
  D <- matrix(c(100, 200, 300, 400,
                 10,  20,  30,  40),
              nrow = 4, byrow = FALSE,
              dimnames = list(paste0("t", 1:4), c("IPC_1", "s1")))

  expect_error(
    interPlateNorm(list(D), IPC_wells = list(c("IPC_1", "IPC_2"))),
    "not found in the data"
  )
})

test_that("interPlateNorm IPC_method='mean' computes row means instead of erroring", {
  # rowMeans(IPC_cols, median, na.rm=TRUE) passed `median` positionally into
  # rowMeans' `dims` argument and errored for any input.
  D <- matrix(c(100, 200, 300, 400,
                 300, 400, 500, 600,
                  10,  20,  30,  40),
              nrow = 4, byrow = FALSE,
              dimnames = list(paste0("t", 1:4), c("IPC_1", "IPC_2", "s1")))

  res <- expect_no_error(interPlateNorm(list(D), IPC_wells = list(c("IPC_1", "IPC_2")), IPC_method = "mean"))
  expect_equal(res$interNormData[[1]]["t1", "IPC_1"], 1e4 / 2)
})
