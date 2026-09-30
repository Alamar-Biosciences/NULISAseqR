# Regression tests for one-target inputs, where base R subsetting and apply()
# silently drop a matrix to a vector, and for the case where every per-target
# model fit fails.

make_counts <- function() {
  set.seed(1)
  m <- matrix(rexp(30 * 3, 1 / 100), nrow = 3,
              dimnames = list(c("IL6", "TNF", "CRP"), paste0("S", 1:30)))
  colnames(m)[1:3] <- c("NC1", "NC2", "NC3")
  m[, 1:3] <- 5
  m
}

# Runs targetBoxplot() on a null device and returns the stats of the boxplot
# it drew, so tests can count boxes and read their labels.
capture_target_boxplot <- function(...) {
  drawn <- NULL
  local_mocked_bindings(boxplot = function(x, ...) {
    drawn <<- graphics::boxplot(x, ...)
    invisible(drawn)
  })
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  suppressMessages(capture.output(targetBoxplot(...)))
  drawn
}

test_that("targetBoxplot draws one labelled box for a single target", {
  m <- make_counts()
  one <- m["IL6", , drop = FALSE]

  vertical <- capture_target_boxplot(one)
  expect_equal(ncol(vertical$stats), 1)
  expect_equal(vertical$names, "IL6")

  horizontal <- capture_target_boxplot(one, horizontal = TRUE)
  expect_equal(ncol(horizontal$stats), 1)
  expect_equal(horizontal$names, "IL6")
})

test_that("targetBoxplot keeps a matrix when excludeTargets leaves one target", {
  m <- make_counts()
  res <- capture_target_boxplot(m, excludeTargets = c("TNF", "CRP"))
  expect_equal(ncol(res$stats), 1)
  expect_equal(res$names, "IL6")
})

test_that("targetBoxplot keeps a matrix when excludeSamples leaves one sample", {
  m <- make_counts()
  res <- capture_target_boxplot(m, excludeSamples = setdiff(colnames(m), "S10"))
  expect_equal(ncol(res$stats), 3)
  expect_setequal(res$names, rownames(m))
})

test_that("targetBoxplot with subtractLOD draws one box for a single target", {
  m <- make_counts()
  nc <- c("NC1", "NC2", "NC3")
  one <- m["IL6", , drop = FALSE]

  res <- capture_target_boxplot(one, subtractLOD = TRUE, blanks = nc)
  expect_equal(ncol(res$stats), 1)
  expect_equal(res$names, "IL6")

  res_h <- capture_target_boxplot(one, subtractLOD = TRUE, blanks = nc,
                                  horizontal = TRUE, sortBy = "detect")
  expect_equal(ncol(res_h$stats), 1)
  expect_equal(res_h$names, "IL6")
})

test_that("targetBoxplot still orders multiple targets by median", {
  m <- make_counts()
  res <- capture_target_boxplot(m)
  expect_equal(ncol(res$stats), 3)
  expect_true(all(diff(res$stats[3, ]) <= 0))
  res_h <- capture_target_boxplot(m, horizontal = TRUE)
  expect_equal(res_h$names, rev(res$names))
})

test_that("lod returns an aboveLOD matrix for a single target", {
  m <- make_counts()
  out <- suppressMessages(lod(m["IL6", , drop = FALSE], blanks = c("NC1", "NC2", "NC3")))
  expect_true(is.matrix(out$aboveLOD))
  expect_equal(dim(out$aboveLOD), c(1L, ncol(m)))
  expect_equal(dimnames(out$aboveLOD), dimnames(m["IL6", , drop = FALSE]))
})

test_that("lod aboveLOD compares each target to its own LOD", {
  m <- make_counts()
  out <- suppressMessages(lod(m, blanks = c("NC1", "NC2", "NC3")))
  expected <- sweep(m, 1, out$LOD, ">")
  expect_identical(out$aboveLOD, expected)
})

make_model_inputs <- function() {
  set.seed(2)
  data <- log2(make_counts()[, -(1:3)] + 1)
  info <- data.frame(sampleName = colnames(data),
                     y = rnorm(ncol(data)),
                     age = rnorm(ncol(data)),
                     sid = rep(1:9, 3),
                     one_level = "a",
                     stringsAsFactors = FALSE)
  list(data = data, info = info)
}

test_that("lmNULISAseq_predict handles a single target_subset", {
  d <- make_model_inputs()
  out <- lmNULISAseq_predict(d$data, d$info, "sampleName", "y", "age",
                             target_subset = "IL6")
  expect_equal(nrow(out$modelStats), 1)
  expect_equal(out$modelStats$target, "IL6")
  expect_true("age_pval_FDR" %in% names(out$modelStats))
  expect_equal(out$modelStats$age_pval_FDR, out$modelStats$age_pval_unadj)
})

test_that("lmNULISAseq_predict reports the first error when every fit fails", {
  d <- make_model_inputs()
  expect_error(
    capture.output(
      lmNULISAseq_predict(d$data, d$info, "sampleName", "y", "one_level")
    ),
    "All targets failed model fitting.*IL6: contrasts can be applied only to factors with 2 or more levels"
  )
})

test_that("lmNULISAseq handles a single target", {
  d <- make_model_inputs()
  out <- lmNULISAseq(d$data["IL6", , drop = FALSE], d$info, "sampleName", "age")
  expect_equal(nrow(out$modelStats), 1)
  expect_equal(out$modelStats$target, "IL6")
  expect_true("age_pval_bonf" %in% names(out$modelStats))
})

test_that("lmNULISAseq reports the first error when every fit fails", {
  d <- make_model_inputs()
  expect_error(
    suppressWarnings(
      lmNULISAseq(d$data, d$info, "sampleName", "one_level",
                  analysis_context = "ctx")
    ),
    "^ctx: All targets failed model fitting.*IL6: contrasts can be applied"
  )
})

test_that("lmerNULISAseq_predict handles a single target and all-failed fits", {
  skip_if_not_installed("lmerTest")
  d <- make_model_inputs()
  out <- suppressMessages(
    lmerNULISAseq_predict(d$data, d$info, "sampleName", "y", "age", "(1|sid)",
                          target_subset = "IL6")
  )
  expect_equal(nrow(out$modelStats), 1)
  expect_equal(out$modelStats$target, "IL6")

  expect_error(
    capture.output(
      lmerNULISAseq_predict(d$data, d$info, "sampleName", "y", "one_level", "(1|sid)")
    ),
    "All targets failed model fitting.*IL6:"
  )
})

test_that("glmNULISAseq_predict handles a single target_subset", {
  d <- make_model_inputs()
  d$info$y <- rep(c(0, 1), length.out = nrow(d$info))
  out <- suppressMessages(
    glmNULISAseq_predict(d$data, d$info, "sampleName", "y", "age",
                         target_subset = "IL6")
  )
  expect_equal(nrow(out$modelStats), 1)
  expect_equal(out$modelStats$target, "IL6")
})

test_that("safe_extract_matrix returns a zero-column matrix when nothing was fitted", {
  stats_list <- list(IL6 = NULL, TNF = NULL)
  out <- safe_extract_matrix(stats_list, "coefs", NULL)
  expect_true(is.matrix(out))
  expect_equal(dim(out), c(2L, 0L))
  expect_equal(rownames(out), c("IL6", "TNF"))
})

test_that("p_adjust_columns keeps a one-row matrix", {
  p <- matrix(c(0.01, 0.2), nrow = 1, dimnames = list("IL6", c("a", "b")))
  out <- p_adjust_columns(p, "BH")
  expect_equal(dim(out), c(1L, 2L))
  expect_equal(dimnames(out), dimnames(p))
  p3 <- matrix(c(0.01, 0.02, 0.03, 0.5, 0.2, 0.1), nrow = 3)
  expect_equal(p_adjust_columns(p3, "BH"), apply(p3, 2, p.adjust, method = "BH"))
})
