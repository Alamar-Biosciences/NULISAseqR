# Test fill_predictors function
test_that("fill_predictors works with complete predictor set", {
  x <- c(a = 1, b = 2, c = 3)
  all <- c("a", "b", "c")
  
  result <- fill_predictors(x, all)
  
  expect_equal(result, x)
  expect_equal(names(result), all)
})

test_that("fill_predictors fills missing predictors with NA", {
  x <- c(a = 1, b = 2)
  all <- c("a", "b", "c", "d")
  
  result <- fill_predictors(x, all)
  
  expect_length(result, 4)
  expect_equal(names(result), all)
  expect_equal(result["a"], c(a = 1))
  expect_equal(result["b"], c(b = 2))
  expect_true(is.na(result["c"]))
  expect_true(is.na(result["d"]))
})

test_that("fill_predictors handles all missing predictors", {
  x <- c(a = 1, b = 2)
  all <- c("c", "d", "e")
  
  result <- fill_predictors(x, all)
  
  expect_length(result, 5)
  expect_equal(result["a"], c(a = 1))
  expect_equal(result["b"], c(b = 2))
  expect_true(all(is.na(result[c("c", "d", "e")])))
})

test_that("fill_predictors handles empty x", {
  x <- numeric(0)
  all <- c("a", "b", "c")
  
  result <- fill_predictors(x, all)
  
  expect_length(result, 3)
  expect_true(all(is.na(result)))
  expect_equal(names(result), all)
})

test_that("fill_predictors handles empty all", {
  x <- c(a = 1, b = 2)
  all <- character(0)
  
  result <- fill_predictors(x, all)
  
  expect_equal(result, x)
})

test_that("fill_predictors preserves order of original x", {
  x <- c(z = 3, a = 1, m = 2)
  all <- c("a", "b", "m", "z")
  
  result <- fill_predictors(x, all)
  
  # Original elements should come first in their original order
  expect_equal(names(result)[1:3], c("z", "a", "m"))
})

test_that("fill_predictors handles NULL x", {
  x <- NULL
  all <- c("a", "b", "c")
  
  # This might fail - test current behavior or fix function
  expect_error(fill_predictors(x, all), NA)  # or expect specific behavior
})


# Test safe_extract_matrix function
test_that("safe_extract_matrix extracts valid data correctly", {
  stats_list <- list(
    gene1 = list(coefs = c(a = 1, b = 2), t_vals = c(a = 1.5, b = 2.5), p_vals = c(a = 0.05, b = 0.01)),
    gene2 = list(coefs = c(a = 3, b = 4), t_vals = c(a = 3.5, b = 4.5), p_vals = c(a = 0.001, b = 0.0001))
  )
  all_predictors <- c("a", "b")
  
  result <- safe_extract_matrix(stats_list, "coefs", all_predictors)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
  expect_equal(rownames(result), c("gene1", "gene2"))
  expect_equal(colnames(result), c("a", "b"))
  expect_equal(result["gene1", "a"], 1)
  expect_equal(result["gene2", "b"], 4)
})

test_that("safe_extract_matrix handles missing predictors", {
  stats_list <- list(
    gene1 = list(coefs = c(a = 1), t_vals = c(a = 1.5), p_vals = c(a = 0.05)),
    gene2 = list(coefs = c(b = 4), t_vals = c(b = 4.5), p_vals = c(b = 0.0001))
  )
  all_predictors <- c("a", "b", "c")
  
  result <- safe_extract_matrix(stats_list, "coefs", all_predictors)
  
  expect_equal(ncol(result), 3)
  expect_true(is.na(result["gene1", "b"]))
  expect_true(is.na(result["gene1", "c"]))
  expect_true(is.na(result["gene2", "a"]))
  expect_true(is.na(result["gene2", "c"]))
})

test_that("safe_extract_matrix handles NULL list elements", {
  stats_list <- list(
    gene1 = list(coefs = c(a = 1, b = 2), t_vals = c(a = 1.5, b = 2.5), p_vals = c(a = 0.05, b = 0.01)),
    gene2 = NULL,
    gene3 = list(coefs = c(a = 3, b = 4), t_vals = c(a = 3.5, b = 4.5), p_vals = c(a = 0.001, b = 0.0001))
  )
  all_predictors <- c("a", "b")
  
  result <- safe_extract_matrix(stats_list, "coefs", all_predictors)
  
  expect_equal(nrow(result), 3)
  expect_true(all(is.na(result["gene2", ])))
  expect_equal(result["gene1", "a"], 1)
  expect_equal(result["gene3", "b"], 4)
})

test_that("safe_extract_matrix handles missing fields", {
  stats_list <- list(
    gene1 = list(coefs = c(a = 1, b = 2), t_vals = c(a = 1.5, b = 2.5)),
    gene2 = list(coefs = c(a = 3, b = 4), t_vals = c(a = 3.5, b = 4.5), p_vals = c(a = 0.001, b = 0.0001))
  )
  all_predictors <- c("a", "b")
  
  # gene1 is missing p_vals
  result <- safe_extract_matrix(stats_list, "p_vals", all_predictors)
  
  expect_equal(nrow(result), 2)
  expect_true(all(is.na(result["gene1", ])))
  expect_equal(result["gene2", "a"], 0.001)
})

test_that("safe_extract_matrix handles NULL field values", {
  stats_list <- list(
    gene1 = list(coefs = c(a = 1, b = 2), t_vals = NULL, p_vals = c(a = 0.05, b = 0.01)),
    gene2 = list(coefs = c(a = 3, b = 4), t_vals = c(a = 3.5, b = 4.5), p_vals = c(a = 0.001, b = 0.0001))
  )
  all_predictors <- c("a", "b")
  
  result <- safe_extract_matrix(stats_list, "t_vals", all_predictors)
  
  expect_true(all(is.na(result["gene1", ])))
  expect_equal(result["gene2", "a"], 3.5)
})

test_that("safe_extract_matrix handles empty stats_list", {
  stats_list <- list()
  all_predictors <- c("a", "b")
  
  result <- safe_extract_matrix(stats_list, "coefs", all_predictors)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 0)
  expect_equal(ncol(result), 2)
})

test_that("safe_extract_matrix handles single element", {
  stats_list <- list(
    gene1 = list(coefs = c(a = 1, b = 2), t_vals = c(a = 1.5, b = 2.5), p_vals = c(a = 0.05, b = 0.01))
  )
  all_predictors <- c("a", "b")
  
  result <- safe_extract_matrix(stats_list, "coefs", all_predictors)
  
  expect_equal(nrow(result), 1)
  expect_equal(rownames(result), "gene1")
})

test_that("safe_extract_matrix handles all fields (coefs, t_vals, p_vals)", {
  stats_list <- list(
    gene1 = list(coefs = c(a = 1, b = 2), t_vals = c(a = 1.5, b = 2.5), p_vals = c(a = 0.05, b = 0.01)),
    gene2 = list(coefs = c(a = 3, b = 4), t_vals = c(a = 3.5, b = 4.5), p_vals = c(a = 0.001, b = 0.0001))
  )
  all_predictors <- c("a", "b")
  
  coef_result <- safe_extract_matrix(stats_list, "coefs", all_predictors)
  tval_result <- safe_extract_matrix(stats_list, "t_vals", all_predictors)
  pval_result <- safe_extract_matrix(stats_list, "p_vals", all_predictors)
  
  expect_equal(dim(coef_result), dim(tval_result))
  expect_equal(dim(tval_result), dim(pval_result))
  expect_equal(rownames(coef_result), rownames(tval_result))
  expect_equal(rownames(tval_result), rownames(pval_result))
})

test_that("safe_extract_matrix preserves numeric types", {
  stats_list <- list(
    gene1 = list(coefs = c(a = 1L, b = 2L), t_vals = c(a = 1.5, b = 2.5), p_vals = c(a = 0.05, b = 0.01))
  )
  all_predictors <- c("a", "b")
  
  result <- safe_extract_matrix(stats_list, "coefs", all_predictors)
  
  expect_true(is.numeric(result))
})

test_that("safe_extract_matrix handles different predictor orders", {
  stats_list <- list(
    gene1 = list(coefs = c(b = 2, a = 1, c = 3)),
    gene2 = list(coefs = c(c = 6, a = 4, b = 5))
  )
  all_predictors <- c("a", "b", "c", "d")
  
  result <- safe_extract_matrix(stats_list, "coefs", all_predictors)
  
  expect_equal(colnames(result), all_predictors)
  expect_true(is.na(result["gene1", "d"]))
  expect_true(is.na(result["gene2", "d"]))
})


# Edge case tests
test_that("Integration: complete workflow with mixed data", {
  stats_list <- list(
    gene1 = list(
      coefs = c(intercept = 0.5, treatment = 1.2),
      t_vals = c(intercept = 2.1, treatment = 3.5),
      p_vals = c(intercept = 0.04, treatment = 0.001)
    ),
    gene2 = NULL,
    gene3 = list(
      coefs = c(treatment = 2.5),
      t_vals = c(treatment = 4.2),
      p_vals = c(treatment = 0.0001)
    ),
    gene4 = list(
      coefs = c(intercept = 1.1, treatment = 0.8, age = 0.3),
      t_vals = NULL,
      p_vals = c(intercept = 0.3, treatment = 0.5, age = 0.7)
    )
  )
  
  all_predictors <- c("intercept", "treatment", "age")
  
  coef <- safe_extract_matrix(stats_list, "coefs", all_predictors)
  t_val <- safe_extract_matrix(stats_list, "t_vals", all_predictors)
  p_val <- safe_extract_matrix(stats_list, "p_vals", all_predictors)
  
  # Check dimensions
  expect_equal(dim(coef), c(4, 3))
  expect_equal(dim(t_val), c(4, 3))
  expect_equal(dim(p_val), c(4, 3))
  
  # Check row names
  expect_equal(rownames(coef), c("gene1", "gene2", "gene3", "gene4"))
  
  # Check specific values
  expect_equal(coef["gene1", "treatment"], 1.2)
  expect_true(all(is.na(coef["gene2", ])))
  expect_true(is.na(coef["gene3", "intercept"]))
  expect_equal(coef["gene4", "age"], 0.3)
  expect_true(all(is.na(t_val["gene4", ])))
})
