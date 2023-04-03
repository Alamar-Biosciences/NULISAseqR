# Define a test for QCSampleCriteria
test_that("QCSampleCriteria returns expected output", {
  expected_thresholds <- c(Detectability=0.8, ICReads=1000, NumReads=500000, IC_Median="-0.3,0.3")
  expected_operators <- c(Detectability=">", ICReads=">", NumReads=">", IC_Median=">,<")
  
  qc_criteria <- QCSampleCriteria()
  
  expect_equal(qc_criteria$thresholds, expected_thresholds, 
               info = "Expected thresholds do not match actual thresholds")
  expect_equal(qc_criteria$operators, expected_operators, 
               info = "Expected operators do not match actual operators")
})
