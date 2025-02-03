# Define a test for QCSampleCriteria
test_that("skeleton.Rmd can be rendered into HTML", {
  test_path <- paste0(testthat::test_path(), "./../inst/rmarkdown/templates/nulisaseq/skeleton/")
  input <- paste0(test_path, "skeleton.Rmd")
  output <- paste0(test_path, "skeleton.html")
  withr::with_tempfile(output, {
    rmarkdown::render(input)  
    expect_true(file.exists(output), output) 
  })
})
#test_that("skeleton.Rmd can be rendered (INF panel with QCS and SN QC criteria)", {
#  test_path <- paste0(testthat::test_path(), "./../inst/rmarkdown/templates/nulisaseq/skeleton/")
#  input <- paste0(test_path, "skeleton.Rmd")
#  output <- paste0(test_path, "skeleton.html")
#  withr::with_tempfile(output, {
#    rmarkdown::render(input, params=list(xmlFiles=c("example_INF_QCS_SN.xml")))  
#    expect_true(file.exists(output), output) 
#  })
#})
#test_that("skeleton.Rmd can be rendered (CNS, panel with QCS and SN QC criteria)", {
#  test_path <- paste0(testthat::test_path(), "./../inst/rmarkdown/templates/nulisaseq/skeleton/")
#  input <- paste0(test_path, "skeleton.Rmd")
#  output <- paste0(test_path, "skeleton.html")
#  withr::with_tempfile(output, {
#    rmarkdown::render(input, params=list(xmlFiles=c("example_CNS_QCS_SN.xml")))  
#    expect_true(file.exists(output), output) 
#  })
#})
