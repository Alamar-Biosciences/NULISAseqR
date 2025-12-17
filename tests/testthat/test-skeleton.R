# Define a test for QCSampleCriteria
test_that("skeleton.Rmd can be rendered into HTML", {
  test_path <- paste0(testthat::test_path(), "./../inst/rmarkdown/templates/nulisaseq/skeleton/")
  input <- paste0(test_path, "skeleton.Rmd")
  files <- c("Analysis_INF250_Lot4_AQ_LC_R3_20241229.xml", "Analysis_LC_R3_INF250_Lot3_20240727.xml")
  output <- paste0(test_path, "skeleton.html")
  withr::with_tempfile(output, {
    rmarkdown::render(input, params=list(dataDir=".", xmlFiles=files))  
    expect_true(file.exists(output), output) 
  })
})
test_that("skeleton.Rmd can be rendered into HTML for XMLs with differing targets", {
  test_path <- paste0(testthat::test_path(), "./../inst/rmarkdown/templates/nulisaseq/skeleton/")
  input <- paste0(test_path, "skeleton.Rmd")
  output <- paste0(test_path, "skeleton.html")
  files <- c("Analysis_INF250_Lot4_AQ_LC_R3_20241229.xml", "Analysis_LC_R3_INF250_Lot3_20240727.xml")
  withr::with_tempfile(output, {
    rmarkdown::render(input, params=list(dataDir=".", xmlFiles=files))  
    expect_true(file.exists(output), output)
  })
})
test_that("skeleton.Rmd can be rendered into HTML with an alternative IC", {
  test_path <- paste0(testthat::test_path(), "./../inst/rmarkdown/templates/nulisaseq/skeleton/")
  input <- paste0(test_path, "skeleton.Rmd")
  output <- paste0(test_path, "skeleton.html")
  files <- c("detectability_P1_Tr03_typemCherry_CCL7_hide2.xml")
  withr::with_tempfile(output, {
    rmarkdown::render(input, params=list(dataDir=".", xmlFiles=files))  
    expect_true(file.exists(output), output)
    orig_content <- readLines(output, warn=FALSE)
    patterns <- c("-mCherry-summary", files)
    html_content <- gsub(paste(patterns, collapse = "|"), "", orig_content)
    expect_true(any(grepl("CCL7", html_content)), info="Output does not contain 'CCL7'")
    expect_false(any(grepl("mCherry", html_content)), info="Output does contain 'mCherry'")
    expect_false(any(grepl("WNT7A", html_content)), info="Output does contain 'WNT7A'")
  })
})
test_that("skeleton.Rmd can be rendered into HTML, NAS version", { 
  test_path <- paste0(testthat::test_path(), "./../inst/rmarkdown/templates/nulisaseq/skeleton/")
  input <- paste0(test_path, "skeleton.Rmd")
  files <- c("Analysis_INF250_Lot4_AQ_LC_R3_20241229.xml", "Analysis_LC_R3_INF250_Lot3_20240727.xml")
  output <- paste0(test_path, "skeleton.html")
  withr::with_tempfile(output, {
    rmarkdown::render(input, params=list(reportType="webApp", outPlateEffect=FALSE, dataDir=".", xmlFiles=files ))  
    expect_true(file.exists(output), output) 
  })
})
test_that("skeleton.Rmd can be rendered into HTML, NAS version, advancedQC", { 
  test_path <- paste0(testthat::test_path(), "./../inst/rmarkdown/templates/nulisaseq/skeleton/")
  input <- paste0(test_path, "skeleton.Rmd")
  files <- c("Analysis_INF250_Lot4_AQ_LC_R3_20241229.xml", "Analysis_LC_R3_INF250_Lot3_20240727.xml")
  output <- paste0(test_path, "skeleton.html")
  withr::with_tempfile(output, {
    rmarkdown::render(input, params=list(reportType="webApp", outPlateEffect=TRUE, advancedQC=TRUE, dataDir=".", xmlFiles=files ))  
    expect_true(file.exists(output), output) 
  })
})
test_that("skeleton.Rmd can be rendered into HTML, NAS version, XML_v1.3.0", { 
  test_path <- paste0(testthat::test_path(), "./../inst/rmarkdown/templates/nulisaseq/skeleton/")
  input <- paste0(test_path, "skeleton.Rmd")
  files <- c("XML_v1.3.0.xml")
  output <- paste0(test_path, "skeleton.html")
  withr::with_tempfile(output, {
    rmarkdown::render(input, params=list(reportType="webApp", outPlateEffect=TRUE, advancedQC=FALSE, dataDir=".", xmlFiles=files ))  
    expect_true(file.exists(output), output) 
  })
})
