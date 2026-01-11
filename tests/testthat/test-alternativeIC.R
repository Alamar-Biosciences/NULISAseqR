# Define a test processXML
test_that("Test that loadNULISAseq reads in an XML file and can use an alternate IC", {
  
  test_path <- paste0(testthat::test_path(), "./../inst/rmarkdown/templates/nulisaseq/skeleton/")

  # Test older XMLs that do not have target types, but manually allow setting the IC
  input1 <- paste0(test_path, "detectability_P1_Tr03.xml")
  data1 <- loadNULISAseq(input1, IPC=NULL, IC='mCherry', SC=NULL)
  data2 <- loadNULISAseq(input1, IPC=NULL, IC='CCL7', SC=NULL)
  expect_true(data1$IC_normed$normData[1,1] != data2$IC_normed$normData[1,1])
  expect_true(data1$IC == "mCherry")
  expect_true(data2$IC == "CCL7")

  # Test a newer XML that has targetType but no modifiers 
  input3 <- paste0(test_path, "detectability_P1_Tr03_typemCherry.xml")
  data3 <- loadNULISAseq(input3, IPC=NULL, IC=NULL, SC=NULL)
  expect_true(data3$IC == "mCherry")
  expect_true(all(data3$IC_normed$normData == data1$IC_normed$normData))

  input4 <- paste0(test_path, "detectability_P1_Tr03_typeCCL7.xml")
  data4 <- loadNULISAseq(input4, IPC=NULL, IC=NULL, SC=NULL)
  expect_true(data4$IC == "CCL7")
  expect_true(all(data4$IC_normed$normData == data2$IC_normed$normData))

  # Test a newer XML that has targetType and has multiple ICs (geometric mean)
  input5 <- paste0(test_path, "detectability_P1_Tr03_typemCherry_CCL7.xml")
  data5 <- loadNULISAseq(input5, IPC=NULL, IC=NULL, SC=NULL)
  expect_true("mCherry" %in% data5$IC) 
  expect_true("CCL7" %in% data5$IC)
  expect_false(all(data5$IC_normed$normData == data4$IC_normed$normData))
  expect_false(all(data5$IC_normed$normData == data3$IC_normed$normData))
  
  # Test that when multiple ICs exist, that only one can be selected
  data5 <- loadNULISAseq(input5, IPC=NULL, IC="mCherry", SC=NULL)
  expect_true("mCherry" %in% data5$IC) 
  expect_false("CCL7" %in% data5$IC)
  expect_true(all(data5$IC_normed$normData == data3$IC_normed$normData))
  
  data5 <- loadNULISAseq(input5, IPC=NULL, IC="CCL7", SC=NULL)
  expect_false("mCherry" %in% data5$IC) 
  expect_true("CCL7" %in% data5$IC)
  expect_true(all(data5$IC_normed$normData == data4$IC_normed$normData))

  # Test a newer XML that has targetType and uses hide on a control
  # Hidden targets are completely removed from data
  input6 <- paste0(test_path, "detectability_P1_Tr03_typemCherry_CCL7_hide.xml")
  data6 <- loadNULISAseq(input6, IPC=NULL, IC=NULL, SC=NULL)
  expect_true("CCL7" %in% data6$IC)
  expect_false("mCherry" %in% data6$IC)
  expect_false("mCherry" %in% data6$targets$targetName, info = "Hidden target mCherry should be removed from targets")

  # Test a newer XML that has targetType and uses noDetectability
  input7 <- paste0(test_path, "detectability_P1_Tr03_typemCherry_noDetectability_AGER.xml")
  data7 <- loadNULISAseq(input7, IPC=NULL, IC=NULL, SC=NULL)
  expect_false("AGER" %in% data7$detectability$all$detectable)
  expect_false("AGER" %in% data7$detectability$all$detectablity)

  # Test a newer XML that has targetType and uses hide on a control and a target
  # Hidden targets are completely removed from data
  input6 <- paste0(test_path, "detectability_P1_Tr03_typemCherry_CCL7_hide2.xml")
  data6 <- loadNULISAseq(input6, IPC=NULL, IC=NULL, SC=NULL)
  expect_true("CCL7" %in% data6$IC)
  expect_false("mCherry" %in% data6$IC)
  expect_false("mCherry" %in% data6$targets$targetName, info = "Hidden target mCherry should be removed from targets")
  expect_false("WNT7A" %in% data6$targets$targetName, info = "Hidden target WNT7A should be removed from targets")
  expect_false("mCherry" %in% rownames(data6$Data), info = "Hidden target mCherry should be removed from Data matrix")
  expect_false("WNT7A" %in% rownames(data6$Data), info = "Hidden target WNT7A should be removed from Data matrix")
})
