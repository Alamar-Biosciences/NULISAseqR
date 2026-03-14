# helper function for reading xlsx files
read_all_sheets <- function(path, skip_by_sheet = list()) {
  sheets <- readxl::excel_sheets(path)
  out <- lapply(sheets, function(s) {
    skip <- skip_by_sheet[[s]] %||% 0
    readxl::read_excel(path, sheet = s, skip = skip, col_types = "text")
  })
  names(out) <- sheets
  out
}

compare_xlsx <- function(actual_path, expected_path, skip_by_sheet = list(),
                         numeric_tolerance = 1e-4) {
  actual   <- read_all_sheets(actual_path, skip_by_sheet)
  expected <- read_all_sheets(expected_path, skip_by_sheet)
  
  expect_equal(names(actual), names(expected),
               label = "Sheet names should match")
  
  for (sheet in names(expected)) {
    act <- actual[[sheet]]
    exp <- expected[[sheet]]
    
    expect_equal(colnames(act), colnames(exp),
                 label = paste("Sheet:", sheet, "| Column names"))
    expect_equal(nrow(act), nrow(exp),
                 label = paste("Sheet:", sheet, "| Row count"))
    
    for (col in colnames(exp)) {
      act_num  <- suppressWarnings(as.numeric(act[[col]]))
      exp_num  <- suppressWarnings(as.numeric(exp[[col]]))
      is_numeric <- !is.na(act_num) & !is.na(exp_num)
      
      if (any(is_numeric)) {
        # compare numeric values with tolerance
        expect_equal(act_num[is_numeric], exp_num[is_numeric],
                     tolerance = numeric_tolerance,
                     label = paste("Sheet:", sheet, "| Column:", col, "(numeric)"))
        # compare non-numeric values exactly
        expect_equal(act[[col]][!is_numeric], exp[[col]][!is_numeric],
                     label = paste("Sheet:", sheet, "| Column:", col, "(character)"))
      } else {
        expect_equal(act[[col]], exp[[col]],
                     label = paste("Sheet:", sheet, "| Column:", col))
      }
    }
  }
}


test_that("writeNULISAseq RQ output matches expected", {
  xml_files <- c("20260129_Neuro220_RRD_CSF.xml",
                 "20260224_Neuro220_RRD_PlasmaSerum.xml")
  
  temp_dir <- tempdir()
  for (f in xml_files) {
    file.copy(test_path("fixtures", f), file.path(temp_dir, f), overwrite = TRUE)
  }
  output_file <- tempfile(tmpdir = temp_dir, fileext = ".xlsx")
  
  writeNULISAseq(xml_files = xml_files,
                 dataDir = temp_dir,
                 output_filename = basename(output_file),
                 PanelLotNumber = "panelLot040",
                 metadata = list(NULISAseqR = 'dlfjqp13948eo',
                                 NULISAseqAQ = 'e98249lkfjqke'),
                 include_unnorm_counts = TRUE,
                 include_IPC = TRUE,
                 include_SC = TRUE,
                 include_NC = TRUE,
                 include_IC_counts = TRUE,
                 output_TAP_AQ = FALSE)
  
  compare_xlsx(
    actual_path   = output_file,
    expected_path = test_path("fixtures", "writeNULISAseq_RQ_Neuro220_RRD.xlsx")
  )
  
  unlink(c(file.path(temp_dir, xml_files), output_file))
})

test_that("writeNULISAseq AQ output matches expected", {
  skip_if_not(requireNamespace("NULISAseqAQ", quietly = TRUE),
              "NULISAseqAQ package not available")
  
  xml_file <- "Analysis_INF250_Lot4_AQ_LC_R3_20241229.xml"
  
  temp_dir <- tempdir()
  file.copy(test_path("fixtures", xml_file), file.path(temp_dir, xml_file), overwrite = TRUE)
  output_file <- tempfile(tmpdir = temp_dir, fileext = ".xlsx")
  
  writeNULISAseq(xml_files = xml_file,
                 dataDir = temp_dir,
                 output_filename = basename(output_file),
                 PanelLotNumber = "panelLot040",
                 metadata = list(NULISAseqR = 'dlfjqp13948eo',
                                 NULISAseqAQ = 'e98249lkfjqke'),
                 include_unnorm_counts = TRUE,
                 include_IPC = TRUE,
                 include_SC = TRUE,
                 include_NC = TRUE,
                 include_IC_counts = TRUE,
                 output_TAP_AQ = TRUE)
  
  compare_xlsx(
    actual_path   = output_file,
    expected_path = test_path("fixtures", "writeNULISAseq_AQ_Analysis_INF250_Lot4_AQ_LC_R3_20241229.xlsx"),
    skip_by_sheet = list("AQ (pg per mL and aM)" = 2)
  )
  
  unlink(c(file.path(temp_dir, xml_file), output_file))
})
