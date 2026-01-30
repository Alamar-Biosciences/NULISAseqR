# Tests for importNULISAseq() validation logic

# =============================================================================
# Tests for AUTO_PLATE ID duplicate detection with named parameters
# =============================================================================

test_that("importNULISAseq detects duplicate IDs with named excludeSamples", {
  # Create two temporary XML files with the same AUTO_PLATE ID
  xml_content <- '<?xml version="1.0" encoding="UTF-8"?>
<NULISAseq AUTO_PLATE="DUPLICATE_ID" version="1.0">
  <Data>Test</Data>
</NULISAseq>'

  temp_file1 <- tempfile(fileext = ".xml")
  temp_file2 <- tempfile(fileext = ".xml")
  writeLines(xml_content, temp_file1)
  writeLines(xml_content, temp_file2)

  # Attempt to import with named excludeSamples should fail
  expect_error(
    importNULISAseq(
      files = c(temp_file1, temp_file2),
      excludeSamples = list(DUPLICATE_ID = c("sample1"))
    ),
    "Duplicate AUTO_PLATE ID.*DUPLICATE_ID.*excludeSamples"
  )

  unlink(c(temp_file1, temp_file2))
})

test_that("importNULISAseq detects duplicate IDs with named IPC", {
  xml_content <- '<?xml version="1.0" encoding="UTF-8"?>
<NULISAseq AUTO_PLATE="DUP_PLATE" version="1.0">
  <Data>Test</Data>
</NULISAseq>'

  temp_file1 <- tempfile(fileext = ".xml")
  temp_file2 <- tempfile(fileext = ".xml")
  writeLines(xml_content, temp_file1)
  writeLines(xml_content, temp_file2)

  # Attempt to import with named IPC should fail
  expect_error(
    importNULISAseq(
      files = c(temp_file1, temp_file2),
      IPC = list(DUP_PLATE = c("IPC1", "IPC2"))
    ),
    "Duplicate AUTO_PLATE ID.*DUP_PLATE.*IPC"
  )

  unlink(c(temp_file1, temp_file2))
})

test_that("importNULISAseq detects duplicate IDs with named SC", {
  xml_content <- '<?xml version="1.0" encoding="UTF-8"?>
<NULISAseq AUTO_PLATE="SAME_ID" version="1.0">
  <Data>Test</Data>
</NULISAseq>'

  temp_file1 <- tempfile(fileext = ".xml")
  temp_file2 <- tempfile(fileext = ".xml")
  writeLines(xml_content, temp_file1)
  writeLines(xml_content, temp_file2)

  # Attempt to import with named SC should fail
  expect_error(
    importNULISAseq(
      files = c(temp_file1, temp_file2),
      SC = list(SAME_ID = c("SC1"))
    ),
    "Duplicate AUTO_PLATE ID.*SAME_ID.*SC"
  )

  unlink(c(temp_file1, temp_file2))
})

test_that("importNULISAseq detects duplicate IDs with multiple named parameters", {
  xml_content <- '<?xml version="1.0" encoding="UTF-8"?>
<NULISAseq AUTO_PLATE="MULTI_DUP" version="1.0">
  <Data>Test</Data>
</NULISAseq>'

  temp_file1 <- tempfile(fileext = ".xml")
  temp_file2 <- tempfile(fileext = ".xml")
  writeLines(xml_content, temp_file1)
  writeLines(xml_content, temp_file2)

  # Attempt to import with multiple named parameters should fail
  # Error message should list all named parameters
  expect_error(
    importNULISAseq(
      files = c(temp_file1, temp_file2),
      SC = list(MULTI_DUP = c("SC1")),
      NC = list(MULTI_DUP = c("NC1")),
      excludeTargets = list(MULTI_DUP = c("target1"))
    ),
    "Duplicate AUTO_PLATE ID.*MULTI_DUP"
  )

  # Check that error message mentions the specific parameters
  tryCatch(
    importNULISAseq(
      files = c(temp_file1, temp_file2),
      SC = list(MULTI_DUP = c("SC1")),
      NC = list(MULTI_DUP = c("NC1")),
      excludeTargets = list(MULTI_DUP = c("target1"))
    ),
    error = function(e) {
      expect_match(e$message, "SC")
      expect_match(e$message, "NC")
      expect_match(e$message, "excludeTargets")
    }
  )

  unlink(c(temp_file1, temp_file2))
})

test_that("importNULISAseq allows duplicate IDs with user-provided plateName", {
  xml_content <- '<?xml version="1.0" encoding="UTF-8"?>
<NULISAseq AUTO_PLATE="DUPLICATE" version="1.0">
  <Data>Test</Data>
</NULISAseq>'

  temp_file1 <- tempfile(fileext = ".xml")
  temp_file2 <- tempfile(fileext = ".xml")
  writeLines(xml_content, temp_file1)
  writeLines(xml_content, temp_file2)

  # Should NOT error about duplicate IDs when plateName is provided
  # Will error during XML parsing, but we check that it's NOT a duplicate ID error
  result <- tryCatch(
    importNULISAseq(
      files = c(temp_file1, temp_file2),
      plateName = c("CustomPlate1", "CustomPlate2"),
      excludeSamples = list(CustomPlate1 = c("sample1"))
    ),
    error = function(e) e$message
  )

  # Should NOT contain duplicate ID error message
  expect_false(grepl("Duplicate AUTO_PLATE ID", result))

  unlink(c(temp_file1, temp_file2))
})

test_that("importNULISAseq allows duplicate IDs without named parameters", {
  xml_content <- '<?xml version="1.0" encoding="UTF-8"?>
<NULISAseq AUTO_PLATE="DUPLICATE" version="1.0">
  <Data>Test</Data>
</NULISAseq>'

  temp_file1 <- tempfile(fileext = ".xml")
  temp_file2 <- tempfile(fileext = ".xml")
  writeLines(xml_content, temp_file1)
  writeLines(xml_content, temp_file2)

  # Should NOT error about duplicate IDs when using vector parameters
  result <- tryCatch(
    importNULISAseq(
      files = c(temp_file1, temp_file2),
      excludeSamples = c("sample1")  # Vector, not named list
    ),
    error = function(e) e$message
  )

  # Should NOT contain duplicate ID error message
  expect_false(grepl("Duplicate AUTO_PLATE ID", result))

  unlink(c(temp_file1, temp_file2))
})

test_that("importNULISAseq handles unique AUTO_PLATE IDs with named parameters", {
  xml_content1 <- '<?xml version="1.0" encoding="UTF-8"?>
<NULISAseq AUTO_PLATE="UNIQUE_ID_1" version="1.0">
  <Data>Test</Data>
</NULISAseq>'

  xml_content2 <- '<?xml version="1.0" encoding="UTF-8"?>
<NULISAseq AUTO_PLATE="UNIQUE_ID_2" version="1.0">
  <Data>Test</Data>
</NULISAseq>'

  temp_file1 <- tempfile(fileext = ".xml")
  temp_file2 <- tempfile(fileext = ".xml")
  writeLines(xml_content1, temp_file1)
  writeLines(xml_content2, temp_file2)

  # Should NOT error about duplicate IDs with unique IDs and named parameters
  result <- tryCatch(
    importNULISAseq(
      files = c(temp_file1, temp_file2),
      excludeSamples = list(UNIQUE_ID_1 = c("sample1"),
                            UNIQUE_ID_2 = c("sample2"))
    ),
    error = function(e) e$message
  )

  # Should NOT contain duplicate ID error message
  expect_false(grepl("Duplicate AUTO_PLATE ID", result))

  unlink(c(temp_file1, temp_file2))
})

test_that("importNULISAseq error message includes all named parameters", {
  xml_content <- '<?xml version="1.0" encoding="UTF-8"?>
<NULISAseq AUTO_PLATE="TEST_DUP" version="1.0">
  <Data>Test</Data>
</NULISAseq>'

  temp_file1 <- tempfile(fileext = ".xml")
  temp_file2 <- tempfile(fileext = ".xml")
  writeLines(xml_content, temp_file1)
  writeLines(xml_content, temp_file2)

  # Test that error message is informative
  expect_error(
    importNULISAseq(
      files = c(temp_file1, temp_file2),
      Bridge = list(TEST_DUP = c("Bridge1")),
      Calibrator = list(TEST_DUP = c("Cal1"))
    ),
    "manually define 'plateName'"
  )

  unlink(c(temp_file1, temp_file2))
})

test_that("importNULISAseq handles missing AUTO_PLATE IDs gracefully", {
  # XML without AUTO_PLATE attribute
  xml_content <- '<?xml version="1.0" encoding="UTF-8"?>
<NULISAseq version="1.0">
  <Data>Test</Data>
</NULISAseq>'

  temp_file1 <- tempfile(fileext = ".xml")
  temp_file2 <- tempfile(fileext = ".xml")
  writeLines(xml_content, temp_file1)
  writeLines(xml_content, temp_file2)

  # Should fall back to Plate_XX naming and not error at validation stage
  result <- tryCatch(
    importNULISAseq(
      files = c(temp_file1, temp_file2),
      excludeSamples = list(Plate_01 = c("sample1"))
    ),
    error = function(e) e$message
  )

  # Should NOT contain duplicate ID error message
  expect_false(grepl("Duplicate AUTO_PLATE ID", result))

  unlink(c(temp_file1, temp_file2))
})

test_that("importNULISAseq handles multiple NA AUTO_PLATE IDs correctly", {
  # Multiple files without AUTO_PLATE attributes (all return NA)
  # Should NOT incorrectly report NA as a duplicate ID
  xml_content <- '<?xml version="1.0" encoding="UTF-8"?>
<NULISAseq version="1.0">
  <Data>Test</Data>
</NULISAseq>'

  temp_file1 <- tempfile(fileext = ".xml")
  temp_file2 <- tempfile(fileext = ".xml")
  temp_file3 <- tempfile(fileext = ".xml")
  writeLines(xml_content, temp_file1)
  writeLines(xml_content, temp_file2)
  writeLines(xml_content, temp_file3)

  # With multiple NA IDs and named parameters, should fall back to Plate_XX
  # Should NOT error about "Duplicate AUTO_PLATE ID 'NA'"
  result <- tryCatch(
    importNULISAseq(
      files = c(temp_file1, temp_file2, temp_file3),
      excludeSamples = list(Plate_01 = c("sample1"), Plate_02 = c("sample2"))
    ),
    error = function(e) e$message
  )

  # Should NOT contain duplicate ID error message
  expect_false(grepl("Duplicate AUTO_PLATE ID", result))
  # Should definitely NOT mention NA as duplicate
  expect_false(grepl("Duplicate AUTO_PLATE ID.*NA", result))

  unlink(c(temp_file1, temp_file2, temp_file3))
})

test_that("importNULISAseq detects actual duplicates with some NAs present", {
  # Mix of files: some with duplicate IDs, some without AUTO_PLATE (NA)
  xml_with_id <- '<?xml version="1.0" encoding="UTF-8"?>
<NULISAseq AUTO_PLATE="DUP_ID" version="1.0">
  <Data>Test</Data>
</NULISAseq>'

  xml_without_id <- '<?xml version="1.0" encoding="UTF-8"?>
<NULISAseq version="1.0">
  <Data>Test</Data>
</NULISAseq>'

  temp_file1 <- tempfile(fileext = ".xml")
  temp_file2 <- tempfile(fileext = ".xml")
  temp_file3 <- tempfile(fileext = ".xml")
  writeLines(xml_with_id, temp_file1)
  writeLines(xml_with_id, temp_file2)  # Duplicate of file1
  writeLines(xml_without_id, temp_file3)  # No AUTO_PLATE (NA)

  # Should detect the duplicate "DUP_ID" and NOT mention NA
  expect_error(
    importNULISAseq(
      files = c(temp_file1, temp_file2, temp_file3),
      excludeSamples = list(DUP_ID = c("sample1"))
    ),
    "Duplicate AUTO_PLATE ID.*DUP_ID"
  )

  # Verify error message does NOT contain NA
  tryCatch(
    importNULISAseq(
      files = c(temp_file1, temp_file2, temp_file3),
      excludeSamples = list(DUP_ID = c("sample1"))
    ),
    error = function(e) {
      expect_false(grepl("\\bNA\\b", e$message))
    }
  )

  unlink(c(temp_file1, temp_file2, temp_file3))
})