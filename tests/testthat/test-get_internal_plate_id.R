# Tests for get_internal_plate_id() function
# Note: get_internal_plate_id is an internal function, accessed via :::

# =============================================================================
# Tests for get_internal_plate_id()
# =============================================================================

test_that("get_internal_plate_id extracts AUTO_PLATE ID from valid XML", {
  # Create a temporary XML file with AUTO_PLATE attribute
  xml_content <- '<?xml version="1.0" encoding="UTF-8"?>
<NULISAseq AUTO_PLATE="TEST_PLATE_12345" version="1.0">
  <PlateInfo>
    <Sample>Test Sample</Sample>
  </PlateInfo>
</NULISAseq>'

  temp_file <- tempfile(fileext = ".xml")
  writeLines(xml_content, temp_file)

  result <- NULISAseqR:::get_internal_plate_id(temp_file)

  expect_equal(result, "TEST_PLATE_12345")

  # Cleanup
  unlink(temp_file)
})

test_that("get_internal_plate_id handles AUTO_PLATE with special characters", {
  # Test with hyphens, underscores, and alphanumeric characters
  xml_content <- '<?xml version="1.0" encoding="UTF-8"?>
<NULISAseq AUTO_PLATE="PLATE-2024_A1-B2" version="1.0">
  <Data>Test</Data>
</NULISAseq>'

  temp_file <- tempfile(fileext = ".xml")
  writeLines(xml_content, temp_file)

  result <- NULISAseqR:::get_internal_plate_id(temp_file)

  expect_equal(result, "PLATE-2024_A1-B2")

  unlink(temp_file)
})

test_that("get_internal_plate_id handles AUTO_PLATE with spaces", {
  # Test with spaces in the plate ID
  xml_content <- '<?xml version="1.0" encoding="UTF-8"?>
<NULISAseq AUTO_PLATE="PLATE 001 TEST" version="1.0">
  <Data>Test</Data>
</NULISAseq>'

  temp_file <- tempfile(fileext = ".xml")
  writeLines(xml_content, temp_file)

  result <- NULISAseqR:::get_internal_plate_id(temp_file)

  expect_equal(result, "PLATE 001 TEST")

  unlink(temp_file)
})

test_that("get_internal_plate_id returns NA when AUTO_PLATE is missing", {
  # Create XML without AUTO_PLATE attribute
  xml_content <- '<?xml version="1.0" encoding="UTF-8"?>
<NULISAseq version="1.0">
  <PlateInfo>
    <Sample>Test Sample</Sample>
  </PlateInfo>
</NULISAseq>'

  temp_file <- tempfile(fileext = ".xml")
  writeLines(xml_content, temp_file)

  result <- NULISAseqR:::get_internal_plate_id(temp_file)

  expect_true(is.na(result))

  unlink(temp_file)
})

test_that("get_internal_plate_id returns NA for empty AUTO_PLATE", {
  # Test with empty AUTO_PLATE attribute
  xml_content <- '<?xml version="1.0" encoding="UTF-8"?>
<NULISAseq AUTO_PLATE="" version="1.0">
  <Data>Test</Data>
</NULISAseq>'

  temp_file <- tempfile(fileext = ".xml")
  writeLines(xml_content, temp_file)

  result <- NULISAseqR:::get_internal_plate_id(temp_file)

  # Empty AUTO_PLATE won't match the regex pattern [^"]+ (requires at least one char)
  expect_true(is.na(result))

  unlink(temp_file)
})

test_that("get_internal_plate_id handles AUTO_PLATE appearing later in file", {
  # Create XML with AUTO_PLATE not in the first few lines
  xml_content <- paste0(
    '<?xml version="1.0" encoding="UTF-8"?>\n',
    paste(rep('<Comment>Filler line</Comment>\n', 50), collapse = ''),
    '<NULISAseq AUTO_PLATE="LATE_PLATE_ID" version="1.0">\n',
    '  <Data>Test</Data>\n',
    '</NULISAseq>'
  )

  temp_file <- tempfile(fileext = ".xml")
  writeLines(xml_content, temp_file)

  result <- NULISAseqR:::get_internal_plate_id(temp_file)

  expect_equal(result, "LATE_PLATE_ID")

  unlink(temp_file)
})

test_that("get_internal_plate_id handles AUTO_PLATE beyond 1000 lines", {
  # Create XML with AUTO_PLATE after 1000 lines (should return NA due to n=1000 limit)
  xml_content <- paste0(
    '<?xml version="1.0" encoding="UTF-8"?>\n',
    paste(rep('<Comment>Filler line</Comment>\n', 1100), collapse = ''),
    '<NULISAseq AUTO_PLATE="VERY_LATE_PLATE" version="1.0">\n',
    '  <Data>Test</Data>\n',
    '</NULISAseq>'
  )

  temp_file <- tempfile(fileext = ".xml")
  writeLines(xml_content, temp_file)

  result <- NULISAseqR:::get_internal_plate_id(temp_file)

  # Should return NA because AUTO_PLATE is beyond line 1000
  expect_true(is.na(result))

  unlink(temp_file)
})

test_that("get_internal_plate_id returns NA for non-existent file", {
  # Test with a file path that doesn't exist
  fake_path <- tempfile(fileext = ".xml")

  # The function catches errors internally with tryCatch
  result <- NULISAseqR:::get_internal_plate_id(fake_path)

  # Should return NA when file cannot be read
  expect_true(is.na(result))
})

test_that("get_internal_plate_id handles malformed XML gracefully", {
  # Create a file with malformed XML but valid AUTO_PLATE attribute
  xml_content <- '<NULISAseq AUTO_PLATE="MALFORMED_PLATE"
  <unclosed_tag>
  <Data>Test'

  temp_file <- tempfile(fileext = ".xml")
  writeLines(xml_content, temp_file)

  # Should still extract AUTO_PLATE even if XML is malformed
  result <- NULISAseqR:::get_internal_plate_id(temp_file)

  expect_equal(result, "MALFORMED_PLATE")

  unlink(temp_file)
})

test_that("get_internal_plate_id handles file with only plain text", {
  # Test with a non-XML file
  text_content <- "This is just plain text with no AUTO_PLATE attribute"

  temp_file <- tempfile(fileext = ".txt")
  writeLines(text_content, temp_file)

  result <- NULISAseqR:::get_internal_plate_id(temp_file)

  expect_true(is.na(result))

  unlink(temp_file)
})

test_that("get_internal_plate_id matches first AUTO_PLATE when multiple exist", {
  # Test with multiple AUTO_PLATE attributes (should match the first one)
  xml_content <- '<?xml version="1.0" encoding="UTF-8"?>
<NULISAseq AUTO_PLATE="FIRST_PLATE" version="1.0">
  <SubNode AUTO_PLATE="SECOND_PLATE">
    <Data>Test</Data>
  </SubNode>
</NULISAseq>'

  temp_file <- tempfile(fileext = ".xml")
  writeLines(xml_content, temp_file)

  result <- NULISAseqR:::get_internal_plate_id(temp_file)

  # Should return the first match
  expect_equal(result, "FIRST_PLATE")

  unlink(temp_file)
})

test_that("get_internal_plate_id handles AUTO_PLATE with single quotes", {
  # Test with single quotes (should not match, as regex looks for double quotes)
  xml_content <- "<?xml version='1.0' encoding='UTF-8'?>
<NULISAseq AUTO_PLATE='SINGLE_QUOTE_PLATE' version='1.0'>
  <Data>Test</Data>
</NULISAseq>"

  temp_file <- tempfile(fileext = ".xml")
  writeLines(xml_content, temp_file)

  result <- NULISAseqR:::get_internal_plate_id(temp_file)

  # Regex specifically looks for double quotes, so should return NA
  expect_true(is.na(result))

  unlink(temp_file)
})

test_that("get_internal_plate_id handles AUTO_PLATE with escaped characters", {
  # Test with XML entities and escaped characters
  xml_content <- '<?xml version="1.0" encoding="UTF-8"?>
<NULISAseq AUTO_PLATE="PLATE&lt;123&gt;" version="1.0">
  <Data>Test</Data>
</NULISAseq>'

  temp_file <- tempfile(fileext = ".xml")
  writeLines(xml_content, temp_file)

  result <- NULISAseqR:::get_internal_plate_id(temp_file)

  # Should capture the raw attribute value as written (with entities)
  expect_equal(result, "PLATE&lt;123&gt;")

  unlink(temp_file)
})

test_that("get_internal_plate_id handles empty file", {
  # Test with completely empty file
  temp_file <- tempfile(fileext = ".xml")
  writeLines("", temp_file)

  result <- NULISAseqR:::get_internal_plate_id(temp_file)

  expect_true(is.na(result))

  unlink(temp_file)
})

test_that("get_internal_plate_id handles very long plate IDs", {
  # Test with a very long plate ID (e.g., 500 characters)
  long_id <- paste(rep("A", 500), collapse = "")
  xml_content <- sprintf('<?xml version="1.0" encoding="UTF-8"?>
<NULISAseq AUTO_PLATE="%s" version="1.0">
  <Data>Test</Data>
</NULISAseq>', long_id)

  temp_file <- tempfile(fileext = ".xml")
  writeLines(xml_content, temp_file)

  result <- NULISAseqR:::get_internal_plate_id(temp_file)

  expect_equal(result, long_id)
  expect_equal(nchar(result), 500)

  unlink(temp_file)
})

test_that("get_internal_plate_id is case-sensitive for attribute name", {
  # Test with lowercase auto_plate (should not match)
  xml_content <- '<?xml version="1.0" encoding="UTF-8"?>
<NULISAseq auto_plate="lowercase_attr" version="1.0">
  <Data>Test</Data>
</NULISAseq>'

  temp_file <- tempfile(fileext = ".xml")
  writeLines(xml_content, temp_file)

  result <- NULISAseqR:::get_internal_plate_id(temp_file)

  # Should return NA because regex looks for uppercase AUTO_PLATE
  expect_true(is.na(result))

  unlink(temp_file)
})

test_that("get_internal_plate_id handles numeric plate IDs", {
  # Test with purely numeric plate ID
  xml_content <- '<?xml version="1.0" encoding="UTF-8"?>
<NULISAseq AUTO_PLATE="123456789" version="1.0">
  <Data>Test</Data>
</NULISAseq>'

  temp_file <- tempfile(fileext = ".xml")
  writeLines(xml_content, temp_file)

  result <- NULISAseqR:::get_internal_plate_id(temp_file)

  expect_equal(result, "123456789")
  expect_type(result, "character")  # Should still return as character

  unlink(temp_file)
})