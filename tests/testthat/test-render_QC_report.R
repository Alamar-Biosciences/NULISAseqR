test_that("render_QC_report fails early with a clear message on an empty/missing template path (#691)", {
  # Empty string: the failure mode seen in production (system.file() -> "")
  expect_error(
    render_QC_report(xml_files = character(0), Rmd_input_file = ""),
    "QC report template not found"
  )

  # Non-existent path
  expect_error(
    render_QC_report(xml_files = character(0),
                     Rmd_input_file = file.path(tempdir(), "does-not-exist.Rmd")),
    "QC report template not found"
  )

  # NULL / zero-length path
  expect_error(
    render_QC_report(xml_files = character(0), Rmd_input_file = character(0)),
    "QC report template not found"
  )

  # NA path
  expect_error(
    render_QC_report(xml_files = character(0), Rmd_input_file = NA_character_),
    "QC report template not found"
  )

  # length > 1 vector
  expect_error(
    render_QC_report(xml_files = character(0),
                     Rmd_input_file = c("a.Rmd", "b.Rmd")),
    "QC report template not found"
  )

  # a directory (exists, but is not a file) must be rejected
  expect_error(
    render_QC_report(xml_files = character(0), Rmd_input_file = tempdir()),
    "QC report template not found"
  )

  # non-character input must not crash inside nzchar()/file_test() with an
  # opaque error; it should hit the same friendly guard message (#691 review)
  expect_error(
    render_QC_report(xml_files = character(0), Rmd_input_file = factor("a.Rmd")),
    "QC report template not found"
  )
})

test_that("render_QC_report reports NA_character_ template path as NA, not '' (#691 review)", {
  expect_error(
    render_QC_report(xml_files = character(0), Rmd_input_file = NA_character_),
    "Rmd_input_file='NA'",
    fixed = TRUE
  )
})

test_that("render_QC_report rejects an unreadable template file (#691 review)", {
  skip_on_os("windows") # file.access() read-permission semantics are unreliable on Windows
  skip_if(identical(Sys.getenv("USER"), "root") || .Platform$OS.type != "unix",
          "permission bits are not enforced for root")

  unreadable_file <- tempfile(fileext = ".Rmd")
  writeLines("test", unreadable_file)
  Sys.chmod(unreadable_file, mode = "0200") # write-only, not readable
  on.exit(unlink(unreadable_file), add = TRUE)

  if (unname(file.access(unreadable_file, mode = 4)) != 0) {
    expect_error(
      render_QC_report(xml_files = character(0), Rmd_input_file = unreadable_file),
      "QC report template not found"
    )
  } else {
    skip("could not make file unreadable in this environment")
  }
})
