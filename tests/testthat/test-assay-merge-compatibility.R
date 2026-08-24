# Tests for the CNS/Inflammation panel-variant merge compatibility logic
# (issue #3263 / NAS PR #3621). Covers assay_merge_key(), canonical_assay_label(),
# assay_identity() and check_assay_type().
#
# All five are unexported internals, called via NULISAseqR::: throughout - matching
# this package's own convention (see test-detectability-format.R) and required for
# CI, which runs `library(NULISAseqR); testthat::test_dir(...)` against the
# installed package rather than devtools::test()/test_check(), so bare names are
# not visible.

mk_execution_details <- function(assays) {
  lapply(assays, function(a) list(Assay = a))
}

test_that("assay_merge_key strips version qualifiers but keeps panel identity", {
  expect_equal(
    NULISAseqR:::assay_merge_key(c("CNS Disease Panel 120 (V2)", "CNS Disease Panel 120 (V3)")),
    c("cns disease panel", "cns disease panel")
  )
  # Bare dotted version (v1.9) and parenthesised major version (V2) reduce to the same key.
  expect_equal(
    NULISAseqR:::assay_merge_key(c("CNS Disease Panel v1.9", "CNS Disease Panel 120 (V2)")),
    c("cns disease panel", "cns disease panel")
  )
  # DBS protocol suffix does not affect the key.
  expect_equal(
    NULISAseqR:::assay_merge_key(c("CNS Disease Panel 120 (V2)", "CNS Disease Panel 120 (V2) (DBS Compatible)")),
    c("cns disease panel", "cns disease panel")
  )
  # NULISAseq platform prefix is not panel identity.
  expect_equal(
    NULISAseqR:::assay_merge_key(c("NULISAseq Inflammation Panel 250", "Inflammation Panel v1.8")),
    c("inflammation panel", "inflammation panel")
  )
})

test_that("assay_merge_key keeps plex count as a discriminator for families not in REDUNDANT_ASSAY_PLEX", {
  # CNS Disease Panel and Inflammation Panel have a known single plex, so the number is dropped.
  # A family with no entry (e.g. Immune) keeps its plex number as part of the key.
  expect_equal(
    NULISAseqR:::assay_merge_key(c("Immune 340 Panel (V2)", "Immune 340 Panel (V3)")),
    c("immune 340 panel", "immune 340 panel")
  )
  expect_false(
    NULISAseqR:::assay_merge_key("Immune 340 Panel (V2)") == NULISAseqR:::assay_merge_key("Immune 250 Panel (V2)")
  )
})

test_that("assay_merge_key keeps panel family and AQ/RQ as hard discriminators", {
  expect_false(
    NULISAseqR:::assay_merge_key("CNS Disease Panel 120 (V2)") == NULISAseqR:::assay_merge_key("Inflammation Panel 250 (V2)")
  )
  expect_false(
    NULISAseqR:::assay_merge_key("Inflammation Panel 250 - AQ") == NULISAseqR:::assay_merge_key("Inflammation Panel 250")
  )
})

test_that("canonical_assay_label picks the longest de-qualified name", {
  expect_equal(
    NULISAseqR:::canonical_assay_label(c("CNS Disease Panel 120 (V2)", "CNS Disease Panel v1.9")),
    "CNS Disease Panel 120"
  )
  expect_equal(
    NULISAseqR:::canonical_assay_label(c("NULISAseq Inflammation Panel 250", "Inflammation Panel v1.8")),
    "Inflammation Panel 250"
  )
  # A single name is returned untouched (existing single-run badges must not shift).
  expect_equal(
    NULISAseqR:::canonical_assay_label("CNS Disease Panel 120 (V2)"),
    "CNS Disease Panel 120 (V2)"
  )
  # No usable names at all.
  expect_null(NULISAseqR:::canonical_assay_label(NA_character_))
})

test_that("the 'NULISAseq' RQ default is a real assay name, not a stand-in for missing", {
  # readNULISAseq.R defaults an RQ run's Assay to the literal string "NULISAseq" when its XML has
  # no <Assay> node. That default is an intentional, real value - it is compared like any other
  # assay name, not normalised away as unknown. (Genuinely absent Assay - e.g. NULL on an AQ run -
  # is the separate, already-established "unknown" case covered below.)

  # Two runs that both hit the default carry the identical (generic) name: compatible with
  # each other, same as any pair of identically-named runs.
  ed_both_default <- mk_execution_details(c("NULISAseq", "NULISAseq"))
  expect_true(NULISAseqR:::check_assay_type(ed_both_default))

  # A run with the default and a run with a real, different panel name are not established to be
  # the same panel - the default carries no panel information, so it must not pass silently.
  ed_mixed <- mk_execution_details(c("NULISAseq", "CNS Disease Panel 120 (V2)"))
  expect_error(NULISAseqR:::check_assay_type(ed_mixed))

  # Both runs report the identical string, so that string is the (only sensible) label - the
  # single-name short-circuit in canonical_assay_label() applies here exactly as it would for two
  # runs that happen to share any other real name.
  identity_both_default <- NULISAseqR:::assay_identity(ed_both_default)
  expect_equal(identity_both_default$label, "NULISAseq")
  expect_null(identity_both_default$variants)

  # A real panel name that merely starts with "NULISAseq " (the vendor prefix, not the bare
  # default) is unaffected - only the bare token is stripped by NON_IDENTIFYING_ASSAY_TOKENS.
  ed_real <- mk_execution_details(c("NULISAseq Inflammation Panel 250", "Inflammation Panel v1.8"))
  expect_true(NULISAseqR:::check_assay_type(ed_real))
  expect_equal(NULISAseqR:::assay_identity(ed_real)$label, "Inflammation Panel 250")
})

test_that("a genuinely missing Assay (NULL, e.g. an AQ run) also normalises to 'NULISAseq'", {
  # extract_assay_names() applies the same "no assay reported -> NULISAseq" default uniformly at
  # this package's own boundary with loadNULISAseq()'s output, regardless of why the value is
  # absent (NULL here, vs. readNULISAseq.R's own RQ-only default). There is no separate "unknown,
  # compatible with anything" bucket - NULL becomes NULISAseq and follows the same name rules as
  # any other value.

  # Two NULL runs both become "NULISAseq": identical, compatible.
  ed_both_null <- mk_execution_details(list(NULL, NULL))
  expect_true(NULISAseqR:::check_assay_type(ed_both_null))
  identity_both_null <- NULISAseqR:::assay_identity(ed_both_null)
  expect_equal(identity_both_null$label, "NULISAseq")
  expect_null(identity_both_null$variants)

  # A NULL run and a run with a real, different panel name are not established to be the same
  # panel - blocked, same as the explicit "NULISAseq" placeholder case above.
  ed_null_and_real <- mk_execution_details(list(NULL, "CNS Disease Panel 120 (V2)"))
  expect_error(NULISAseqR:::check_assay_type(ed_null_and_real))

  # A whitespace-only Assay (readNULISAseq() trims <Assay> but only converts zero-length values
  # to NULL, so "   " arrives here as "") normalises to "NULISAseq" like every other absent form -
  # it is never passed through as "" nor recorded as a variant.
  ed_blank <- mk_execution_details(list("   ", NULL))
  expect_equal(NULISAseqR:::extract_assay_names(ed_blank), c("NULISAseq", "NULISAseq"))
  identity_blank <- NULISAseqR:::assay_identity(ed_blank)
  expect_equal(identity_blank$label, "NULISAseq")
  expect_null(identity_blank$variants)
  expect_error(NULISAseqR:::check_assay_type(mk_execution_details(list("  ", "CNS Disease Panel 120 (V2)"))))
})

test_that("assay_identity reports variants only when more than one distinct name merged", {
  single <- NULISAseqR:::assay_identity(mk_execution_details("CNS Disease Panel 120 (V2)"))
  expect_equal(single$label, "CNS Disease Panel 120 (V2)")
  expect_null(single$variants)

  merged <- NULISAseqR:::assay_identity(mk_execution_details(
    c("CNS Disease Panel 120 (V2)", "CNS Disease Panel v1.9")
  ))
  expect_equal(merged$label, "CNS Disease Panel 120")
  expect_setequal(merged$variants, c("CNS Disease Panel 120 (V2)", "CNS Disease Panel v1.9"))

  # All Assay fields missing: both normalise to "NULISAseq", which becomes the label (identical to
  # any other pair of identically-named runs) - not NULL, and not recorded as a variant either
  # (there is only one distinct name, not "more than one merged").
  none <- NULISAseqR:::assay_identity(mk_execution_details(list(NULL, NULL)))
  expect_equal(none$label, "NULISAseq")
  expect_null(none$variants)
})

test_that("check_assay_type: compatibility truth table", {
  compatible <- list(
    "identical names"                     = c("CNS Disease Panel 120 (V2)", "CNS Disease Panel 120 (V2)"),
    "paren-version variants (V2 + V3)"    = c("CNS Disease Panel 120 (V2)", "CNS Disease Panel 120 (V3)"),
    "paren-version + DBS protocol"        = c("CNS Disease Panel 120 (V2)", "CNS Disease Panel 120 (V2) (DBS Compatible)"),
    "bare dotted version + paren version" = c("CNS Disease Panel v1.9", "CNS Disease Panel 120 (V2)"),
    "NULISAseq prefix vs no prefix"       = c("NULISAseq Inflammation Panel 250", "Inflammation Panel v1.8"),
    "Inflammation V2 + V3"                = c("Inflammation Panel 250 (V2)", "Inflammation Panel 250 (V3)"),
    "family with no redundant-plex entry" = c("Immune 340 Panel (V2)", "Immune 340 Panel (V3)"),
    "all plates missing Assay (both -> NULISAseq)" = c(NA, NA),
    "three-way, all same merge key"       = c("CNS Disease Panel 120 (V2)", "CNS Disease Panel v1.9", "CNS Disease Panel 120 (V3)")
  )
  for (name in names(compatible)) {
    ed <- mk_execution_details(compatible[[name]])
    expect_true(NULISAseqR:::check_assay_type(ed), info = name)
  }

  incompatible <- list(
    "same family, different plex"       = c("CNS Disease Panel 120 (V2)", "CNS Disease Panel 250 (V2)"),
    "different panel families"          = c("CNS Disease Panel 120 (V2)", "Inflammation Panel 250 (V2)"),
    "AQ vs RQ, same family"              = c("Inflammation Panel 250 - AQ", "Inflammation Panel 250"),
    "different plex, family not in REDUNDANT_ASSAY_PLEX" = c("Immune 340 Panel (V2)", "Immune 250 Panel (V2)"),
    "three-way, one incompatible"        = c("CNS Disease Panel 120 (V2)", "CNS Disease Panel v1.9", "Inflammation Panel 250"),
    "one plate missing Assay (-> NULISAseq) vs a real name" = c(NA, "CNS Disease Panel 120 (V2)")
  )
  for (name in names(incompatible)) {
    ed <- mk_execution_details(incompatible[[name]])
    expect_error(NULISAseqR:::check_assay_type(ed), info = name)
  }
})

test_that("check_assay_type reports incompatible assays with file names in the error message", {
  ed <- mk_execution_details(c("CNS Disease Panel 120 (V2)", "Inflammation Panel 250 (V2)"))
  expect_error(
    NULISAseqR:::check_assay_type(ed, fileNames = c("plate_a.xml", "plate_b.xml")),
    "plate_a\\.xml.*plate_b\\.xml"
  )
})

test_that("combine_targets retains non-shared targets with NA rather than dropping them", {
  dataList <- list(
    list(targets = data.frame(targetName = c("T1", "T2"), targetType = "protein")),
    list(targets = data.frame(targetName = c("T1", "T2", "T3"), targetType = "protein"))
  )
  res <- NULISAseqR:::combine_targets(dataList, plateID = c("p1", "p2"))
  expect_equal(res$excluded, "T3")
  # The excluded target stays in the merged targets frame - it is not filtered out.
  expect_true("T3" %in% res$targets$targetName)
})
