# tests/testthat/test-fixture-generation.R
# ===========================================================================
# SPEC-RC-004 -- THE MIGRATION GUARD.
#
# The anchor fixture and the kappa oracle are now GENERATED from inst/extdata
# rather than read from tests/testthat/fixtures/, so that R CMD check on the
# built tarball runs them instead of skipping them.
#
# A generator can drift where a stored file cannot.  These tests are the fence:
# while the stored copies still exist in the source tree, the generated objects
# must be identical() to them.  This file must pass BEFORE any stored fixture
# is removed, and it keeps passing afterwards in the source tree.
#
# I8 NOTE (REQ-RC004-020, attributed re-baseline).  SPEC-RC-006's token-level
# gene deduplication inside merge_enrichment_results() has LANDED on this tree.
# Measured 2026-08-26: I8 moved 74018316 -> 51962367 and the merged frame
# gained a DatasetCount column, so the pre-dedup fx-input.rds ceased to be
# identical() to the generated fixture.  Per REQ-RC004-020 the literal below is
# the NEW value, recorded with its attribution in progress.md section E.2, and
# fx-input.rds was re-frozen in the same change so that REQ-RC004-009 holds.
# I1-I7 and I9 are unmoved; the oracle invariants O1-O6 are unmoved.
# ===========================================================================

test_that("the generated anchor fixture is identical to the stored one", {
  stored_path <- file.path(rc_fixture_dir(), "fx-input.rds")
  if (!file.exists(stored_path)) {
    testthat::succeed("no stored fixture in this configuration; nothing to compare")
    return(invisible(NULL))
  }
  expect_identical(rc_fixture_generate(), readRDS(stored_path))
})

test_that("the generated oracle matrix is identical to the stored one", {
  stored_path <- rc_oracle_crossteam_file()
  if (!file.exists(stored_path)) {
    testthat::succeed("no stored oracle in this configuration; nothing to compare")
    return(invisible(NULL))
  }
  expect_identical(rc_oracle_generate()$kappa, readRDS(stored_path)$kappa)
})

test_that("generation satisfies every SPEC-RC-004 invariant", {
  fx <- rc_fixture_frozen()
  expect_identical(length(fx$Term), 580L)
  expect_identical(length(unique(unlist(strsplit(fx$GeneID, ",")))), 1819L)
  expect_identical(nrow(fx$merged), 580L)
  # I8: the post-SPEC-RC-006 value.  See the file header for the attribution.
  expect_identical(sum(utf8ToInt(paste(fx$GeneID, collapse = "|"))), 51962367L)
  expect_identical(sum(utf8ToInt(paste(fx$Term,   collapse = "|"))), 1998341L)
})

test_that("generation is loud, not silent, when an input is missing", {
  # REQ-RC004-007 / REQ-RC004-015: a broken instrument FAILS, it does not skip.
  expect_error(rc_fixture_generate(extdata_dir = tempdir()))
})

# --- The v110 family, per the author ruling of 2026-08-26 ------------------
#
# SPEC-RC-009 re-anchored test-david-kappa-matrix.R's V1(a) / V1(b) / V3 to the
# v110-* family, which post-dates this SPEC's authored generator shape list.
# The author ruled that the generator SHALL also emit v110-david and v110-dc at
# check time, so those blocks stay in the shipped suite and run at CRAN.  v110
# is "what the current build produces" and is regenerable by construction --
# unlike the frozen v102-* historical record, which REQ-RC004-011 forbids
# regenerating.  This guard is the v110 half of the migration fence.

test_that("the generated v110 artifacts are identical to the stored ones", {
  for (tag in c("v110-david", "v110-dc")) {
    stored <- rc_artifact_path(tag)
    if (!file.exists(stored)) {
      testthat::succeed(paste("no stored", tag, "in this configuration"))
      next
    }
    expect_identical(rc_artifact_generate(tag)$object, readRDS(stored)$object)
  }
})

test_that("v110 generation is loud, not silent, for an ungeneratable tag", {
  # Only the v110 family is generatable.  A v102-* tag must NOT be regenerable
  # from the corrected build -- REQ-RC004-011, spec.md section 4.1.
  expect_false(rc_artifact_generatable("v102-david"))
  expect_error(rc_artifact_generate("v102-david"))
})
