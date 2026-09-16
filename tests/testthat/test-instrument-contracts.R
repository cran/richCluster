# tests/testthat/test-instrument-contracts.R
# ===========================================================================
# SPEC-RC-005 -- CONTRACTS THE WAVE-0 INSTRUMENT MUST HOLD.
#
# The plan documents that originally carried these contracts (PLAN_V2_*.md) are
# irrecoverably lost.  Every ruling that used to live in a plan therefore lives
# HERE, where it is executed rather than merely written down.
#
# Each block cites the reviewer finding it closes.
# ===========================================================================


# --- A1: the rc_emit signature is FROZEN at four arguments -----------------
#
# Helpers ship rc_emit(tag, object, entry_point, overwrite); the lost plan
# specified rc_emit(tag, object).  Reviewer A's ruling -- the implementation is
# right and the plan is wrong -- is adopted, and asserted here.

test_that("A1: rc_emit() has the frozen four-argument signature", {
  f <- formals(rc_emit)
  expect_identical(names(f), c("tag", "object", "entry_point", "overwrite"))
  expect_identical(f$overwrite, FALSE)
  # entry_point must carry NO default: a default makes the T3-01 V6
  # entry-point guard bypassable by omission.
  expect_true(identical(f$entry_point, quote(expr = )))
})

test_that("A1: rc_expect_identical() has the frozen signature", {
  f <- formals(rc_expect_identical)
  expect_identical(names(f), c("tag", "object", "entry_point", "args"))
  expect_null(f$args)
  expect_true(identical(f$entry_point, quote(expr = )))
})

test_that("A1: rc_emit_v102() has the frozen signature", {
  f <- formals(rc_emit_v102)
  expect_identical(names(f), c("tag", "object", "entry_point", "overwrite"))
  expect_identical(f$overwrite, FALSE)
  expect_true(identical(f$entry_point, quote(expr = )))
})



test_that("A4: the oracle record names its input_digest expression", {
  rc_skip_if_no_fixture()
  fx  <- rc_fixture_frozen()
  rec <- rc_oracle_record(matrix(1, 1, 1), fx$GeneID[1], clamp = TRUE)
  expect_true(is.character(rec$input_digest_expr) && nzchar(rec$input_digest_expr),
              info = "record carries an input_digest but does not name the expression digested")
  expect_true(grepl("utf8ToInt", rec$input_digest_expr, fixed = TRUE))
})


# --- A5: a missing frozen fixture is an ERROR, not a skip ------------------

# SPEC-RC-004 RE-ANCHOR (REQ-RC004-002 / -003 / -007 / -015).
# A5 previously asserted that rc_fixture_frozen() HARD-STOPS when the fixture
# directory is empty.  SPEC-RC-004 makes "fixture absent" unreachable: the
# fixture is GENERATED from inst/extdata, which ships, so an empty fixture
# directory now yields the generated fixture rather than an error.  Asserting
# the old hard stop would fail on correct work.
#
# The property A5 exists to protect is UNCHANGED in kind -- a missing
# instrument must never disappear quietly -- so it is re-anchored one level
# down, onto the input the generator cannot manufacture: with inst/extdata
# unreachable, generation ERRORS and does not skip.  SPEC-RC-005 REQ-RC-073/074
# ruled the hard stop correct precisely BECAUSE SPEC-RC-004 was going to make
# absence unreachable; this is that ruling arriving.
test_that("A5: an absent fixture is generated, and an absent INPUT still errors", {
  empty <- file.path(tempdir(), "rc005-empty-fixture-dir")
  dir.create(empty, showWarnings = FALSE, recursive = TRUE)
  old <- Sys.getenv("RC_FIXTURE_DIR", unset = NA_character_)
  on.exit({
    if (is.na(old)) Sys.unsetenv("RC_FIXTURE_DIR") else Sys.setenv(RC_FIXTURE_DIR = old)
  }, add = TRUE)

  Sys.setenv(RC_FIXTURE_DIR = empty)
  # No stored fixture: generation supplies it, and the object is the real one.
  fx <- rc_fixture_frozen()
  expect_identical(length(fx$Term), 580L)
  expect_identical(length(rc_david_fixture()$Term), 580L)

  # The remaining genuine absence: the bundled input itself.  This must ERROR,
  # naming the invariant or the file -- it must never skip (REQ-RC004-015).
  expect_error(rc_fixture_generate(extdata_dir = empty))
})


# --- A7: dice is not admissible in the stock (v102) family -----------------

test_that("A7: the v102 family rejects the dice metric", {
  expect_error(rc_artifact_registry("v102-dm-dice-average"), "dice")
  expect_error(rc_artifact_registry("v102-clusters-dice-single"), "dice")
})

test_that("A7: post-fix families still accept dice, and v102 still accepts kappa/jaccard", {
  expect_silent(rc_artifact_registry("w1-dm-dice-average"))
  expect_silent(rc_artifact_registry("w2-clusters-dice-complete"))
  expect_silent(rc_artifact_registry("t0-05-dm-dice-average"))
  expect_silent(rc_artifact_registry("v102-dm-kappa-average"))
  expect_silent(rc_artifact_registry("v102-dm-jaccard-average"))
  expect_true("dice" %in% RC_METRICS)
})


# --- B2: the Wave-0 instrument manifest ------------------------------------
#
# helper-characterization.R was outside its item's declared `files` list.  That
# list lived in a plan document that no longer exists, so the declaration lives
# here instead, where it is checked.
#
# SUBSET check, deliberately NOT set-equality: SPEC-RC-002, -003 and -004 each
# add test files, and a set-equality assertion would fail on every sibling SPEC.

RC_WAVE0_MANIFEST <- c(
  "helper-artifacts.R",
  "helper-characterization.R",
  "helper-fixture.R",
  "helper-oracle.R",
  "test-cpp-units.R",
  "test-david-kappa-matrix.R",
  "test-kappa-oracle.R",
  "test-instrument-contracts.R",
  "test-news-disclosure.R",
  "test-review-gate.R",
  "review-ledger.csv"
)

# SPEC-RC-004 M4 split the Wave-0 manifest in two.  test-characterization.R was
# renamed test-src-characterization.R and is excluded from the tarball by
# .Rbuildignore (REQ-RC004-017), so asserting its presence unconditionally would
# FAIL in the shipped configuration -- converting a currently-passing block into
# a failure, which is the opposite of this SPEC's purpose.  Splitting the
# manifest keeps B2 running in BOTH configurations without a skip
# (REQ-RC004-015): the shipped entries are asserted always, and the source-only
# entries are asserted whenever we are in a source tree, detected by the
# presence of any test-src-*.R file rather than by probing the very file whose
# presence is in question.
RC_WAVE0_MANIFEST_SRC <- c(
  "test-src-characterization.R"
)

test_that("B2: every declared Wave-0 instrument file is present", {
  missing <- RC_WAVE0_MANIFEST[
    !file.exists(vapply(RC_WAVE0_MANIFEST, testthat::test_path, character(1)))]
  expect_identical(missing, character(0),
                   info = paste("declared but absent:", paste(missing, collapse = ", ")))

  in_source_tree <- length(list.files(testthat::test_path("."),
                                      pattern = "^test-src-.*\\.R$")) > 0L
  if (in_source_tree) {
    missing_src <- RC_WAVE0_MANIFEST_SRC[
      !file.exists(vapply(RC_WAVE0_MANIFEST_SRC, testthat::test_path, character(1)))]
    expect_identical(missing_src, character(0),
                     info = paste("declared source-only but absent:",
                                  paste(missing_src, collapse = ", ")))
  } else {
    # Shipped configuration: the source-only files are absent BY DECLARATION.
    expect_identical(
      RC_WAVE0_MANIFEST_SRC[
        file.exists(vapply(RC_WAVE0_MANIFEST_SRC, testthat::test_path, character(1)))],
      character(0))
  }
})


# --- B3: build state is DERIVED, never asserted in prose -------------------
#
# Records claimed a "stock" build in prose while deriving
# package_version = 1.1.0 from the installed package -- the record contradicted
# itself.  The version now comes only from the derived field.
#
# Scope is the `build` field ONLY.  component / note / args / the
# characterization header legitimately mention 1.0.2 and are untouched.

test_that("B3: no registry build string asserts a version number", {
  builds <- rc_artifact_families()$build
  offenders <- builds[grepl("[0-9]+\\.[0-9]+\\.[0-9]+", builds)]
  expect_identical(unique(offenders), character(0),
                   info = paste("build prose asserting a version:",
                                paste(unique(offenders), collapse = " | ")))
})

