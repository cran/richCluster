# tests/testthat/helper-fixture.R
# ===========================================================================
# T3-01 (WAVE 0) -- THE ANCHOR FIXTURE.  Both teams import this file verbatim.
#
# Revision 1 anchored ~12 exact-number acceptance assertions on "the shipped
# demo fixture (n=580)" without ever defining it (critique C1).  This file is
# that definition, and it is an ARTIFACT both teams import -- not a constant
# either team re-derives.  A hand-copied constant is exactly what failed.
#
# Facts an implementer must know (all verified):
#   * min_value = 1e-4 is NOT the cluster() default.  The documented default
#     0.1 filters nothing here (max Pvalue = 0.04990348) and yields n = 3403.
#   * Row order is base::merge's Term-sorted order; there are 0 duplicated
#     terms in the 3403 rows.
#   * N = 1819 distinct genes; nrow = 580.
#   * The fixture is consumed by runRichCluster() DIRECTLY with these vectors.
#     It does not pass through cluster()'s filter, so T1-17 / D-06
#     (filter_on = "Padj" takes the 0.1 filter from 3403 terms to 2863)
#     cannot perturb it.
#
# Measurement discipline: run ONE linkage method and ONE permutation per R
# process.  Repeated .Call plus nested sink() in a single process produced a
# spurious "non-determinism" that the revision-1 critic had to retract.
#
# Depends on rc_fixture_dir(), defined in helper-artifacts.R.  testthat sources
# every helper-*.R before running any test, so call-time lookup always resolves.
# ===========================================================================


# --- (a) The one canonical fixture definition ------------------------------
#
# VERBATIM from PLAN_V2_T3.md T3-01(a).  This is the DEFINITION and it
# RECOMPUTES.  Tests must not call it; they load the frozen copy via
# rc_fixture_frozen() -- see (b) below and rc_freeze_fixture().
# --- SPEC-RC-004: the assertion helper -------------------------------------

#' Assert a generator invariant, naming what was expected and what was seen.
#' Errors -- never skips.  REQ-RC004-006 / REQ-RC004-007.
rc_assert_invariant <- function(label, observed, expected, tolerance = 0) {
  ok <- if (tolerance > 0) {
    is.numeric(observed) && length(observed) == length(expected) &&
      all(abs(observed - expected) <= tolerance)
  } else {
    identical(observed, expected)
  }
  if (!isTRUE(ok)) {
    stop(sprintf(paste0(
      "SPEC-RC-004 generator invariant '%s' VIOLATED.\n",
      "  expected: %s\n  observed: %s\n",
      "The generated fixture is NOT the frozen anchor fixture.  Do not proceed: ",
      "a fixture that has moved silently invalidates every downstream artifact ",
      "identity.  If a sibling SPEC changed merge_enrichment_results() on ",
      "purpose, see SPEC-RC-004 REQ-RC004-020 for the attributed re-baseline ",
      "procedure -- do NOT just update the literal."), label,
      format(expected), format(observed)), call. = FALSE)
  }
  invisible(TRUE)
}

# --- SPEC-RC-004: the anchor fixture generator -----------------------------

#' Generate the anchor fixture from inst/extdata.
#'
#' The ONE canonical definition.  stringsAsFactors is pinned explicitly
#' (REQ-RC004-008): the default flipped at R 4.0.0 and DESCRIPTION declares
#' R (>= 3.5.0), so an unpinned read.delim() yields a factor Term column on old
#' R and the object stops being identical() to the frozen fixture.  Row names
#' are NOT reset (REQ-RC004-014) -- subsetting leaves them non-sequential
#' ("11", "13", "18", ...) and they are part of the object's identity.
#'
#' I8 is the ONE invariant SPEC-RC-006's token-level gene deduplication moves.
#' Measured 2026-08-26 on this tree, with SPEC-RC-006 landed: 74018316 ->
#' 51962367.  The literal below is the post-dedup value, recorded with its
#' attribution in SPEC-RC-004 progress.md section E.2 per REQ-RC004-020, and
#' fx-input.rds was re-frozen in the same change so REQ-RC004-009 still holds.
rc_fixture_generate <- function(extdata_dir = NULL) {
  get1 <- function(f) {
    p <- if (is.null(extdata_dir)) {
      system.file("extdata", f, package = "richCluster")
    } else {
      file.path(extdata_dir, f)
    }
    if (!nzchar(p) || !file.exists(p)) {
      stop(sprintf(paste0(
        "SPEC-RC-004: bundled input '%s' not found.  The anchor fixture is ",
        "generated from inst/extdata, which ships in the tarball; its absence ",
        "is a real failure, not a reason to skip."), f), call. = FALSE)
    }
    p
  }

  d1 <- read.delim(get1("HF36wk_vs_HF12wk.txt"), stringsAsFactors = FALSE)
  d2 <- read.delim(get1("HF36wk_vs_WT12wk.txt"), stringsAsFactors = FALSE)
  rc_assert_invariant("I1 nrow(HF36wk_vs_HF12wk)", nrow(d1), 3059L)
  rc_assert_invariant("I2 nrow(HF36wk_vs_WT12wk)", nrow(d2), 2206L)

  m <- merge_enrichment_results(list(d1, d2))
  rc_assert_invariant("I3 nrow(merged)",       nrow(m), 3403L)
  rc_assert_invariant("I4 duplicated terms",   sum(duplicated(m$Term)), 0L)
  # C12 (2026-09-05): the merged order is canonical BYTE order, which is the
  # same in every locale.  A bare sort() consults LC_COLLATE and so cannot
  # express the invariant -- under a case-folding collation it calls the
  # canonical order unsorted.
  rc_assert_invariant("I5 merged Term-sorted",
                      identical(m$Term, sort(m$Term, method = "radix")), TRUE)

  fx <- m[m$Pvalue < 1e-4, ]                       # row names deliberately kept
  rc_assert_invariant("I6 nrow(fixture)", nrow(fx), 580L)
  rc_assert_invariant("I7 gene universe",
    length(unique(unlist(strsplit(fx$GeneID, ",")))), 1819L)
  rc_assert_invariant("I8 GeneID digest",
    sum(utf8ToInt(paste(fx$GeneID, collapse = "|"))), 51962367L)
  rc_assert_invariant("I9 Term digest",
    sum(utf8ToInt(paste(fx$Term, collapse = "|"))), 1998341L)

  list(Term = fx$Term, GeneID = fx$GeneID, merged = fx)
}

# Backward-compatible alias: rc_fixture() was the Wave-0 name.
rc_fixture <- function() rc_fixture_generate()

# --- SPEC-RC-004: the memo -------------------------------------------------
#
# REQ-RC004-012: each generated object is computed at most once per R process.
.rc_memo <- new.env(parent = emptyenv())


# --- (b) Freeze it ---------------------------------------------------------
#
# rc_fixture() is run ONCE in Wave 0 and the result saved to
# tests/testthat/fixtures/fx-input.rds.  EVERY LATER TEST LOADS THE FROZEN
# FILE, not the recomputed function.  merge_enrichment_results() is edited by
# T1-10, T1-11, T1-12 and T1-13; a recomputed fixture would silently move n,
# move N, move every kappa and move every artifact identity the instant one of
# those lands.

#' Path of the frozen anchor fixture.
rc_fixture_file <- function() {
  file.path(rc_fixture_dir(), "fx-input.rds")
}

#' TRUE when the frozen anchor fixture is present.
#'
#' tests/testthat/fixtures is .Rbuildignore'd (the frozen fixture is 0.65 MB and
#' the T3-03 / w1 / w2 artifact snapshots are far larger), so a suite run from an
#' installed tarball has no fixtures.  Downstream tests guard with
#' rc_skip_if_no_fixture().
rc_have_fixture <- function() {
  file.exists(rc_fixture_file())
}

#' Load the frozen anchor fixture.  This is what every later test uses.
#'
#' @section Why this ERRORS rather than skipping (SPEC-RC-005, finding A5):
#' A reviewer observed that this hard-stops when the fixture is absent and
#' suggested a skip.  That suggestion is SUPERSEDED and must not be applied.
#'
#' SPEC-RC-004 makes fixture GENERATION mandatory at check time: the two
#' signal-carrying fixtures are generated from inst/extdata rather than shipped.
#' Once fixtures are generated rather than stored, "fixtures absent" ceases to
#' be a reachable state -- the generator ran, or it errored, and either way
#' there is no legitimate absence left to skip on.  A hard failure is then
#' correct, and a skip would create exactly the silent-disappearance defect
#' that reviewer finding B1 names.
#' SPEC-RC-004 REQ-RC004-002 / -003 / -012: the stored file wins where it
#' exists (source tree, bit-for-bit Wave-0 semantics preserved); otherwise the
#' fixture is GENERATED from inst/extdata, which ships.  "Fixture absent" is
#' therefore no longer a reachable state -- the generator ran, or it errored.
rc_fixture_frozen <- function() {
  path <- rc_fixture_file()
  if (file.exists(path)) return(readRDS(path))          # source tree wins
  if (is.null(.rc_memo$fixture)) .rc_memo$fixture <- rc_fixture_generate()
  .rc_memo$fixture
}

#' REQ-RC004-016: the old skip guard becomes a hard requirement.
rc_require_fixture <- function() {
  invisible(rc_fixture_frozen())
}

#' Retained name, now FAIL-NOT-SKIP.  REQ-RC004-015 / REQ-RC004-016.
#'
#' The name is KEPT deliberately so that no call site changes: all 30 existing
#' callers keep working and simply stop being skips.  This single edit is what
#' recovers the fixture-gated blocks that R CMD check used to skip.
rc_skip_if_no_fixture <- function() rc_require_fixture()

#' Freeze the anchor fixture.  Wave 0, run once.
#'
#' Refuses to overwrite an existing frozen fixture whose content differs, unless
#' overwrite = TRUE is passed explicitly: silently re-freezing after T1-10 /
#' T1-11 / T1-12 / T1-13 have edited merge_enrichment_results() is precisely the
#' failure (b) exists to prevent.
rc_freeze_fixture <- function(overwrite = FALSE) {
  dir <- rc_fixture_dir()
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  path <- rc_fixture_file()
  fx <- rc_fixture()
  if (file.exists(path) && !isTRUE(overwrite)) {
    if (!identical(readRDS(path), fx)) {
      stop(sprintf(
        paste0("Refusing to overwrite the frozen anchor fixture '%s': the ",
               "recomputed fixture differs from the frozen one.  This is the ",
               "T3-01(b) guard.  Pass overwrite = TRUE only if you intend to ",
               "move every downstream artifact identity."), path), call. = FALSE)
    }
    return(invisible(path))
  }
  saveRDS(fx, path)
  invisible(path)
}


# --- (c) The DAVID-path fixture contract -----------------------------------
#
# The DAVID path needs its own fixture because david_cluster() applies NO
# p-value filter, so "david_cluster() on the fixture" is ambiguous between two
# very different inputs, and the DAVID result is highly sensitive to all four of
# its parameters.  Both the term vector AND the parameters are fixed here.
#
# Measured at exactly these documented defaults via runDavidClustering() on the
# frozen 580-term vectors: 90 CLUSTERS COVERING 523 TERMS.
#
# Two other readings that have circulated, recorded so neither is mistaken for
# the baseline:
#   igm = fgm = 4            -> 77 / 501  (the source of the previously quoted
#                                          "77 clusters covering 501 terms",
#                                          which is NOT a documented-defaults
#                                          figure)
#   igm = 3, fgm = 5         -> 67 / 486
#   david_cluster(list(d1,d2)) at documented defaults
#                            -> 557 clusters covering 2922 of 3403 terms,
#                               because it filters nothing.
# Of these, two are frozen under two different tags because the DAVID path has
# two entry points and items assert at both:
#   v102-david  runDavidClustering() on the 580-term / 0.5-3-3-0.5 reading
#               (the reading THIS contract pins, and the one the oracle gate uses)
#   v102-dc     david_cluster(list(d1, d2)) at documented defaults (557 / 2922)
# The igm = fgm = 4 and igm = 3, fgm = 5 readings are frozen as NOTHING.
#
# NOTE what this returns: Term / GeneID VECTORS plus four scalars.  That is an
# argument set for runDavidClustering(terms, geneIDs, 0.5, 3, 3, 0.5), and it
# CANNOT BE PASSED TO david_cluster() AT ALL, whose first formal
# `enrichment_results` is a list of data frames.  Any verification phrased as
# "david_cluster() on the DAVID fixture" is not runnable as written.
rc_david_fixture <- function() {
  fx <- rc_fixture_frozen()                # the SAME frozen 580-term vectors
  list(Term = fx$Term, GeneID = fx$GeneID,
       similarity_threshold       = 0.50,  # documented default
       initial_group_membership   = 3,     # documented default
       final_group_membership     = 3,     # documented default
       multiple_linkage_threshold = 0.50)  # documented default
}
