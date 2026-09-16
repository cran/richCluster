# tests/testthat/helper-t3-04a.R
# ===========================================================================
# SPEC-RC-004 M4 -- SHARED T3-04a RUNS, lifted out of test-david-kappa-matrix.R
# so that both the shipped file and the source-tree-only test-src-david-v102.R
# can reach them.  Helpers are sourced before every test file, so this resolves
# in both configurations.  The block below is moved VERBATIM.
#
# It also backs the v110-* check-time generator (author ruling 2026-08-26):
# rc_v110_recompute() reuses these memos, so generating v110-david / v110-dc
# costs NO extra computation beyond what the surviving V1/V3 blocks already
# perform in the same process.
# ===========================================================================


# --- Shared, memoised runs -------------------------------------------------
#
# runDavidClustering() on the 580-term fixture costs ~7 s and david_cluster()
# on the demo pair ~60 s, so each is computed at most once per process.
# Output suppression is done in exactly ONE place and is never nested:
# repeated .Call plus NESTED sink() in a single process is what produced the
# spurious "non-determinism" the revision-1 critic had to retract (T3-01).

rc_t3_04a_cache <- new.env(parent = emptyenv())

rc_t3_04a_quiet <- function(expr) {
  out <- NULL
  utils::capture.output(out <- expr)   # one sink level, never nested
  out
}

rc_t3_04a_memo <- function(key, compute) {
  if (!exists(key, envir = rc_t3_04a_cache, inherits = FALSE)) {
    assign(key, compute(), envir = rc_t3_04a_cache)
  }
  get(key, envir = rc_t3_04a_cache, inherits = FALSE)
}

# runDavidClustering(), frozen 580-term vectors, 0.5 / 3 / 3 / 0.5 -- the
# UNTOUCHED entry point, and the argument set the v102-david tag names.
rc_t3_04a_david <- function() {
  rc_t3_04a_memo("david", function() {
    dfx <- rc_david_fixture()
    rc_t3_04a_quiet(richCluster:::runDavidClustering(
      dfx$Term, dfx$GeneID,
      dfx$similarity_threshold, dfx$initial_group_membership,
      dfx$final_group_membership, dfx$multiple_linkage_threshold))
  })
}

# The ADDITIVE entry point on the same arguments: `clusters` (exactly what
# runDavidClustering() returns) plus `kappa_matrix`.
rc_t3_04a_david_with_kappa <- function() {
  rc_t3_04a_memo("david_with_kappa", function() {
    dfx <- rc_david_fixture()
    rc_t3_04a_quiet(richCluster:::runDavidClusteringWithKappa(
      dfx$Term, dfx$GeneID,
      dfx$similarity_threshold, dfx$initial_group_membership,
      dfx$final_group_membership, dfx$multiple_linkage_threshold))
  })
}

# david_cluster(list(d1, d2)) at documented defaults -- the OTHER entry point,
# and the argument set the v102-dc tag names.  rc_v102_demo_pair() is
# T3-03's definition of list(d1, d2); reused rather than restated so the two
# cannot drift.
rc_t3_04a_dc <- function() {
  rc_t3_04a_memo("dc", function() {
    rc_t3_04a_quiet(richCluster::david_cluster(rc_v102_demo_pair()))
  })
}


# --- SPEC-RC-004: the v110-* check-time generator --------------------------
#
# AUTHOR RULING 2026-08-26.  SPEC-RC-009 created the v110-* family (v110-david,
# v110-dc) and re-anchored test-david-kappa-matrix.R's V1(a) / V1(b) / V3 to it.
# That family post-dates SPEC-RC-004's authored generator shape list.  The
# author ruled that the generator SHALL also emit v110-david and v110-dc at
# check time, so those blocks stay in the SHIPPED suite and run at CRAN.
#
# Why this is admissible where regenerating v102-* is not: v110 is "what the
# current build produces", so it is regenerable BY CONSTRUCTION -- verified,
# both tags recompute identical() to their stored records.  v102-* records the
# 1.0.2 build and is not regenerable from corrected code at all (REQ-RC004-011,
# spec.md section 4.1).
#
# Honest limitation, recorded rather than hidden: in the TARBALL configuration
# the stored v110 files are absent, so V1(a)'s and V1(b)'s identical() clause
# compares the current build against itself.  There it degrades to a shape and
# determinism check; the names() clauses and V3's perturbation half still bind
# on real content.  In the SOURCE tree the stored record wins (stored-file
# precedence in rc_artifact_record()), so the assertion remains a genuine
# regression gate for the developer.  See SPEC-RC-004 progress.md section E.2.

#' Recompute the object a v110-* tag names, from its registry argument set.
#'
#' The ONE place that says how a v110 tag's object is produced.
rc_v110_recompute <- function(tag) {
  switch(rc_check_tag_string(tag),
    "v110-david" = rc_t3_04a_david(),
    "v110-dc"    = rc_t3_04a_dc(),
    stop(sprintf(paste0(
      "SPEC-RC-004: no v110 generator is defined for tag '%s'.  The v110 ",
      "family carries exactly two shapes, 'david' and 'dc' ",
      "(helper-artifacts.R rc_artifact_families())."), tag), call. = FALSE))
}
