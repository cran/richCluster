# tests/testthat/helper-t3-02.R
# ===========================================================================
# SPEC-RC-004 M4 -- SHARED T3-02 ORACLE RUNS, lifted out of test-kappa-oracle.R
# so that both the shipped file and the source-tree-only test-src-oracle-v102.R
# can reach them.  Helpers are sourced before every test file, so this resolves
# in both configurations.  The block below is moved VERBATIM.
#
# plan.md M4 anticipated exactly this hazard for rc_t3_04a_*() and prescribed a
# shared helper; the same hazard applies to rc_t3_02_*() because the two
# NULL-ACTION blocks that move to test-src-oracle-v102.R call rc_t3_02_oracle()
# and rc_t3_02_david_matrix().  Measured: without this file those two blocks
# error with "could not find function".
# ===========================================================================


# --- Shared, memoised runs -------------------------------------------------
#
# One process, one nesting level of output suppression, never nested:
# repeated .Call plus NESTED sink() in a single process is what produced the
# spurious "non-determinism" the revision-1 critic had to retract (T3-01).

rc_t3_02_cache <- new.env(parent = emptyenv())

rc_t3_02_memo <- function(key, compute) {
  if (!exists(key, envir = rc_t3_02_cache, inherits = FALSE)) {
    assign(key, compute(), envir = rc_t3_02_cache)
  }
  get(key, envir = rc_t3_02_cache, inherits = FALSE)
}

rc_t3_02_genes <- function() {
  rc_t3_02_memo("genes", function() rc_fixture_frozen()$GeneID)
}

# The shipping oracle: OD-5 clamp applied.
rc_t3_02_oracle <- function() {
  rc_t3_02_memo("oracle", function()
    rc_oracle_kappa_matrix(rc_t3_02_genes(), clamp = TRUE))
}

# The same oracle with the clamp lifted -- the raw chance-corrected kappa.
# This is what the independence check against DavidClustering compares, since
# DavidClustering predates OD-5 and applies no floor.
rc_t3_02_oracle_raw <- function() {
  rc_t3_02_memo("oracle_raw", function()
    rc_oracle_kappa_matrix(rc_t3_02_genes(), clamp = FALSE))
}

# The UNTOUCHED shipped computation, observed behind T3-04a's additive
# accessor.  DavidClustering::calculateKappaScores() counts the gene
# intersection correctly with a manual membership loop, so david_cluster() has
# been computing valid kappa all along -- it is the in-repo reference this
# oracle is validated against, and it predates every edit in this plan.
rc_t3_02_david_matrix <- function() {
  rc_t3_02_memo("david", function() {
    dfx <- rc_david_fixture()
    out <- NULL
    utils::capture.output(out <- richCluster:::david_kappa_matrix(
      dfx$Term, dfx$GeneID,
      dfx$similarity_threshold, dfx$initial_group_membership,
      dfx$final_group_membership, dfx$multiple_linkage_threshold))
    out
  })
}
