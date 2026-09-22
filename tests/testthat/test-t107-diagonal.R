# tests/testthat/test-t107-diagonal.R
# ===========================================================================
# SPEC-RC-007 / T1-07 -- the EXPORTED distance-matrix diagonal is 1.
#
# OQ-1 was answered by the author on 2026-08-25 as READING 1: the value is
# written at export time, inside DistanceMatrix::export_r().  The in-memory
# matrix keeps richCluster::SAME_TERM_DISTANCE, the constant stays in
# src/RichCluster.h, and tools/cpp-units/rc_cpp_units.cpp's StubMatrix -- which
# is built around that in-memory convention -- stays valid and untouched, as
# test-src-cpp-units.R requires.
#
# Why 1 (spec.md SS-E.4): the oracle computes kappa(A, A) == 1; Jaccard and Dice
# self-similarity are 1 by definition; the downstream consumer in
# R/cluster_correlation.R patched to exactly 1 at two independent sites (both
# now removed, REQ-007-4); and NA is excluded because full_network() feeds the
# raw matrix to igraph::graph_from_adjacency_matrix().
#
# The off-diagonal no-op proof (AC-007-2) and the cluster identity (AC-007-3)
# assert against the pre-T1-07-* tags, which were emitted from the tree
# immediately BEFORE this change landed.  They are run-phase criteria, gated
# here behind RC_RUN_SLOW=1 because each recomputes the 580-term fixture.
# ===========================================================================

rc_t107_frames <- function(n_terms = 30, block = 6, overlap = 3) {
  pool  <- sprintf("GENE%03d", seq_len(n_terms * block))
  terms <- sprintf("TERM_%03d", seq_len(n_terms))
  genes <- vapply(seq_len(n_terms), function(i) {
    start <- (i - 1) * (block - overlap) + 1
    paste(pool[start:(start + block - 1)], collapse = ",")
  }, character(1))
  df <- data.frame(Term = terms, GeneID = genes,
                   Pvalue = rep(1e-6, n_terms), Padj = rep(1e-6, n_terms),
                   stringsAsFactors = FALSE)
  list(df, df)
}

rc_t107_quiet <- function(expr) {
  out <- NULL
  utils::capture.output(out <- expr)
  out
}

RC_T107_METRICS  <- c("kappa", "jaccard", "dice")
RC_T107_LINKAGES <- c("single", "complete", "average", "ward")


test_that("V1: the exported diagonal is exactly 1 for every metric x linkage", {
  x <- rc_t107_frames()
  for (metric in RC_T107_METRICS) {
    for (linkage in RC_T107_LINKAGES) {
      res <- rc_t107_quiet(richCluster::cluster(
        x, distance_metric = metric, linkage_method = linkage))
      dm <- res$distance_matrix
      expect_identical(unique(diag(dm)), 1,
                       info = paste(metric, linkage))
    }
  }
})


test_that("V2: the diagonal is not NA, and full_network() still builds", {
  x   <- rc_t107_frames()
  cl  <- rc_t107_quiet(richCluster::cluster(x))
  dm  <- cl$distance_matrix

  expect_false(anyNA(diag(dm)))
  expect_no_error(richCluster::full_network(cl))
})


test_that("V3: the change is confined to the diagonal", {
  # The off-diagonal is untouched at export: every entry still equals what the
  # in-memory matrix holds.  Asserted structurally here (the matrix is
  # symmetric and its off-diagonal carries real distances, not the sentinel);
  # the bit-identity form against the pre-T1-07 baseline is V5.
  x  <- rc_t107_frames()
  dm <- rc_t107_quiet(richCluster::cluster(x))$distance_matrix

  expect_true(isSymmetric(unname(dm)))
  ut <- dm[upper.tri(dm)]
  expect_false(any(ut == -99))
  expect_true(any(abs(ut) > 1e-8))
})


# SPEC-RC-004 REQ-RC004-021 -- DECLARED SKIP EXCEPTION (site S4).
# The RC_RUN_SLOW gate is ruled legitimate by SPEC-RC-007, which put these
# criteria behind it deliberately: they recompute the full 580-term fixture (or
# the whole pipeline over thousands of terms) and are minutes, not seconds.  It
# is an OPT-IN cost gate, not a fixture-availability gate, so REQ-RC004-015 does
# not bind it -- that requirement covers fixture, oracle and inst/extdata
# availability only.  The criteria still run for the developer with
# RC_RUN_SLOW=1 and in any check host willing to pay for them.
test_that("V4 (slow): the frozen fixture's diagonal is 1 at runRichCluster level", {
  testthat::skip_if_not(identical(Sys.getenv("RC_RUN_SLOW"), "1"),
                        "AC-007-1 over the 580-term fixture; set RC_RUN_SLOW=1")
  rc_skip_if_no_fixture()
  fx <- rc_fixture_frozen()
  for (metric in c("kappa", "jaccard")) {
    for (linkage in RC_T107_LINKAGES) {
      dm <- rc_t107_quiet(richCluster::runRichCluster(
        fx$Term, fx$GeneID, metric, 0.5, linkage, 0.5))$distance_matrix
      expect_identical(unique(diag(dm)), 1, info = paste(metric, linkage))
    }
  }
})


# SPEC-RC-004 REQ-RC004-021 -- DECLARED SKIP EXCEPTION (site S4).
# The RC_RUN_SLOW gate is ruled legitimate by SPEC-RC-007, which put these
# criteria behind it deliberately: they recompute the full 580-term fixture (or
# the whole pipeline over thousands of terms) and are minutes, not seconds.  It
# is an OPT-IN cost gate, not a fixture-availability gate, so REQ-RC004-015 does
# not bind it -- that requirement covers fixture, oracle and inst/extdata
# availability only.  The criteria still run for the developer with
# RC_RUN_SLOW=1 and in any check host willing to pay for them.
test_that("V5 (slow): off-diagonal and clusters are bit-identical to pre-T1-07", {
  testthat::skip_if_not(identical(Sys.getenv("RC_RUN_SLOW"), "1"),
                        "AC-007-2 / AC-007-3 recompute the 580-term fixture; set RC_RUN_SLOW=1")
  rc_skip_if_no_fixture()
  rc_skip_if_no_artifact("pre-T1-07-dm-kappa-average")
  rc_skip_if_no_artifact("pre-T1-07-clusters-kappa-average")
  fx  <- rc_fixture_frozen()
  res <- rc_t107_quiet(richCluster::runRichCluster(
    fx$Term, fx$GeneID, "kappa", 0.5, "average", 0.5))

  dm  <- res$distance_matrix
  old <- rc_artifact_object("pre-T1-07-dm-kappa-average")

  # The baseline is a genuine "before": its diagonal is still the sentinel.
  expect_identical(unique(diag(old)), -99)
  expect_identical(dm[upper.tri(dm)], old[upper.tri(old)])
  expect_identical(dm[lower.tri(dm)], old[lower.tri(old)])

  rc_expect_identical("pre-T1-07-clusters-kappa-average", res$all_clusters,
                      "runRichCluster()")
})
