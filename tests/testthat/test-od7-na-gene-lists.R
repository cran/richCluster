# tests/testthat/test-od7-na-gene-lists.R
# ===========================================================================
# OD-7 -- what a MISSING gene list means. Author decision (2026-08-30):
# a missing gene list is an EMPTY SET, and its exclusion is announced.
#
# An empty set is similar to nothing, so a term carrying one cannot join a
# cluster and min_terms drops it naturally. Nothing is removed from the
# caller's data and no row count changes; a warning names how many terms were
# affected so the outcome is never silent.
#
# SCOPE, MEASURED BEFORE FIXING. The converge ledger describes this as "two
# terms with missing gene lists report similarity 1.0 and cluster together",
# blocked on src/. That overstates the reach and the fix both:
#
#   cluster() / david_cluster()   merge_enrichment_results() already collapses
#                                 an all-NA gene list to "" via na.omit(), and
#                                 "" already yields similarity 0 in the C++.
#                                 Correct before this change.
#   NA in one dataset only        the other dataset's genes are used. Correct.
#   runRichCluster() direct       NA reaches Rcpp::as and becomes the literal
#                                 string "NA", so two such terms share their
#                                 single "gene" and score 1.0. THE defect.
#
# So the C++ needs no change: the empty-set behaviour it already implements is
# simply not reached for NA. Normalising at the R boundary is the whole fix.
# ===========================================================================

rc_od7_quiet <- function(expr) {
  out <- NULL
  utils::capture.output(out <- expr)
  out
}

rc_od7_frames <- function(g) {
  data.frame(Term = c("T1", "T2", "T3"), GeneID = g,
             Pvalue = rep(1e-6, 3), Padj = rep(1e-6, 3),
             stringsAsFactors = FALSE)
}


test_that("W1: two NA gene lists are NOT similar -- the empty set matches nothing", {
  r <- rc_od7_quiet(suppressWarnings(richCluster::runRichCluster(
    c("A", "B", "C"), c(NA_character_, NA_character_, "X,Y"),
    "jaccard", 0.5, "average", 0.5)))
  # Pre-fix this is 1: NA becomes the literal "NA", so A and B share a "gene".
  expect_identical(r$distance_matrix["A", "B"], 0)
})


test_that("W2: NA and the empty string are treated identically", {
  na_run <- rc_od7_quiet(suppressWarnings(richCluster::runRichCluster(
    c("A", "B", "C"), c(NA_character_, NA_character_, "X,Y"),
    "jaccard", 0.5, "average", 0.5)))
  empty_run <- rc_od7_quiet(suppressWarnings(richCluster::runRichCluster(
    c("A", "B", "C"), c("", "", "X,Y"),
    "jaccard", 0.5, "average", 0.5)))
  expect_identical(na_run$distance_matrix, empty_run$distance_matrix)
})


test_that("W3: the exclusion is announced, not silent, and names the count", {
  w <- tryCatch(
    rc_od7_quiet(richCluster::runRichCluster(
      c("A", "B", "C"), c(NA_character_, NA_character_, "X,Y"),
      "jaccard", 0.5, "average", 0.5)),
    warning = function(w) w)
  expect_s3_class(w, "richCluster_empty_gene_list")
  expect_match(conditionMessage(w), "2")            # the count, not a vague notice
  expect_match(conditionMessage(w), "gene list")
})


test_that("W4: cluster() announces it too, through the merge path", {
  a <- rc_od7_frames(c(NA, NA, "X,Y,Z"))
  b <- rc_od7_frames(c(NA, NA, "X,Y,W"))
  w <- tryCatch({ richCluster::merge_enrichment_results(list(a, b)); NULL },
                warning = function(w) w)
  expect_s3_class(w, "richCluster_empty_gene_list")
  expect_match(conditionMessage(w), "2")
})


test_that("W5: data with no missing gene lists stays silent", {
  a <- rc_od7_frames(c("P,Q", "R,S", "X,Y,Z"))
  b <- rc_od7_frames(c("P,Q", "R,S", "X,Y,W"))
  expect_no_warning(richCluster::merge_enrichment_results(list(a, b)))

  expect_no_warning(rc_od7_quiet(richCluster::runRichCluster(
    c("A", "B", "C"), c("P,Q", "R,S", "X,Y"),
    "jaccard", 0.5, "average", 0.5)))
})


test_that("W6: NA in one dataset only still uses the other dataset's genes", {
  a <- rc_od7_frames(c(NA, "P,Q", "X,Y"))
  b <- rc_od7_frames(c("P,Q", "P,Q", "X,Y"))
  m <- suppressWarnings(richCluster::merge_enrichment_results(list(a, b)))
  expect_identical(m$GeneID[m$Term == "T1"], "P,Q")
})
