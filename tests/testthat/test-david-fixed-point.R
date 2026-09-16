# tests/testthat/test-david-fixed-point.R
# SPEC-RC-009 / DS-03: mergeSeeds() must reach a fixed point of its own
# predicate.  Data: bundled extdata (ships in the tarball); no frozen fixtures.
test_that("DS-03: no pair of final DAVID clusters satisfies the merge predicate", {
  f1 <- read.delim(system.file("extdata", "HF36wk_vs_HF12wk.txt", package = "richCluster"),
                   stringsAsFactors = FALSE)
  f2 <- read.delim(system.file("extdata", "HF36wk_vs_WT12wk.txt", package = "richCluster"),
                   stringsAsFactors = FALSE)
  res <- NULL
  utils::capture.output(res <- david_cluster(list(head(f1, 120), head(f2, 120))))
  sets <- lapply(res$clusters$TermIndices,
                 function(s) as.integer(strsplit(s, ", ", fixed = TRUE)[[1]]))
  n <- length(sets)
  violations <- 0L
  for (i in seq_len(n - 1)) for (j in (i + 1):n) {
    m <- length(intersect(sets[[i]], sets[[j]]))
    if (m > 0 && 2 * m / (length(sets[[i]]) + length(sets[[j]])) >
          res$cluster_options$multiple_linkage_threshold) {
      violations <- violations + 1L
    }
  }
  # The merge predicate is STRICT (score > threshold, DavidClustering.cpp:149).
  # Pairs at exactly Dice == 0.5 are legitimate and remain; do NOT tighten to
  # >= (boundary conventions are DS-11, deferred -- SPEC-RC-009 REQ-009-04).
  expect_identical(violations, 0L)                       # as-found tree: 2
  expect_identical(length(sets), 25L)                    # as-found tree: 27
  # The fix merges clusters; it never adds or drops clustered terms.
  expect_identical(length(unique(unlist(sets))), 145L)   # unchanged by the fix
})
