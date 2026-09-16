# SPEC-RC-011 REQ-011-007 -- permutation invariance of the C++ clustering core.
#
# Positive-control provenance: this exact configuration (580 bundled terms,
# kappa/average 0.5/0.5) FAILED on the pre-RC-011 build -- reverse, shuffle1 and
# shuffle2 all returned clusterings non-identical to the identity ordering
# (measured 2026-08-17; DS08-INVESTIGATION.md section 2.4, SPEC-RC-011 spec.md
# section B.4). A 156-row variant did NOT catch that defect (3/3 permutations
# identical on the defective build) and must not be substituted for this one.

test_that("runRichCluster is invariant to input row order (DS-08 / SPEC-RC-011)", {
  # Four 580-term runs, ~21 s locally; CRAN caps a check at 10 min.  Skipped there
  # rather than shrunk, because the smaller variant above misses the defect.
  skip_on_cran()
  quiet <- function(expr) { invisible(capture.output(r <- expr)); r }
  f1 <- read.delim(system.file("extdata", "HF36wk_vs_HF12wk.txt", package = "richCluster"),
                   stringsAsFactors = FALSE)
  f2 <- read.delim(system.file("extdata", "HF36wk_vs_WT12wk.txt", package = "richCluster"),
                   stringsAsFactors = FALSE)
  md <- quiet(merge_enrichment_results(list(f1, f2)))
  md <- md[md$Pvalue < 1e-4, ]
  tv <- md$Term
  gv <- md$GeneID

  key <- function(ord) {
    t2 <- tv[ord]
    r <- quiet(runRichCluster(t2, gv[ord], "kappa", 0.5, "average", 0.5))
    sets <- lapply(r$all_clusters$TermIndices,
                   function(s) sort(t2[as.integer(strsplit(s, ", ")[[1]]) + 1L]))
    sort(vapply(sets, paste, character(1), collapse = " ~ "))
  }

  ref <- key(seq_along(tv))
  expect_gt(length(ref), 0)

  set.seed(20260817)
  perms <- list(reverse  = rev(seq_along(tv)),
                shuffle1 = sample(seq_along(tv)),
                shuffle2 = sample(seq_along(tv)))
  for (nm in names(perms)) {
    expect_identical(key(perms[[nm]]), ref, info = paste("permutation:", nm))
  }
})
