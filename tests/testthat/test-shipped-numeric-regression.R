# End-to-end numeric regression on the SHIPPED example data.
#
# WHY THIS FILE EXISTS.  Every other test in this suite checks a mechanism: a
# metric's closed form, a guard's error text, an invariant of the instrument.
# None of them assert what the package actually RETURNS on the data it ships.
# That gap is not hypothetical -- it has been measured twice:
#
#   * SPEC-RC-011 recorded that its seed-iteration change moved the 580-term
#     result from 40 final clusters to 42, and the suite passed unchanged
#     because no test named those figures.
#   * On 2026-08-27 NEWS.md was found stating 41 final clusters where the tree
#     produced 30.  The release notes of a release whose entire subject is
#     "the old numbers were wrong" carried a wrong number, and nothing failed.
#
# So this file pins the numbers themselves.  It uses inst/extdata, which ships,
# rather than a generated fixture, so it runs identically in the source tree and
# inside R CMD check on the tarball.  The three clustering tests skip on CRAN,
# which caps a whole check at 10 min (together ~40 s locally); set NOT_CRAN=true
# to run them under R CMD check.
#
# WHEN THIS FILE FAILS, IT HAS DONE ITS JOB.  A failure means the clustering
# output moved.  That is sometimes correct -- a deliberate fix will move it --
# but it must never move SILENTLY.  The response is:
#   1. establish WHY the number moved, and
#   2. update the expectation here AND every other place the figure is quoted,
#      NEWS.md and tools/baseline/results/ included, in the same commit.
# Never edit an expectation here to make a red suite green without doing (1).

demo_inputs <- function() {
  ed <- system.file("extdata", package = "richCluster")
  paths <- file.path(ed, c("HF36wk_vs_HF12wk.txt", "HF36wk_vs_WT12wk.txt"))
  skip_if_not(all(file.exists(paths)), "bundled extdata not available")
  lapply(paths, utils::read.delim, stringsAsFactors = FALSE)
}

# Quiet: the C++ core narrates every merge iteration to stdout.
run_quiet <- function(...) {
  res <- NULL
  invisible(utils::capture.output(res <- richCluster::cluster(...)))
  res
}

test_that("the shipped inputs merge to a known shape", {
  d <- demo_inputs()
  expect_equal(nrow(d[[1]]), 3059L)
  expect_equal(nrow(d[[2]]), 2206L)
  m <- merge_enrichment_results(d)
  expect_equal(nrow(m), 3403L)
  # C6's negative control: 0 duplicated terms, so the cross-join guard is silent
  # here.  If this ever fires, the C6 warning path is being exercised and the
  # figures below are measuring a different frame than the one they were frozen
  # against.
  expect_equal(sum(duplicated(m$Term)), 0L)
})

test_that("filter_on selects the term populations the release notes quote", {
  skip_on_cran()
  d <- demo_inputs()
  # These two counts are quoted verbatim in NEWS.md ("580 x 580 terms to
  # 376 x 376").  Pinning them here is what keeps that sentence honest.
  r_pv <- run_quiet(d, df_names = c("A", "B"), min_value = 1e-4,
                    filter_on = "Pvalue", distance_metric = "kappa",
                    linkage_method = "average")
  r_pa <- run_quiet(d, df_names = c("A", "B"), min_value = 1e-4,
                    filter_on = "Padj", distance_metric = "kappa",
                    linkage_method = "average")
  expect_equal(nrow(r_pv$distance_matrix), 580L)
  expect_equal(nrow(r_pa$distance_matrix), 376L)
  # ... and so are these two cluster counts.
  expect_equal(nrow(r_pv$final_clusters), 42L)
  expect_equal(nrow(r_pa$final_clusters), 30L)
})

test_that("every metric x linkage combination returns its frozen counts", {
  skip_on_cran()
  d <- demo_inputs()
  # min_value = 1e-4 with the Padj default: the 376-term configuration the
  # living baseline (tools/baseline/results/fixed_580_summary.csv) records.
  # dice has no baseline row -- it postdates that file (SPEC-RC-002) -- so its
  # figures are frozen here from a measured run on 2026-08-27.
  #
  # The two ward rows were re-measured on 2026-08-28 after the ESS increment was
  # corrected for OVERLAPPING clusters: the Lance-Williams closed form
  # nA*nB/(nA+nB) assumes A and B are disjoint, but mergeClusters unions them.
  # The other six rows and both distance-matrix fingerprints were unaffected,
  # which is the check that the change stayed inside ward().
  #
  # The jaccard/average row was re-measured on 2026-09-04 (C9, converge ledger
  # session 2026-09-04): cluster member sets became std::set, so average()
  # sums in canonical term order on every platform instead of the standard
  # library's hash order.  101/33 -> 102/34.  An exact-rational replay of the
  # seed-and-merge algorithm (no rounding) produces the 102-cluster partition,
  # so the old 101 was a floating-point artifact of the hash-order summation.
  # The other seven rows are unchanged, which is the check that the change is
  # confined to exact ties.
  expected <- list(
    "kappa/single"     = c(all =  42L, final =  12L),
    "kappa/complete"   = c(all = 131L, final =  59L),
    "kappa/average"    = c(all =  73L, final =  30L),
    "kappa/ward"       = c(all =  77L, final =  32L),
    "jaccard/single"   = c(all =  74L, final =  18L),
    "jaccard/complete" = c(all = 152L, final =  41L),
    "jaccard/average"  = c(all = 102L, final =  34L),
    "jaccard/ward"     = c(all = 117L, final =  33L)
  )
  for (key in names(expected)) {
    parts <- strsplit(key, "/", fixed = TRUE)[[1]]
    r <- run_quiet(d, df_names = c("A", "B"), min_value = 1e-4,
                   distance_metric = parts[1], linkage_method = parts[2])
    expect_equal(nrow(r$all_clusters),   expected[[key]][["all"]],
                 info = paste(key, "all_clusters"))
    expect_equal(nrow(r$final_clusters), expected[[key]][["final"]],
                 info = paste(key, "final_clusters"))
  }
})

test_that("the distance matrix carries a stable fingerprint", {
  skip_on_cran()
  d <- demo_inputs()
  # A checksum over the whole matrix catches a numeric drift that leaves the
  # cluster COUNTS intact -- the shape of change a count-only assertion misses.
  # Values are the measured sums at the 376-term configuration.
  expected_sum <- c(kappa   = 16360.804889050281,
                    jaccard = 16399.996191789854)
  for (mt in names(expected_sum)) {
    r <- run_quiet(d, df_names = c("A", "B"), min_value = 1e-4,
                   distance_metric = mt, linkage_method = "average")
    expect_equal(sum(r$distance_matrix), expected_sum[[mt]],
                 tolerance = 1e-9, info = paste(mt, "distance-matrix sum"))
    # Range invariants: every metric exports on [0, 1] with a unit diagonal.
    expect_gte(min(r$distance_matrix), 0)
    expect_lte(max(r$distance_matrix), 1)
    expect_equal(unique(diag(r$distance_matrix)), 1)
  }
})
