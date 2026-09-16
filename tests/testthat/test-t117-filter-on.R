# tests/testthat/test-t117-filter-on.R
# ===========================================================================
# SPEC-RC-007 / T1-17 -- cluster() selects terms on a named column.
#
# At 1.0.2 cluster() filtered on the RAW Pvalue, hardcoded, while its own
# documentation named Padj as a required column.  Selecting enrichment terms on
# unadjusted p-values across thousands of tested terms is the wrong default.
#
# THE CRITICAL CONSTRAINT (spec.md SS-A.2, as authored): nothing that ran
# against 1.0.2 may error against the next release.  The release was renumbered
# to 2.0.0 on 2026-08-27, which lifts the semver obligation but NOT the
# behaviour: frames carrying Pvalue but no Padj ran fine at 1.0.2, so the
# missing-column path WARNS and completes because erroring on working input
# buys nothing.  V3 pins that, and rewriting it as expect_error() would encode
# the exact defect this item exists to prevent.
#
# The two demo-pair row counts (3403 -> 2863 at min_value = 0.1) are AC-017-1
# and AC-017-2.  They each run the full pipeline over thousands of terms, so
# they are gated behind RC_RUN_SLOW=1 rather than run on every suite pass; both
# were executed and recorded during the run phase.
# ===========================================================================

rc_t117_frames <- function(cols = c("Term", "GeneID", "Pvalue", "Padj"),
                           n_terms = 40, block = 6, overlap = 3) {
  pool  <- sprintf("GENE%03d", seq_len(n_terms * block))
  terms <- sprintf("TERM_%03d", seq_len(n_terms))
  genes <- vapply(seq_len(n_terms), function(i) {
    start <- (i - 1) * (block - overlap) + 1
    paste(pool[start:(start + block - 1)], collapse = ",")
  }, character(1))
  # Padj is deliberately coarser than Pvalue, so a filter on one is visibly a
  # different filter from the other.
  full <- data.frame(Term = terms, GeneID = genes,
                     Pvalue = rep(1e-6, n_terms),
                     Padj   = rep(c(1e-6, 0.5), length.out = n_terms),
                     stringsAsFactors = FALSE)
  list(full[, cols, drop = FALSE], full[, cols, drop = FALSE])
}

rc_t117_quiet <- function(expr) {
  out <- NULL
  utils::capture.output(out <- expr)
  out
}


test_that("V1: filter_on exists and defaults to Padj", {
  expect_true("filter_on" %in% names(formals(richCluster::cluster)))
  expect_identical(formals(richCluster::cluster)$filter_on, "Padj")
})


test_that("V2: the chosen column is the column actually filtered on", {
  x <- rc_t117_frames()
  padj   <- rc_t117_quiet(richCluster::cluster(x))                        # default
  pvalue <- rc_t117_quiet(richCluster::cluster(x, filter_on = "Pvalue"))

  # Padj is 1e-6 on every other term, so the default keeps half the rows.
  expect_identical(nrow(padj$merged_df),   20L)
  expect_identical(nrow(pvalue$merged_df), 40L)
})


test_that("V3: a missing filter column WARNS and completes (the API gate)", {
  nop <- rc_t117_frames(cols = c("Term", "GeneID", "Pvalue"))   # no Padj

  res <- NULL
  expect_warning(res <- rc_t117_quiet(richCluster::cluster(nop)),
                 class = "richCluster_filter_on_fallback")

  # AC-017-3 writes this as expect_s3_class(res, "list"); cluster() returns a
  # BARE list with no class attribute (measured: "Actual OO type: none"), at
  # 1.0.2 as well as now, so expect_s3_class can never hold and expect_type is
  # the assertion that carries the criterion's intent -- it completed and
  # returned the documented list.
  expect_type(res, "list")
  expect_false(is.null(res))
  expect_true(all(c("distance_matrix", "all_clusters") %in% names(res)))
  # It fell back to Pvalue, which keeps every row.
  expect_identical(nrow(res$merged_df), 40L)
})


test_that("V4: with neither column present it errors informatively", {
  neither <- rc_t117_frames(cols = c("Term", "GeneID"))
  expect_error(rc_t117_quiet(richCluster::cluster(neither)), "Pvalue")
})


test_that("V5: cluster_options records filter_on", {
  x <- rc_t117_frames()
  res <- rc_t117_quiet(richCluster::cluster(x))
  expect_identical(res$cluster_options$filter_on, "Padj")
  res2 <- rc_t117_quiet(richCluster::cluster(x, filter_on = "Pvalue"))
  expect_identical(res2$cluster_options$filter_on, "Pvalue")
})


test_that("V6: filter_on is validated against a closed set", {
  x <- rc_t117_frames()
  expect_error(richCluster::cluster(x, filter_on = "FDR"), "filter_on")
})


test_that("V7: david_cluster() gains no filter_on (REQ-017-9)", {
  expect_false("filter_on" %in% names(formals(richCluster::david_cluster)))
})


test_that("V8: the replacement filter is equivalent to the 1.0.2 dplyr form", {
  # 1.0.2 used dplyr::filter(Pvalue < min_value), which drops NA rows and
  # renumbers row names.  Base R `[` does neither, so the implementation adds
  # both explicitly.  This asserts the equivalence directly, including on NA.
  df <- data.frame(Term = sprintf("T%02d", 1:6),
                   Pvalue = c(0.01, NA, 0.2, 0.05, NA, 0.4),
                   stringsAsFactors = FALSE)
  min_value <- 0.1

  dplyr_form <- as.data.frame(dplyr::filter(df, Pvalue < min_value))
  keep <- !is.na(df[["Pvalue"]]) & df[["Pvalue"]] < min_value
  base_form <- df[keep, , drop = FALSE]
  rownames(base_form) <- NULL

  expect_identical(base_form, dplyr_form)
})


# SPEC-RC-004 REQ-RC004-021 -- DECLARED SKIP EXCEPTION (site S4).
# The RC_RUN_SLOW gate is ruled legitimate by SPEC-RC-007, which put these
# criteria behind it deliberately: they recompute the full 580-term fixture (or
# the whole pipeline over thousands of terms) and are minutes, not seconds.  It
# is an OPT-IN cost gate, not a fixture-availability gate, so REQ-RC004-015 does
# not bind it -- that requirement covers fixture, oracle and inst/extdata
# availability only.  The criteria still run for the developer with
# RC_RUN_SLOW=1 and in any check host willing to pay for them.
test_that("V9 (slow): the demo pair filters 3403 -> 2863 on Padj", {
  testthat::skip_if_not(identical(Sys.getenv("RC_RUN_SLOW"), "1"),
                        "AC-017-1 / AC-017-2 run the full pipeline over thousands of terms; set RC_RUN_SLOW=1")
  d <- list(
    read.delim(system.file("extdata", "HF36wk_vs_HF12wk.txt", package = "richCluster")),
    read.delim(system.file("extdata", "HF36wk_vs_WT12wk.txt", package = "richCluster")))

  expect_identical(nrow(rc_t117_quiet(
    richCluster::cluster(d, min_value = 0.1))$merged_df), 2863L)
  expect_identical(nrow(rc_t117_quiet(
    richCluster::cluster(d, min_value = 0.1, filter_on = "Pvalue"))$merged_df), 3403L)
})
