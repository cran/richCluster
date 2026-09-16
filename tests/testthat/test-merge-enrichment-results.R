# SPEC-RC-006 -- harden merge_enrichment_results()
#
# Covers REQ-006-001 (C1, single-dataset merge), REQ-006-002 (C2, missing-GeneID
# guard), REQ-006-003 (C3, single-column averaging), REQ-006-004 (C4,
# DatasetCount), REQ-006-005 (C5, column-name standardisation) and REQ-006-006
# (T1-10, gene-level dedup).
#
# The two demo frames are the package's own shipped extdata pair; every hard
# number below is the measured figure recorded in spec.md section B.

rc006_demo <- function(which_file) {
  read.delim(
    system.file("extdata", which_file, package = "richCluster"),
    stringsAsFactors = FALSE
  )
}

test_that("REQ-006-001: single-dataset merge returns the unsuffixed columns", {
  d1 <- rc006_demo("HF36wk_vs_HF12wk.txt")
  m <- merge_enrichment_results(list(d1))
  expect_true(all(c("GeneID", "Pvalue", "Padj") %in% names(m)))
  expect_true(all(c("GeneID_1", "Pvalue_1", "Padj_1") %in% names(m)))
  expect_identical(nrow(m), 3059L)
  # Clustering all 3059 terms takes ~30 s locally; CRAN caps a whole check at 10 min.
  skip_on_cran()
  expect_no_error(cluster(list(d1)))
})

test_that("REQ-006-001: the single-dataset return is row-ordered by Term", {
  d1 <- rc006_demo("HF36wk_vs_HF12wk.txt")
  m <- merge_enrichment_results(list(d1))
  # C12 (2026-09-05): the contract is canonical BYTE order, not the session's
  # collation.  is.unsorted() consults LC_COLLATE and so cannot express it --
  # under a case-folding collation it calls the canonical order unsorted.
  expect_identical(m$Term, sort(m$Term, method = "radix"))
  expect_setequal(m$Term, d1$Term)
})

test_that("REQ-006-004: DatasetCount reports how many datasets reported each term", {
  d1 <- rc006_demo("HF36wk_vs_HF12wk.txt")
  d2 <- rc006_demo("HF36wk_vs_WT12wk.txt")
  m <- merge_enrichment_results(list(d1, d2))
  expect_true("DatasetCount" %in% names(m))
  expect_identical(names(m)[ncol(m)], "DatasetCount")     # positioned last
  expect_type(m$DatasetCount, "integer")
  expect_false(anyNA(m$DatasetCount))
  expect_identical(sum(m$DatasetCount == 1L), 1541L)
  expect_identical(sum(m$DatasetCount == 2L), 1862L)
  expect_true(all(m$DatasetCount >= 1L & m$DatasetCount <= 2L))

  m1 <- merge_enrichment_results(list(d1))
  expect_true(all(m1$DatasetCount == 1L))
})

test_that("REQ-006-005: 'value' is not mapped, and collisions warn without duplicating", {
  expect_identical(unname(richCluster:::format_colnames(c("Term", "GeneID", "value"))),
                   c("Term", "GeneID", "value"))
  expect_warning(out <- richCluster:::format_colnames(c("Term", "Pvalue", "pval")),
                 class = "richCluster_colname_collision")
  expect_false(anyDuplicated(unname(out)) > 0)
  expect_identical(unname(out), c("Term", "Pvalue", "pval"))
})

test_that("REQ-006-005: the demo frames standardise without a collision warning", {
  d1 <- rc006_demo("HF36wk_vs_HF12wk.txt")
  d2 <- rc006_demo("HF36wk_vs_WT12wk.txt")
  expect_no_warning(richCluster:::format_colnames(colnames(d1)))
  expect_no_warning(richCluster:::format_colnames(colnames(d2)))
})

test_that("REQ-006-002: a frame set with no GeneID column errors diagnosably", {
  d1 <- rc006_demo("HF36wk_vs_HF12wk.txt")
  d2 <- rc006_demo("HF36wk_vs_WT12wk.txt")
  a <- d1[, c("Term", "Pvalue", "Padj")]
  b <- d2[, c("Term", "Pvalue", "Padj")]
  expect_error(merge_enrichment_results(list(a, b)), class = "richCluster_missing_geneid")
  err <- tryCatch(merge_enrichment_results(list(a, b)), error = function(e) e)
  expect_match(conditionMessage(err), "GeneID_1", fixed = TRUE)
  expect_false(grepl("undefined columns selected", conditionMessage(err), fixed = TRUE))
})

test_that("REQ-006-002: one available GeneID column is enough", {
  d1 <- rc006_demo("HF36wk_vs_HF12wk.txt")
  d2 <- rc006_demo("HF36wk_vs_WT12wk.txt")
  g1 <- d1[, c("Term", "GeneID", "Pvalue")]
  g2 <- d2[, c("Term", "Pvalue")]
  expect_no_error(m <- merge_enrichment_results(list(g1, g2)))
  expect_true("GeneID" %in% names(m))
})

test_that("REQ-006-006: GeneID carries each gene at most once, in first-appearance order", {
  d1 <- rc006_demo("HF36wk_vs_HF12wk.txt")
  d2 <- rc006_demo("HF36wk_vs_WT12wk.txt")
  m <- merge_enrichment_results(list(d1, d2))
  toks <- strsplit(m$GeneID, ",", fixed = TRUE)
  expect_true(all(vapply(toks, function(g) !anyDuplicated(g), logical(1))))
  expect_identical(max(lengths(toks)), 1639L)
  expect_identical(length(unique(unlist(toks))), 1819L)
  expect_lte(max(lengths(toks)), length(unique(unlist(toks))))
})

test_that("REQ-006-006: dedup keeps first-appearance order and drops empty tokens", {
  a <- data.frame(Term = "T1", GeneID = "A,B,C", Pvalue = 0.1, Padj = 0.1,
                  stringsAsFactors = FALSE)
  b <- data.frame(Term = "T1", GeneID = "B,C,D", Pvalue = 0.2, Padj = 0.2,
                  stringsAsFactors = FALSE)
  expect_identical(merge_enrichment_results(list(a, b))$GeneID, "A,B,C,D")

  e1 <- data.frame(Term = "T1", GeneID = "A,,B", Pvalue = 0.1, Padj = 0.1,
                   stringsAsFactors = FALSE)
  e2 <- data.frame(Term = "T1", GeneID = ",B,", Pvalue = 0.2, Padj = 0.2,
                   stringsAsFactors = FALSE)
  expect_identical(merge_enrichment_results(list(e1, e2))$GeneID, "A,B")
})

test_that("REQ-006-003: exactly one available Pvalue column does not error", {
  d1 <- rc006_demo("HF36wk_vs_HF12wk.txt")
  d2 <- rc006_demo("HF36wk_vs_WT12wk.txt")
  p1 <- d1[, c("Term", "GeneID", "Pvalue")]
  p2 <- d2[, c("Term", "GeneID")]
  expect_no_error(m <- merge_enrichment_results(list(p1, p2)))
  expect_true("Pvalue" %in% names(m))
  expect_identical(m$Pvalue, m$Pvalue_1)
})

test_that("REQ-006-003: exactly one available Padj column does not error", {
  d1 <- rc006_demo("HF36wk_vs_HF12wk.txt")
  d2 <- rc006_demo("HF36wk_vs_WT12wk.txt")
  q1 <- d1[, c("Term", "GeneID", "Padj")]
  q2 <- d2[, c("Term", "GeneID")]
  expect_no_error(m <- merge_enrichment_results(list(q1, q2)))
  expect_identical(m$Padj, m$Padj_1)
})

test_that("REQ-006-003: the two-column mean is the 1.0.2 arithmetic, unmoved", {
  d1 <- rc006_demo("HF36wk_vs_HF12wk.txt")
  d2 <- rc006_demo("HF36wk_vs_WT12wk.txt")
  m <- merge_enrichment_results(list(d1, d2))
  expect_identical(m$Pvalue, rowMeans(m[, c("Pvalue_1", "Pvalue_2")], na.rm = TRUE))
  expect_identical(m$Padj,   rowMeans(m[, c("Padj_1", "Padj_2")],   na.rm = TRUE))
})

test_that("REQ-006-004: the marker columns do not leak into the return", {
  d1 <- rc006_demo("HF36wk_vs_HF12wk.txt")
  d2 <- rc006_demo("HF36wk_vs_WT12wk.txt")
  m <- merge_enrichment_results(list(d1, d2))
  expect_false(any(grepl(".rc_contrib", names(m), fixed = TRUE)))
})
