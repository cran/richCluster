# tests/testthat/test-t006-gene-delim.R
# ===========================================================================
# SPEC-RC-007 / T0-06 -- the gene delimiter is a parameter, not a literal.
#
# At 1.0.2 the delimiter is hardcoded "," in three places: the per-term split
# (src/RichCluster.cpp) twice, and the distinct-gene universe count
# (StringUtils::countUniqueElements).  A user whose GeneID column is ";"- or
# "/"-separated gets ONE "gene" per term, every term disjoint, every kappa ~0,
# and clustering silently returns nothing -- no error, no warning.
#
# The half-threading trap (spec.md SS-C.3) is the reason V3 exists.  Threading
# the delimiter into the split but NOT into the universe count leaves N counted
# under "," while |A| and |B| are counted under the user's delimiter: the 2x2
# agreement table is then inconsistent with its own universe and produces
# plausible-looking but wrong kappa.  V3 is the ONLY assertion here that fails
# in that state.  Never drop or weaken it.
#
# REQ-006-6: gene_delim is deliberately NOT added to david_cluster().  V6 pins
# that exclusion.
# ===========================================================================

# A two-frame enrichment input whose GeneID column uses `delim`.  Pvalue and
# Padj carry the same values so the frames behave identically before and after
# T1-17 changes which column cluster() filters on.
rc_t006_frames <- function(delim = ",", n_terms = 40, block = 6, overlap = 3) {
  pool  <- sprintf("GENE%03d", seq_len(n_terms * block))
  terms <- sprintf("TERM_%03d", seq_len(n_terms))
  genes <- vapply(seq_len(n_terms), function(i) {
    start <- (i - 1) * (block - overlap) + 1
    paste(pool[start:(start + block - 1)], collapse = delim)
  }, character(1))
  mk <- function(shift) {
    data.frame(Term = terms, GeneID = genes,
               Pvalue = rep(1e-6, n_terms), Padj = rep(1e-6, n_terms),
               stringsAsFactors = FALSE)
  }
  list(mk(0), mk(1))
}

rc_t006_quiet <- function(expr) {
  out <- NULL
  utils::capture.output(out <- expr)
  out
}


test_that("V1: cluster() exposes gene_delim, defaulting to a comma", {
  expect_true("gene_delim" %in% names(formals(richCluster::cluster)))
  expect_identical(formals(richCluster::cluster)$gene_delim, ",")
})


test_that("V2: the live runRichCluster carries geneDelim (duplicate-definition trap)", {
  # runRichCluster is defined TWICE -- R/RcppExports.R (generated) and
  # R/cluster.R (hand-written) -- and cluster.R WINS under C-locale collation.
  # Regenerating the bindings without updating the hand-written wrapper makes
  # gene_delim vanish from the R-visible function with no error anywhere.
  expect_true("geneDelim" %in% names(formals(richCluster::runRichCluster)))
  expect_identical(formals(richCluster::runRichCluster)$geneDelim, ",")
})


test_that("V3: the delimiter reaches the universe count too (half-threading catch)", {
  skip_on_cran()   # two 580-term runs, ~11 s locally; CRAN caps a check at 10 min
  rc_skip_if_no_fixture()
  fx <- rc_fixture_frozen()

  semi <- rc_t006_quiet(richCluster::runRichCluster(
    fx$Term, gsub(",", ";", fx$GeneID, fixed = TRUE),
    "kappa", 0.5, "average", 0.5, geneDelim = ";"))
  comma <- rc_t006_quiet(richCluster::runRichCluster(
    fx$Term, fx$GeneID, "kappa", 0.5, "average", 0.5, geneDelim = ","))

  # Same gene sets, same universe, only the separator byte differs -- so the
  # matrices must agree bit for bit.  Under half-threading N is counted under
  # "," (giving 580, one "gene" per term) while |A| and |B| are counted under
  # ";" (giving the real sizes), and these diverge.
  expect_identical(semi$distance_matrix, comma$distance_matrix)
  expect_identical(semi$all_clusters, comma$all_clusters)
})


test_that("V4: gene_delim actually changes the result, and is not degenerate", {
  x_semi <- rc_t006_frames(delim = ";")

  right <- rc_t006_quiet(richCluster::cluster(x_semi, gene_delim = ";"))
  wrong <- rc_t006_quiet(richCluster::cluster(x_semi))

  expect_false(identical(right$distance_matrix, wrong$distance_matrix))

  ut <- right$distance_matrix[upper.tri(right$distance_matrix)]
  expect_true(any(abs(ut) > 1e-8))          # the correctly-split run has signal

  # The mis-split run is the defect this item exists to remove: one "gene" per
  # term, so every term is disjoint from every other.
  utw <- wrong$distance_matrix[upper.tri(wrong$distance_matrix)]
  expect_true(all(abs(utw) < 1e-8))
})


test_that("V5: the default preserves comma behaviour exactly", {
  x_comma <- rc_t006_frames(delim = ",")
  a <- rc_t006_quiet(richCluster::cluster(x_comma))
  b <- rc_t006_quiet(richCluster::cluster(x_comma, gene_delim = ","))
  expect_identical(a$distance_matrix, b$distance_matrix)
  expect_identical(a$all_clusters, b$all_clusters)
})


test_that("V6: david_cluster() gains no gene_delim (REQ-006-6)", {
  expect_false("gene_delim" %in% names(formals(richCluster::david_cluster)))
})


# ===========================================================================
# The remaining hardcoded literal: the MERGE-stage dedup.
#
# R/merge_enrichment_results.R:138-139 splits on "," and re-joins on "," no
# matter what the caller's delimiter is, and cluster() calls it without passing
# gene_delim at all.  With gene_delim=";" two frames carrying "A;B;C" and
# "B;C;D" merge to "A;B;C,B;C;D" -- and the C++, splitting that on ";", yields
# a gene literally named "C,B" that appears in no input.  Dedup also silently
# fails: B and C survive twice.
#
# V4/V5 above cannot see this, because rc_t006_frames() builds two IDENTICAL
# frames -- unique() collapses them to one token before any cross-delimiter
# join can happen.  These assertions use frames that genuinely differ.
# ===========================================================================

# Two frames over the same terms whose gene sets OVERLAP but are not equal --
# the only shape in which the merge-stage join is observable.
rc_t006_distinct_frames <- function(delim = ",", n_terms = 40, block = 6, overlap = 3) {
  pool  <- sprintf("GENE%03d", seq_len(n_terms * block + 8))
  terms <- sprintf("TERM_%03d", seq_len(n_terms))
  mk <- function(shift) {
    genes <- vapply(seq_len(n_terms), function(i) {
      start <- (i - 1) * (block - overlap) + 1 + shift
      paste(pool[start:(start + block - 1)], collapse = delim)
    }, character(1))
    data.frame(Term = terms, GeneID = genes,
               Pvalue = rep(1e-6, n_terms), Padj = rep(1e-6, n_terms),
               stringsAsFactors = FALSE)
  }
  list(mk(0), mk(2))   # shift 2: overlapping but distinct gene sets
}


test_that("V7: merge_enrichment_results() exposes gene_delim, defaulting to a comma", {
  expect_true("gene_delim" %in% names(formals(richCluster::merge_enrichment_results)))
  expect_identical(formals(richCluster::merge_enrichment_results)$gene_delim, ",")
})


test_that("V8: the merge dedups at gene level under a non-comma delimiter", {
  a <- data.frame(Term = "T1", GeneID = "A;B;C", Pvalue = 1e-6, Padj = 1e-6,
                  stringsAsFactors = FALSE)
  b <- data.frame(Term = "T1", GeneID = "B;C;D", Pvalue = 1e-6, Padj = 1e-6,
                  stringsAsFactors = FALSE)
  m <- richCluster::merge_enrichment_results(list(a, b), gene_delim = ";")

  # Gene-level dedup, first-appearance order -- the comma contract, under ";".
  expect_identical(m$GeneID, "A;B;C;D")
  # The fabricated-identifier symptom: no foreign separator may survive.
  expect_false(grepl(",", m$GeneID, fixed = TRUE))
})


test_that("V9: cluster() threads gene_delim into the merge, end to end", {
  x <- rc_t006_distinct_frames(delim = ";")
  r <- rc_t006_quiet(richCluster::cluster(x, gene_delim = ";"))

  # Pre-fix this is "GENE001;...;GENE006,GENE003;...;GENE008" -- a comma the
  # user never supplied, splitting into a gene present in neither input.
  expect_false(any(grepl(",", r$merged_df$GeneID, fixed = TRUE)))

  # And the dedup must actually have deduped: each term's gene list carries no
  # repeated identifier.
  dup_free <- vapply(strsplit(r$merged_df$GeneID, ";", fixed = TRUE),
                     function(g) !anyDuplicated(g), logical(1))
  expect_true(all(dup_free))
})


test_that("V10: the comma default preserves merge behaviour byte for byte", {
  a <- data.frame(Term = "T1", GeneID = "A,B,C", Pvalue = 1e-6, Padj = 1e-6,
                  stringsAsFactors = FALSE)
  b <- data.frame(Term = "T1", GeneID = "B,C,D", Pvalue = 1e-6, Padj = 1e-6,
                  stringsAsFactors = FALSE)
  expect_identical(richCluster::merge_enrichment_results(list(a, b))$GeneID, "A,B,C,D")
  expect_identical(
    richCluster::merge_enrichment_results(list(a, b), gene_delim = ",")$GeneID, "A,B,C,D")
})


test_that("V11: an unusable gene_delim is refused, as the C++ entry already refuses it", {
  a <- data.frame(Term = "T1", GeneID = "A,B", Pvalue = 1e-6, Padj = 1e-6,
                  stringsAsFactors = FALSE)
  # "" would split every identifier into single characters rather than erroring,
  # which is the silent-corruption shape 7c1ec17 blocked at the C++ entry.
  #
  # The regexp is the guard's own wording, NOT the bare token "gene_delim":
  # before the parameter existed these calls errored with "unused argument
  # (gene_delim = ...)", and a loose match would have accepted that -- passing
  # against the unfixed build for entirely the wrong reason.
  expect_error(richCluster::merge_enrichment_results(list(a), gene_delim = ""),
               "gene_delim must be")
  expect_error(richCluster::merge_enrichment_results(list(a), gene_delim = c(",", ";")),
               "gene_delim must be")
  expect_error(richCluster::merge_enrichment_results(list(a), gene_delim = 1),
               "gene_delim must be")
})
