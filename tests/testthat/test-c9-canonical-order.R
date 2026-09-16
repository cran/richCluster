# C9 (converge ledger, session 2026-09-04): cluster output ORDER was a function of
# the C++ standard library's hash-set iteration -- libc++ (macOS) and libstdc++
# (Linux, mingw) numbered the same clusters differently, listed the same members
# in a different order, and a bit-exact Padj tie made get_representative_terms()
# pick a different term.  R CMD check on macOS ended in `Status: 1 ERROR` while
# Linux and Windows passed the identical tarball.  These tests pin the canonical
# order the fix delivers; they FAILED on the pre-fix build on Linux (RED first).
#
# Canonical term order (SPEC-RC-011, runRichCluster): (ascii_fold(Term), Term,
# GeneID), bytewise.  ascii_fold maps ONLY 'A'..'Z' to 'a'..'z'.

c9_quiet <- function(expr) { invisible(capture.output(r <- expr)); r }

c9_input <- function() {
  f1 <- read.delim(system.file("extdata", "HF36wk_vs_HF12wk.txt", package = "richCluster"),
                   stringsAsFactors = FALSE)
  f2 <- read.delim(system.file("extdata", "HF36wk_vs_WT12wk.txt", package = "richCluster"),
                   stringsAsFactors = FALSE)
  md <- c9_quiet(merge_enrichment_results(list(f1, f2)))
  md[md$Pvalue < 1e-4, ]
}

# rank of every row under the canonical term order, C-locale bytewise
c9_canonical_rank <- function(term, geneid) {
  fold <- chartr("ABCDEFGHIJKLMNOPQRSTUVWXYZ", "abcdefghijklmnopqrstuvwxyz", term)
  rank(order(order(fold, term, geneid, method = "radix")), ties.method = "first")
}

c9_members <- function(all_clusters) {
  lapply(all_clusters$TermIndices,
         function(s) as.integer(strsplit(s, ", ", fixed = TRUE)[[1]]) + 1L)
}

# The three tests below each cluster the 580-term input (~5-16 s apiece locally)
# and skip on CRAN, which caps a whole check at 10 min.
test_that("C9: members of every cluster are listed in canonical term order", {
  skip_on_cran()
  md <- c9_input()
  r  <- c9_quiet(runRichCluster(md$Term, md$GeneID, "kappa", 0.5, "average", 0.5))
  rk <- c9_canonical_rank(md$Term, md$GeneID)
  members <- c9_members(r$all_clusters)
  expect_gt(length(members), 0)
  in_order <- vapply(members, function(m) !is.unsorted(rk[m], strictly = TRUE), logical(1))
  expect_true(all(in_order),
              label = sprintf("%d of %d clusters list members out of canonical order",
                              sum(!in_order), length(in_order)))
})

test_that("C9: clusters are numbered by their first member's canonical rank", {
  skip_on_cran()
  md <- c9_input()
  r  <- c9_quiet(runRichCluster(md$Term, md$GeneID, "kappa", 0.5, "average", 0.5))
  rk <- c9_canonical_rank(md$Term, md$GeneID)
  first_rank <- vapply(c9_members(r$all_clusters), function(m) min(rk[m]), numeric(1))
  expect_false(is.unsorted(first_rank))   # non-strict: clusters may overlap on their first member
  expect_equal(r$all_clusters$Cluster, seq_along(first_rank))
})

test_that("C9: cluster numbering AND member order survive an input permutation", {
  skip_on_cran()
  md <- c9_input()
  tv <- md$Term; gv <- md$GeneID
  seq_of <- function(ord) {
    t2 <- tv[ord]
    r <- c9_quiet(runRichCluster(t2, gv[ord], "kappa", 0.5, "average", 0.5))
    lapply(c9_members(r$all_clusters), function(m) t2[m])   # ordered, NOT sorted
  }
  ref <- seq_of(seq_along(tv))
  set.seed(20260904)
  for (nm in c("reverse", "shuffle")) {
    ord <- if (nm == "reverse") rev(seq_along(tv)) else sample(seq_along(tv))
    expect_identical(seq_of(ord), ref, info = paste("permutation:", nm))
  }
})

test_that("C9: get_representative_terms breaks a value tie by the other column, then Term", {
  # two rows share Cluster and Padj bit-for-bit; the larger Pvalue is listed FIRST
  cd <- data.frame(Cluster = c(1L, 1L, 2L, 2L),
                   Term    = c("zeta", "alpha", "beta", "gamma"),
                   Padj    = c(1e-5, 1e-5, 0.5, 0.5),
                   Pvalue  = c(2e-7, 1e-7, 1e-3, 1e-3),
                   stringsAsFactors = FALSE)
  rep_padj <- unname(richCluster:::get_representative_terms(cd, "Padj"))
  expect_identical(rep_padj, c("alpha", "beta"))   # cluster 1: smaller Pvalue; cluster 2: Term
  rep_pv <- unname(richCluster:::get_representative_terms(cd, "Pvalue"))
  expect_identical(rep_pv, c("alpha", "beta"))     # cluster 1: smaller Padj? tie -> Term
})

test_that("C9: the DS-06 tie on the shipped example resolves to the smaller Pvalue", {
  f1 <- read.delim(system.file("extdata", "HF36wk_vs_HF12wk.txt", package = "richCluster"),
                   stringsAsFactors = FALSE)
  f2 <- read.delim(system.file("extdata", "HF36wk_vs_WT12wk.txt", package = "richCluster"),
                   stringsAsFactors = FALSE)
  r <- c9_quiet(cluster(list(head(f1, 120), head(f2, 120)), min_value = 0.05))
  cd <- r$cluster_df
  two <- cd[cd$Term %in% c("locomotion", "movement of cell or subcellular component"), ]
  expect_equal(nrow(two), 2L)
  expect_identical(two$Padj[1], two$Padj[2])            # the tie is bit-exact
  reps <- unname(richCluster:::get_representative_terms(cd, "Padj"))
  expect_true("movement of cell or subcellular component" %in% reps)
  expect_false("locomotion" %in% reps)
})
