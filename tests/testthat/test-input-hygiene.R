# SPEC-RC-008 (DS-01 / DS-02 / DS-07) -- end-to-end input-hygiene tests.
# All inputs are inline; no fixture dependency, so this file ships in the
# tarball.  Do NOT rename to test-src-*: that prefix is tarball-excluded.

test_that("DS-01: gene-less terms never cluster together (richCluster path)", {
  d1 <- data.frame(Term = c("T1","T2","T3","T4","T5"),
                   GeneID = c("a,b,c,d","a,b,c,e","","",""),
                   Pvalue = 1e-5, Padj = 1e-4, stringsAsFactors = FALSE)
  d2 <- data.frame(Term = c("T1","T2","T3","T4","T5"),
                   GeneID = c("a,b,c,d","a,b,c,f","x,y,z","",""),
                   Pvalue = 1e-5, Padj = 1e-4, stringsAsFactors = FALSE)
  sink(nullfile())
  res <- cluster(list(d1, d2), min_value = 1, distance_metric = "jaccard",
                 distance_cutoff = 0.5, linkage_method = "average",
                 linkage_cutoff = 0.5)
  sink()
  m <- res$distance_matrix
  expect_equal(m["T4", "T5"], 0)   # stock scored the two gene-less terms 1.0
  expect_equal(m["T3", "T4"], 0)   # stock scored 0.25 on the shared phantom
  both <- vapply(res$all_clusters$TermNames,
                 function(s) grepl("T4", s) && grepl("T5", s), logical(1))
  expect_false(any(both))          # stock emitted cluster {T4, T5}
})

test_that("DS-02: comma-space gene lists parse to the intended sets", {
  j <- function(genes) {
    sink(nullfile())
    o <- runRichCluster(c("A","B"), genes, "jaccard", 0.5, "average", 0.5)
    sink()
    m <- o$distance_matrix; dimnames(m) <- NULL; m[1, 2]
  }
  expect_equal(j(c("g1, g2, g3", "g1,g2,g3")), 1)    # stock: 0.2
  expect_equal(j(c("g1, g2, g3", "g3, g2, g1")), 1)  # stock: 0.2
  expect_equal(j(c("g1,g2,", "g3,g4,")), 0)          # stock: 0.2
})

test_that("DS-02 boundary: interior whitespace is preserved (trim is ends-only)", {
  sink(nullfile())
  o <- runRichCluster(c("A","B"), c("a b,c", "a b,d"), "jaccard", 0.5,
                      "average", 0.5)
  sink()
  m <- o$distance_matrix; dimnames(m) <- NULL
  expect_equal(m[1, 2], 1/3)   # {"a b","c"} vs {"a b","d"} share "a b"
})

test_that("DS-07: unknown linkage method errors; the four valid methods do not", {
  t <- c("A","B","C","D"); g <- c("1,2,3,4","1,2,3,5","1,2,3,6","7,8,9,10")
  ok <- function(link) {
    sink(nullfile())
    o <- runRichCluster(t, g, "kappa", 0.5, link, 0.5)
    sink(); o
  }
  for (mth in c("single", "complete", "average", "ward"))
    expect_equal(nrow(ok(mth)$all_clusters), 2, info = mth)
  bad <- function(link) {
    sink(nullfile()); on.exit(sink(), add = TRUE)
    runRichCluster(t, g, "kappa", 0.5, link, 0.5)
  }
  expect_error(bad("Average"), "unsupported linkage method")
  expect_error(bad("ward.D2"), "unsupported linkage method")
  expect_error(bad(""),        "unsupported linkage method")
})

test_that("DS-01: DAVID path scores gene-less pairs at 0 and never clusters them", {
  t <- c("T1","T2","T3","T4"); g <- c("a,b,c,d","a,b,c,e","","")
  sink(nullfile())
  r <- richCluster:::runDavidClusteringWithKappa(t, g, 0.5, 2, 2, 0.5)
  sink()
  expect_equal(r[["kappa_matrix"]][3, 4], 0)   # stock: 1.0
  expect_equal(nrow(r[["clusters"]]), 0)       # stock: cluster {T3, T4}
})

test_that("DS-01 edge: all-empty input runs without error on both paths", {
  sink(nullfile())
  o <- runRichCluster(c("A","B"), c("",""), "kappa", 0.5, "average", 0.5)
  sink()
  expect_equal(o$distance_matrix[1, 2], 0)
  expect_equal(nrow(o$all_clusters), 2)
  sink(nullfile())
  r <- richCluster:::runDavidClusteringWithKappa(c("A","B"), c("",""),
                                                 0.5, 2, 2, 0.5)
  sink()
  expect_equal(nrow(r[["clusters"]]), 0)
})
