# tests/testthat/test-david-validation.R
# SPEC-RC-009 / DS-04: david_cluster() validates its inputs the way cluster()
# always has.  Synthetic frames; no file IO, no fixtures.
rc_ds04_frames <- function() {
  d1 <- data.frame(Term = c("T1", "T2", "T3", "T4"),
                   GeneID = c("a,b,c,d", "a,b,c,e", "a,b,f,g", "x,y,z"),
                   Pvalue = 1e-5, Padj = 1e-4, stringsAsFactors = FALSE)
  d2 <- data.frame(Term = c("T1", "T2", "T5", "T6"),
                   GeneID = c("a,b,c,d", "a,b,c,f", "p,q,r", "p,q,s"),
                   Pvalue = 1e-5, Padj = 1e-4, stringsAsFactors = FALSE)
  list(d1, d2)
}

test_that("DS-04: out-of-domain parameters error instead of silently returning zero clusters", {
  fr <- rc_ds04_frames()
  msg_st  <- "similarity_threshold must be between 0 and 1."
  msg_mlt <- "multiple_linkage_threshold must be between 0 and 1."
  msg_igm <- "initial_group_membership must be a whole number >= 1."
  msg_fgm <- "final_group_membership must be a whole number >= 1."
  expect_error(david_cluster(fr, similarity_threshold = 5),  msg_st,  fixed = TRUE)
  expect_error(david_cluster(fr, similarity_threshold = 0),  msg_st,  fixed = TRUE)
  expect_error(david_cluster(fr, similarity_threshold = -5), msg_st,  fixed = TRUE)
  expect_error(david_cluster(fr, initial_group_membership = 0),   msg_igm, fixed = TRUE)
  expect_error(david_cluster(fr, initial_group_membership = -1),  msg_igm, fixed = TRUE)
  expect_error(david_cluster(fr, initial_group_membership = 2.5), msg_igm, fixed = TRUE)
  expect_error(david_cluster(fr, final_group_membership = 0),   msg_fgm, fixed = TRUE)
  expect_error(david_cluster(fr, final_group_membership = -1),  msg_fgm, fixed = TRUE)
  expect_error(david_cluster(fr, multiple_linkage_threshold = 5),  msg_mlt, fixed = TRUE)
  expect_error(david_cluster(fr, multiple_linkage_threshold = 0),  msg_mlt, fixed = TRUE)
  expect_error(david_cluster(fr, multiple_linkage_threshold = -1), msg_mlt, fixed = TRUE)
})

test_that("DS-04: non-list and non-dataframe inputs error with cluster()'s messages", {
  expect_error(david_cluster("not a list"),
               "enrichment_results must be a list of dataframes.", fixed = TRUE)
  d1 <- data.frame(Term = "T1", GeneID = "a", Pvalue = 1, Padj = 1,
                   stringsAsFactors = FALSE)
  expect_error(david_cluster(list(d1, "x")),
               "Each element of enrichment_results must be a dataframe.", fixed = TRUE)
})

test_that("DS-04: valid inputs still run at the documented defaults", {
  fr <- rc_ds04_frames()
  res <- NULL
  utils::capture.output(res <- david_cluster(fr))
  expect_true(is.data.frame(res$clusters))
  expect_identical(res$cluster_options$similarity_threshold, 0.5)
})

test_that("DS-04: the internal C++ entry points reject out-of-domain parameters", {
  expect_error(
    richCluster:::runDavidClustering(c("A", "B"), c("a,b", "a,c"), 5, 3L, 3L, 0.5),
    "similarity_threshold must be between 0 and 1.", fixed = TRUE)
  expect_error(
    richCluster:::runDavidClusteringWithKappa(c("A", "B"), c("a,b", "a,c"), 0.5, 0L, 3L, 0.5),
    "initial_group_membership must be a whole number >= 1.", fixed = TRUE)
})
