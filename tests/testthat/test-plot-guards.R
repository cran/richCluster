# SPEC-RC-011 -- cluster-id guards on the plotting layer (TDD RED first).
#
# Four defects, all on surfaces that accept a user-supplied cluster id:
#   DS-14  term_hmap() rejected numeric ids its own error message called valid.
#   DS-13  term_bar()/term_dot() recycled element-wise when a term name
#          resolved to more than one cluster.
#   DS-15  a nonexistent cluster id produced a silently empty plot.
#   DS-16  1:length() on an empty term vector iterated c(1, 0).
#
# The sibling surfaces cluster_network() / cluster_correlation_hmap() already
# carry the DS-05 (SPEC-RC-010) guard; these tests hold the rest of the
# plotting layer to the same contract.

rc011_result <- load_cluster_result()
rc011_df <- rc011_result$cluster_df

# Measured on the shipped example object (inst/extdata/cluster_result.rds,
# 2026-08-28): 53 of 298 unique terms sit in 2+ clusters, so the DS-13
# recycling path was reached by roughly one term in six, not by a corner case.
rc011_multi_term <- "anatomical structure morphogenesis"   # in 2 clusters (ids read from the object below)
rc011_single_term <- "circulatory system development"      # cluster 1 only
rc011_absent_term <- "no such term exists in this object"
rc011_absent_id <- 9999

# heatmap row count, read off the plotly layout: since C8 (2026-09-04) term_hmap()
# returns a plotly object drawn by heatmaply (iheatmapr dropped), and every row
# label sits on the y axis the matrix and the cluster strip share.
rc011_hmap_rows <- function(h) length(plotly::plotly_build(h)$x$layout$yaxis$ticktext)

test_that("DS-14: term_hmap accepts a numeric cluster id", {
  h <- term_hmap(rc011_result, clusters = 1, terms = NULL,
                 value_type = "Padj", aggr_type = mean)
  expect_s3_class(h, "plotly")
  expect_equal(rc011_hmap_rows(h), sum(rc011_df$Cluster == 1))
})

test_that("DS-14: term_hmap accepts an integer id and a multi-id vector", {
  hi <- term_hmap(rc011_result, clusters = 1L, terms = NULL,
                  value_type = "Padj", aggr_type = mean)
  expect_equal(rc011_hmap_rows(hi), sum(rc011_df$Cluster == 1))

  hv <- term_hmap(rc011_result, clusters = c(1, 2), terms = NULL,
                  value_type = "Padj", aggr_type = mean)
  expect_equal(rc011_hmap_rows(hv), sum(rc011_df$Cluster %in% c(1, 2)))
})

test_that("C8: term_hmap carries the cluster strip and honours a custom title", {
  h <- term_hmap(rc011_result, clusters = c(1, 2), title = "custom title")
  b <- plotly::plotly_build(h)
  types <- vapply(b$x$data, function(tr) tr$type, character(1))
  expect_equal(sum(types == "heatmap"), 2)   # the matrix + a one-column cluster strip
  ncols <- vapply(b$x$data[types == "heatmap"], function(tr) length(tr$x), integer(1))
  expect_setequal(ncols, c(length(rc011_result$df_names), 1L))
  # plotly stores the title as a string or as list(text = ...) depending on version
  title_of <- function(b) { tt <- b$x$layout$title; if (is.list(tt)) tt$text else tt }
  expect_identical(title_of(b), "custom title")
  h2 <- term_hmap(rc011_result, clusters = 1)
  expect_identical(title_of(plotly::plotly_build(h2)), "-log10(Padj)")
})

test_that("C8: a term in several clusters gets one uniquely labelled row per cluster, without a warning", {
  # rc011_multi_term sits in two clusters; with no cluster filter both copies are drawn
  expect_no_warning(h <- term_hmap(rc011_result))
  labs <- plotly::plotly_build(h)$x$layout$yaxis$ticktext
  ids <- sort(unique(rc011_df$Cluster[rc011_df$Term == rc011_multi_term]))
  expect_setequal(grep(rc011_multi_term, labs, fixed = TRUE, value = TRUE),
                  sprintf("%s (cluster %s)", rc011_multi_term, ids))
  expect_false(any(duplicated(labs)))
})

test_that("C11: clusters sharing a representative term get uniquely labelled rows, without a warning", {
  # cluster_hmap() labels one row per cluster with that cluster's
  # representative term -- the minimum merged value_type term.  Two clusters
  # can share one, and on the shipped object two pairs do, so heatmaply was
  # handed duplicate row names: it prefixed them with numbers and warned, and
  # that warning is visible in the vignette.  Same defect the C8 fix closed on
  # term_hmap(); cluster_hmap() was not covered by it.
  reps <- unname(richCluster:::get_representative_terms(rc011_df, "Padj"))
  dup_reps <- unique(reps[duplicated(reps)])
  expect_gt(length(dup_reps), 0)          # the object still exercises this path

  expect_no_warning(h <- cluster_hmap(rc011_result))
  labs <- plotly::plotly_build(h)$x$layout$yaxis$ticktext
  expect_false(any(duplicated(labs)))
  expect_length(labs, length(reps))

  ids <- sort(unique(rc011_df$Cluster))
  for (tm in dup_reps) {
    expect_true(all(sprintf("%s (cluster %s)", tm, ids[reps == tm]) %in% labs))
  }
  # a representative term unique to one cluster keeps its bare name
  solo <- setdiff(reps, dup_reps)
  expect_true(all(solo %in% labs))
})

test_that("DS-14: term_hmap's documented defaults are its actual defaults", {
  # The roxygen block called clusters/terms/value_type/aggr_type optional while
  # the signature declared no defaults, so term_hmap(cluster_result) failed on
  # a missing argument.
  h <- term_hmap(rc011_result)
  expect_equal(rc011_hmap_rows(h), nrow(rc011_df))
})

test_that("DS-15: term_hmap rejects a nonexistent numeric cluster id", {
  expect_error(
    term_hmap(rc011_result, clusters = rc011_absent_id, terms = NULL,
              value_type = "Padj", aggr_type = mean),
    "does not match any cluster id")
})

test_that("DS-13: term_bar/term_dot refuse a term that is in several clusters", {
  # C9 (2026-09-04): the two ids are read from the object rather than pinned --
  # cluster numbering is canonical now and the shipped object was regenerated
  # (the pinned "25, 30" of the pre-C9 object became "3, 5").
  ids <- sort(unique(rc011_df$Cluster[rc011_df$Term == rc011_multi_term]))
  expect_length(ids, 2)
  ids_txt <- paste(ids, collapse = ", ")
  expect_error(term_bar(rc011_result, cluster = rc011_multi_term),
               "is in 2 clusters")
  expect_error(term_bar(rc011_result, cluster = rc011_multi_term), ids_txt, fixed = TRUE)
  expect_error(term_dot(rc011_result, cluster = rc011_multi_term),
               "is in 2 clusters")
  expect_error(term_dot(rc011_result, cluster = rc011_multi_term), ids_txt, fixed = TRUE)
})

test_that("DS-13: term_bar/term_dot refuse a term that is in no cluster", {
  expect_error(term_bar(rc011_result, cluster = rc011_absent_term),
               "is not in any cluster")
  expect_error(term_dot(rc011_result, cluster = rc011_absent_term),
               "is not in any cluster")
})

test_that("DS-13: a term in exactly one cluster still resolves to that cluster", {
  k <- unique(rc011_df$Cluster[rc011_df$Term == rc011_single_term])
  expect_length(k, 1)
  by_term <- plotly::plotly_build(term_bar(rc011_result, cluster = rc011_single_term))
  by_id   <- plotly::plotly_build(term_bar(rc011_result, cluster = k))
  expect_identical(as.character(by_term$x$data[[1]]$y),
                   as.character(by_id$x$data[[1]]$y))
  expect_identical(by_term$x$layout$title, by_id$x$layout$title)
})

test_that("DS-13: a multi-id numeric cluster is refused, not recycled", {
  expect_error(term_bar(rc011_result, cluster = c(1, 2)), "single cluster")
  expect_error(term_dot(rc011_result, cluster = c(1, 2)), "single cluster")
})

test_that("DS-15: term_bar/term_dot reject a nonexistent cluster id", {
  expect_error(term_bar(rc011_result, cluster = rc011_absent_id),
               "does not match any cluster id")
  expect_error(term_dot(rc011_result, cluster = rc011_absent_id),
               "does not match any cluster id")
})

test_that("DS-15: compare_network_graphs_plotly rejects a nonexistent cluster id", {
  expect_error(
    compare_network_graphs_plotly(rc011_result, rc011_absent_id,
                                  c("Padj_1", "Padj_2")),
    "does not match any cluster id")
})

test_that("DS-16: plot_network_graph rejects a nonexistent cluster id", {
  # Before the guard this reached `node_colors[[0]] <- ...` via 1:length(x) on
  # an empty vector and died with "attempt to select less than one element".
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(
    plot_network_graph(rc011_result, rc011_absent_id,
                       rc011_result$distance_matrix, c("Padj_1", "Padj_2")),
    "does not match any cluster id")
})

test_that("DS-15/DS-16: the valid-id paths still draw", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(
    plot_network_graph(rc011_result, 1, rc011_result$distance_matrix,
                       c("Padj_1", "Padj_2")))
  expect_s3_class(
    compare_network_graphs_plotly(rc011_result, 1, c("Padj_1", "Padj_2")),
    "plotly")
  expect_s3_class(term_bar(rc011_result, cluster = 1), "plotly")
  expect_s3_class(term_dot(rc011_result, cluster = 1), "plotly")
})
