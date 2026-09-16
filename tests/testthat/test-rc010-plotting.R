# SPEC-RC-010 — DS-05 / DS-06 / DS-09 defect inversions (TDD RED first).
# One fresh cluster result from the shipped example inputs; no fixtures needed.

rc010_result <- local({
  f1 <- read.delim(system.file("extdata", "HF36wk_vs_HF12wk.txt", package = "richCluster"),
                   stringsAsFactors = FALSE)
  f2 <- read.delim(system.file("extdata", "HF36wk_vs_WT12wk.txt", package = "richCluster"),
                   stringsAsFactors = FALSE)
  # capture.output (not sink/on.exit) so the sink is restored even if cluster() errors
  invisible(capture.output(
    r <- cluster(list(head(f1, 120), head(f2, 120)), min_value = 0.05)))
  r
})

# The 13 intended labels, in cluster order. The SET is exactly the one recorded
# in spec.md §B.2 (re-confirmed 2026-09-04). The ORDER was re-measured on the
# C9 build (2026-09-04): clusters are now numbered in canonical term order --
# by content, identically on every platform -- where before they followed the
# C++ standard library's hash-set iteration (converge ledger, session
# 2026-09-04: macOS numbered them differently and failed this file).
rc010_intended <- c(
  "localization", "single-multicellular organism process",
  "regulation of cell death", "metabolic process",
  "small molecule biosynthetic process",
  "movement of cell or subcellular component",
  "cell projection organization", "protein phosphorylation",
  "response to chemical", "sterol metabolic process",
  "regulation of response to stimulus",
  "positive regulation of biological process",
  "positive regulation of developmental process")

test_that("DS-05: final_clusters ids are canonical 1..k", {
  fc <- rc010_result$final_clusters
  expect_equal(as.integer(fc$Cluster), seq_len(nrow(fc)))
})

test_that("DS-05: the same id names the same cluster on every surface", {
  fc <- rc010_result$final_clusters
  cd <- rc010_result$cluster_df
  for (k in fc$Cluster) {
    n_fc <- length(strsplit(fc$TermIndices[fc$Cluster == k], ", ")[[1]])
    expect_equal(sum(cd$Cluster == k), n_fc, label = paste("id", k))
  }
})

test_that("DS-05: cluster_network / cluster_correlation_hmap select by id, not row position", {
  fc_sparse <- rc010_result$final_clusters
  fc_sparse$Cluster <- fc_sparse$Cluster * 10L   # simulate a legacy sparse-id object
  k <- fc_sparse$Cluster[4]
  n <- cluster_network(fc_sparse, rc010_result$distance_matrix, k, rc010_result$merged_df)
  expect_equal(nrow(n$x$nodes),
               length(strsplit(fc_sparse$TermIndices[4], ", ")[[1]]))
})

test_that("DS-05: a nonexistent cluster id errors clearly", {
  expect_error(
    cluster_network(rc010_result$final_clusters, rc010_result$distance_matrix,
                    999, rc010_result$merged_df),
    "does not match any cluster id")
  expect_error(
    cluster_correlation_hmap(rc010_result$final_clusters, rc010_result$distance_matrix,
                             999, rc010_result$merged_df),
    "does not match any cluster id")
})

test_that("DS-06: get_representative_terms selects the minimum merged value_type", {
  rep_terms <- unname(richCluster:::get_representative_terms(rc010_result$cluster_df, "Padj"))
  expect_equal(rep_terms, rc010_intended)
})

test_that("DS-06: get_representative_terms errors on an absent column", {
  expect_error(
    richCluster:::get_representative_terms(rc010_result$cluster_df, "Qvalue"),
    "not a column of cluster_df")
})

test_that("DS-06: the five plot surfaces carry the intended labels", {
  pb <- plotly::plotly_build(cluster_bar(rc010_result))
  expect_true(all(as.character(pb$x$data[[1]]$y) %in% rc010_intended))
  pd <- plotly::plotly_build(cluster_dot(rc010_result))
  expect_true(all(as.character(pd$x$data[[1]]$y) %in% rc010_intended))
  ph <- plotly::plotly_build(cluster_hmap(rc010_result))
  expect_setequal(as.character(ph$x$layout$yaxis$ticktext), rc010_intended)
  # C9: look the cluster id up from the intended labels (position == id) rather
  # than hard-coding 1 -- the numbering is canonical, not "first merged".
  k_mp <- which(rc010_intended == "metabolic process")
  expect_length(k_mp, 1L)
  ptb <- plotly::plotly_build(term_bar(rc010_result, cluster = k_mp))
  expect_identical(ptb$x$layout$title, sprintf("Padj, metabolic process (cluster %d)", k_mp))
  ptd <- plotly::plotly_build(term_dot(rc010_result, cluster = k_mp))
  expect_identical(ptd$x$layout$title, sprintf("Padj, metabolic process (cluster %d)", k_mp))
})

test_that("DS-09: full_network carries no self-links, no negative values, exact edge count", {
  fn <- full_network(rc010_result)
  links <- fn$x$links
  expect_equal(sum(links$source == links$target), 0)
  expect_gte(min(links$value), 0)
  dm <- rc010_result$distance_matrix
  expect_equal(nrow(links), sum(dm[upper.tri(dm)] > 0))
  expect_equal(nrow(fn$x$nodes), nrow(dm))
})
