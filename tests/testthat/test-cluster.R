test_that("cluster handles ward linkage", {
  cluster_result <- load_cluster_result()
  result <- cluster(
    cluster_result$df_list,
    df_names = cluster_result$df_names,
    min_terms = 3,
    min_value = 0.0001,
    distance_metric = "kappa",
    distance_cutoff = 0.5,
    linkage_method = "ward",
    linkage_cutoff = 0.5
  )
  expect_true(is.data.frame(result$final_clusters))
  expect_gt(nrow(result$final_clusters), 0)
})

test_that("cluster_correlation_hmap returns heatmaply object", {
  cluster_result <- load_cluster_result()
  h <- cluster_correlation_hmap(
    cluster_result$final_clusters,
    cluster_result$distance_matrix,
    1,
    cluster_result$merged_df
  )
  expect_true("plotly" %in% class(h))
})

test_that("cluster_network returns htmlwidget", {
  cluster_result <- load_cluster_result()
  n <- cluster_network(
    cluster_result$final_clusters,
    cluster_result$distance_matrix,
    1,
    cluster_result$merged_df
  )
  expect_true("htmlwidget" %in% class(n))
})
