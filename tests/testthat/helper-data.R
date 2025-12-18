load_cluster_result <- function() {
  path <- system.file("extdata", "cluster_result.rds", package = "richCluster")
  testthat::skip_if(path == "", "Example clustering result not found.")
  readRDS(path)
}
