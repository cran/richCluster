#' @useDynLib richCluster, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats na.omit
NULL

utils::globalVariables(c(
  "Cluster", "ClusterName", "GeneID", "JS", "Pvalue",
  "Term", "TermIndices", "merged_richsets", "n_terms", "row_id"
))
