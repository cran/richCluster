# NULL placeholder for roxygen namespace declarations
NULL

#' Cluster Terms using DAVID's method
#'
#' This function performs clustering on enrichment results using an algorithm
#' inspired by DAVID's functional clustering method.
#'
#' @param enrichment_results A list of dataframes, each containing enrichment results.
#'        Each dataframe should include at least the columns 'Term', 'GeneID', and 'Padj'.
#' @param df_names Optional, a character vector of names for the enrichment result dataframes. Must
#'        match the length of `enrichment_results`. Default is `NULL`.
#' @param similarity_threshold A numeric value for the kappa score cutoff (0 < cutoff <= 1).
#'        The comparison is STRICT: a pair counts as similar when its kappa is
#'        greater than the threshold, not when it equals it.
#' @param initial_group_membership Minimum number of terms to form an initial seed group.
#' @param final_group_membership Minimum number of terms for a final cluster.
#' @param multiple_linkage_threshold A numeric value for the merging threshold.
#'        The comparison is STRICT: groups merge when their shared fraction is
#'        greater than the threshold, not when it equals it.
#'
#' @param verbose Logical; print the C++ core's progress narration to the
#'        console. Default `FALSE` --- the core runs silent. Before 2.0.0 this
#'        narration was unconditional and could not be switched off, and it
#'        scales with the data (one line per merge iteration), so it is now
#'        opt-in.
#'
#' @return A named list containing the clustering results.
#'
#' @examples
#' \donttest{
#' # A small two-dataset input.  The shipped data is far larger, and an example
#' # that runs during R CMD check should stay quick.
#' terms <- sprintf("TERM_%02d", 1:12)
#' genes <- vapply(1:12, function(i) paste0("G", i:(i + 4), collapse = ","),
#'                 character(1))
#' d1 <- data.frame(Term = terms, GeneID = genes,
#'                  Pvalue = 1e-5, Padj = 1e-4, stringsAsFactors = FALSE)
#' d2 <- d1
#' res <- david_cluster(list(d1, d2), similarity_threshold = 0.5)
#' names(res)
#' }
#' @export
david_cluster <- function(enrichment_results, df_names = NULL,
                          similarity_threshold = 0.5,
                          initial_group_membership = 3,
                          final_group_membership = 3,
                          multiple_linkage_threshold = 0.5,
                          verbose = FALSE) {

  validate_david_inputs(enrichment_results, similarity_threshold,
                        initial_group_membership, final_group_membership,
                        multiple_linkage_threshold)

  if (is.null(df_names) || length(enrichment_results) != length(df_names)) {
    df_names <- as.character(seq_along(enrichment_results))
  }

  merged_df <- merge_enrichment_results(enrichment_results)
  term_vec <- merged_df$Term
  geneID_vec <- merged_df$GeneID

  cluster_result <- runDavidClustering(
    term_vec,
    geneID_vec,
    similarity_threshold,
    initial_group_membership,
    final_group_membership,
    multiple_linkage_threshold,
    verbose
  )

  cluster_options <- list(
    similarity_threshold = similarity_threshold,
    initial_group_membership = initial_group_membership,
    final_group_membership = final_group_membership,
    multiple_linkage_threshold = multiple_linkage_threshold
  )

  cluster_result$df_list <- enrichment_results
  cluster_result$merged_df <- merged_df
  cluster_result$cluster_options <- cluster_options
  cluster_result$df_names <- df_names
  cluster_result$final_clusters <- cluster_result$clusters
  cluster_result$cluster_df <- make_full_clusterdf(cluster_result$final_clusters, merged_df)

  return(cluster_result)
}

# SPEC-RC-009 / DS-04 -- mirror of validate_inputs() (R/cluster.R:139-166) for
# the DAVID path.  Domains verified against the C++ in SPEC-RC-009 spec.md §C:
# out-of-domain values previously returned zero clusters (or garbage) silently.
validate_david_inputs <- function(enrichment_results,
                                  similarity_threshold,
                                  initial_group_membership,
                                  final_group_membership,
                                  multiple_linkage_threshold) {
  if (!is.list(enrichment_results)) {
    stop("enrichment_results must be a list of dataframes.")
  }
  if (any(!sapply(enrichment_results, is.data.frame))) {
    stop("Each element of enrichment_results must be a dataframe.")
  }
  check_unit_interval <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x > 1) {
      stop(sprintf("%s must be between 0 and 1.", name))
    }
  }
  check_whole_count <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1 || x %% 1 != 0) {
      stop(sprintf("%s must be a whole number >= 1.", name))
    }
  }
  check_unit_interval(similarity_threshold, "similarity_threshold")
  check_unit_interval(multiple_linkage_threshold, "multiple_linkage_threshold")
  check_whole_count(initial_group_membership, "initial_group_membership")
  check_whole_count(final_group_membership, "final_group_membership")
}



# ===========================================================================
# T3-04a (WAVE 0) -- THE R-SIDE ACCESSOR FOR THE DAVID KAPPA MATRIX.
#
# INTERNAL, and deliberately not exported: it is an observability hook for the
# Wave-0 oracle gate, not new public API.  Reach it as
# `richCluster:::david_kappa_matrix()`, exactly as the frozen-snapshot helper
# already reaches `richCluster:::runDavidClustering()`.
#
# WHY IT EXISTS.  DavidClustering::calculateKappaScores() counts the gene
# intersection correctly with a manual membership loop, so it is the in-repo
# reference implementation the pure-R kappa oracle (T3-02) is validated
# against, and it predates every edit in this plan.  As shipped, that matrix
# was NOT OBSERVABLE FROM R BY ANY MEANS: `kappaMatrix` and
# `calculateKappaScores()` were both private members of the C++ class,
# `run()` returned only `list(clusters = ...)`, and `runDavidClustering` is
# not in NAMESPACE.  The oracle-agreement gate the plan calls its single most
# load-bearing verification was therefore not executable at all.
#
# WHAT CHANGED, AND WHAT DID NOT.  Only observability.  The computation is
# byte-for-byte the one 1.0.2 shipped: `run()`, `calculateKappaScores()`,
# `runDavidClustering()` and `david_cluster()` are all untouched, which is why
# T3-04a V1 can prove neutrality by bit-identity against the v102-david and
# v102-dc tags T3-03 froze on the stock build.  The claim this licenses is
# "agrees with the untouched COMPUTATION, observed behind an additive
# accessor" -- not "agrees with a literally untouched tree".
#
# THE DIAGONAL IS 0, NOT 1.  calculateKappaScores() loops
# `for (j = i + 1; ...)` and never writes the diagonal, so it keeps the
# constructor's 0.0 fill.  EXCLUDE THE DIAGONAL from every oracle comparison
# built on this matrix: a correct Cohen's-kappa oracle gives kappa(A, A) == 1
# there, so an inclusive comparison would fail on n_terms entries for a reason
# that has nothing to do with the intersection defect.  The diagonal is never
# read by the clustering and is NOT "fixed" here.
#
# @param terms character vector of term names.
# @param gene_ids character vector of comma-separated gene ID strings, parallel
#   to `terms`.
# @param similarity_threshold,initial_group_membership,final_group_membership,multiple_linkage_threshold
#   the four DAVID parameters, defaulted exactly as `david_cluster()` defaults
#   them.  The matrix itself depends only on `terms`/`gene_ids`; the other
#   three parameters affect only the clustering stage.
# @return a numeric n_terms x n_terms matrix, un-dimnamed, indexed in the order
#   of `terms`.  Symmetric, finite, zero on the diagonal.
#
# For the matrix TOGETHER WITH the clustering that consumed it -- which is what
# proves the accessor is not reading a stale or separately computed copy --
# call `runDavidClusteringWithKappa()` directly: it returns `clusters` (exactly
# what `runDavidClustering()` returns) plus `kappa_matrix`.
david_kappa_matrix <- function(terms, gene_ids,
                               similarity_threshold = 0.5,
                               initial_group_membership = 3,
                               final_group_membership = 3,
                               multiple_linkage_threshold = 0.5) {

  runDavidClusteringWithKappa(
    terms,
    gene_ids,
    similarity_threshold,
    initial_group_membership,
    final_group_membership,
    multiple_linkage_threshold
  )$kappa_matrix
}
