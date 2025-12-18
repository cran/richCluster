# heatmap creation
#' @importFrom magrittr %>%
#' @importFrom dplyr across bind_rows distinct filter group_by mutate pull select starts_with summarise ungroup where
NULL

# utils
na_to_zero <- function(x) {
  # returns the column vector but replaces all places where it's nan/na/inf
  # to be zero (for hmap)
  x[is.nan(x) | is.na(x) | is.infinite(x)] <- 0
  return(x)
}

# create a heatmap with all clusters

#' Create a Heatmap of Clustered Enrichment Results
#'
#' Generates an interactive heatmap from the given clustering results,
#' visualizing -log10(Padj) values for each cluster. The function aggregates
#' values per cluster and assigns representative terms as row names.
#'
#' @param cluster_result A list containing a data frame (`cluster_df`) with clustering results.
#'   The data frame must contain at least the columns `Cluster`, `Term`, and `value_type_*` values.
#' @param clusters Optional. A numeric or character vector specifying the clusters to include.
#'   If NULL (default), all clusters are included.
#' @param value_type A character string specifying the column name prefix for values to display in hmap cells.
#'   Defaults to `"Padj"`.
#' @param aggr_type A function used to aggregate values across clusters (e.g., `mean` or `median`).
#'   Defaults to `mean`.
#'
#' @return An interactive heatmap object (`plotly`), displaying the -log10(Padj) values
#'   across clusters, with representative terms as row labels.
#'
#' @details
#' The function processes the given cluster data frame (`cluster_df`),
#' aggregating the `value_type_*` values per cluster using the specified `aggr_type` function.
#' The -log10 transformation is applied, and infinite values are replaced with 0.
#'
#' Representative terms are selected by choosing the term with the lowest
#' `value_type` in each cluster.
#'
#' The final heatmap is generated using `heatmaply::heatmaply()`, with
#' an interactive `plotly` visualization.
#'
#' @export
cluster_hmap <- function(cluster_result, clusters=NULL, value_type="Padj", aggr_type=mean){
  # the hmap processing flow
  # for full_hmap
  cluster_df <- cluster_result$cluster_df

  # hmap_matrix
  hmap_matrix <- cluster_df %>%
    group_by(Cluster) %>%
    summarise(across(starts_with(paste0(value_type, "_")), function(x) mean(x, na.rm=TRUE))) %>%
    # mutate(across(where(is.numeric), na_to_zero)) %>%
    select(-Cluster) %>% # remove (don't wanna plot her)
    mutate(across(where(is.numeric), function(x) -log10(x))) %>%
    mutate(across(where(is.numeric), function(x) ifelse(is.infinite(x), 0, x))) %>%
    as.matrix()

  # representative term making
  # use the Pvalue/Padj average column
  representative_terms <- cluster_df %>%
    group_by(Cluster) %>%
    filter(value_type==min(value_type, na.rm=TRUE)) %>%
    slice(1) %>%
    ungroup %>%
    pull(Term)
  rownames(hmap_matrix) <- representative_terms
  colnames(hmap_matrix) <- cluster_result$df_names

  # the hmap object
  hmap <- heatmaply::heatmaply(
    hmap_matrix,
    xlab = "Enrichment Result",
    ylab = "Cluster",
    main = paste0("-log10(", value_type, ")"),
    colors = viridis::viridis(256),
    na.value = "grey",
    margins = c(50, 50, 50, 50),
    cluster_rows=FALSE, cluster_cols=FALSE,
    Rowv=FALSE, Colv=FALSE,
    plot_method = "plotly",
    colorbar_title = paste0("-log10(", value_type, ")")
  )
  return(hmap)
}


# Heatmap displaying all terms in the specified clusters
# optionally accepts explicit list of terms to visualize if specified
# and we display the union of the two clusters/terms vectors in final result
#
# --- Representative Term Handling ---
# Clusters can be specified by cluster #
# or by the name of any term in the cluster

#' Generate a Heatmap of Enrichment Results for Specific Clusters and Terms
#'
#' Creates an interactive heatmap displaying -log10(Padj) values for selected clusters
#' and terms. Users can specify clusters numerically or select them by providing term names.
#' The function ensures that the final heatmap includes all terms from the selected clusters
#' as well as any explicitly provided terms.
#'
#' @param cluster_result A list containing a data frame (`cluster_df`) with clustering results.
#'   The data frame must include at least the columns `Cluster`, `Term`, and `Padj_*` values.
#' @param clusters Optional. A numeric vector specifying the cluster numbers to display,
#'   or a character vector specifying terms whose clusters should be included. Defaults to `NULL`,
#'   which includes all clusters.
#' @param terms Optional. A character vector specifying additional terms to include in the heatmap.
#'   Defaults to `NULL`.
#' @param value_type A character string specifying the column name prefix for adjusted p-values.
#'   Defaults to `"Padj"`.
#' @param aggr_type A function used to aggregate values across clusters (e.g., `mean` or `median`).
#'   Defaults to `mean`.
#' @param title An optional parameter to title the plot something else.
#'
#' @return An interactive heatmap object (`plotly`), displaying the -log10(Padj) values
#'   across clusters, with representative terms as row labels and color-coded cluster annotations.
#'
#' @details
#' The function processes the given `cluster_df`, identifying the clusters and terms to be visualized.
#' If `clusters` is specified as a numeric vector, the function directly filters based on cluster numbers.
#' If `clusters` is given as a character vector, it identifies the clusters associated with those terms
#' and retrieves all terms from the selected clusters.
#'
#' The `Padj_*` values are transformed using `-log10()`, and infinite values are replaced with `0`.
#' The resulting heatmap is generated using `heatmaply::heatmaply()` with fixed row ordering
#' (no hierarchical clustering).
#'
#' @export
term_hmap <- function(cluster_result, clusters, terms, value_type, aggr_type, title=NULL) {

  cluster_df <- cluster_result$cluster_df

  # get numeric vector of cluster numbers to display
  if (is.null(clusters)) {
    # no clusters specified -> by default show all
    clusters <- unique(cluster_df$Cluster)
  }  else if (is.character(clusters)) {
    # --- representative term search ---
    # user supplied terms -> get cluster numbers
    clusters <- cluster_df %>%
      filter(Term %in% clusters) %>%
      pull(Cluster) %>%
      unique()
  } else {
    stop("`clusters` must be either numeric (cluster #s) or character (term names).")
  }
  # use cluster numbers to get all terms in specified clusters
  cluster_terms <- cluster_df %>%
    filter(Cluster %in% clusters) %>%
    select(Cluster, Term, starts_with(paste0(value_type, "_")))

  # search for specific terms if supplied
  if (is.null(terms)) {
    # skip the following checks
    specific_terms <- c()
  } else if (!is.character(terms)) {
    stop("`terms` must be character (Term names).")
  } else {
    specific_terms <- cluster_df %>%
      filter(Term %in% terms) %>%
      select(Cluster, Term, starts_with(paste0(value_type, "_")))
  }

  # get the UNION of all terms in specified clusters
  # and those specified by terms
  final_terms <- bind_rows(specific_terms, cluster_terms) %>%
    distinct()

  # create the hmap_matrix
  # performing final value updates
  hmap_matrix <- final_terms %>%
    group_by(Cluster) %>%
    # mutate(across(where(is.numeric), na_to_zero)) %>%
    # select(-Cluster) %>% # remove (don't wanna plot her)
    mutate(across(where(is.numeric), function(x) -log10(x))) %>%
    mutate(across(where(is.numeric), function(x) ifelse(is.infinite(x), 0, x)))

  # keep these vars for labeling
  cluster_annots <- hmap_matrix$Cluster
  row_names <- hmap_matrix$Term

  hmap_matrix <- hmap_matrix %>%
    ungroup() %>%
    select(starts_with(paste0(value_type, "_"))) %>%
    as.matrix()
  rownames(hmap_matrix) <- row_names
  colnames(hmap_matrix) <- cluster_result$df_names

  # generate default title if none supplied
  if (is.null(title)) {
    cluster_str <- paste(final_terms, ', ')
    title <- paste0("-log10(", value_type, ")")
  }

  h <- iheatmapr::main_heatmap(
    hmap_matrix,
    name=paste0("-log10(", value_type, ")")
  ) %>%
    iheatmapr::add_row_title("Term") %>%
    iheatmapr::add_col_title(title, side=c("top")) %>%
    iheatmapr::add_col_title("Enrichment Result", side=c("bottom")) %>%
    iheatmapr::add_row_annotation(data.frame("Cluster"=cluster_annots))


  # h <- heatmaply::heatmaply(
  #   hmap_matrix,
  #   xlab = "Enrichment Result",
  #   ylab = "Term",
  #   main = paste0("-log10(", value_type, ")"),
  #   colors = viridis::viridis(256),
  #   row_side_colors = cluster_annots,
  #   row_text_angle = 0,
  #   margins = c(60, 120, 40, 10),
  #   plot_method = "plotly",
  #   colorbar_title = paste0("-log10(", value_type, ")"),
  #   cluster_rows=FALSE, cluster_cols=FALSE,
  #   Rowv=FALSE, Colv=FALSE,
  # )
  return(h)
}
