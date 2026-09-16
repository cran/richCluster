#' @importFrom magrittr %>%
#' @importFrom dplyr across filter first group_by mutate n pull slice starts_with summarise ungroup
NULL

#' Cluster-level Dot Plot of Enrichment Significance
#'
#' Creates a dot plot summarizing cluster-level enrichment across datasets.
#' Each point represents a cluster, with its size proportional to the number
#' of terms and its x-position reflecting average significance (e.g., Padj or Pvalue).
#'
#' @param cluster_result A result list returned from \code{\link{cluster}}.
#' @param clusters Optional numeric vector of cluster IDs to include. Defaults to all clusters.
#' @param value_type The name of the value column to visualize (e.g., "Padj" or "Pvalue").
#' @param title Optional title for the plot. If NULL, a default title is generated.
#'
#' @return A \code{plotly} object representing the dot plot.
#' @examples
#' \donttest{
#' # Load example data
#' cluster_result <- readRDS(system.file("extdata", "cluster_result.rds",
#'                                       package = "richCluster"))
#' cdot <- cluster_dot(cluster_result)
#' cdot
#' }
#' @export
cluster_dot <- function(cluster_result, clusters=NULL, value_type="Padj", title=NULL) {
  cluster_df <- cluster_result$cluster_df
  df_names <- cluster_result$df_names

  if (is.null(clusters)) {
    clusters <- unique(cluster_df$Cluster)
  }

  # dot data
  dot_data <- cluster_df %>%
    group_by(Cluster) %>%
    mutate(n_terms = n()) %>% # get n_terms count
    summarise(
      n_terms = first(n_terms), # and assign each cluster with n_terms
      across(starts_with(paste0(value_type, "_")), function(x) mean(x, na.rm=TRUE))) %>%
    mutate(across(starts_with(paste0(value_type, "_")), function(x) -log10(x))) %>%
    mutate(across(starts_with(paste0(value_type, "_")), function(x) ifelse(is.infinite(x), 0, x)))

  # representative term making: minimum of the merged value_type column.
  # DS-06 (SPEC-RC-010): the old dplyr filter compared the argument string to
  # itself and was a no-op.
  representative_terms <- unname(get_representative_terms(cluster_df, value_type))

  dot_data$ClusterName <- representative_terms

  # generate default title if none supplied
  if (is.null(title)) {
    title <- paste0(value_type, ", all clusters")
  }

  # get all dot data across different dfs
  value_cols <- grep(paste0("^", value_type, "_"), names(dot_data), value = TRUE)

  # the plotly object
  dot <- plotly::plot_ly(
    dot_data,
    x = dot_data[[value_cols[1]]],
    y = ~ClusterName,
    type = 'scatter',
    mode = 'markers',
    marker = list(opacity = 0.5),
    size = ~n_terms,
    name = df_names[1],
    text = ~paste0(
      'Cluster ', Cluster, '<br>',
      ClusterName, '<br>',
      '-log10(',value_type,') = ', sprintf("%.2f", dot_data[[value_cols[1]]]), '<br>',
      '# terms = ', n_terms
    ),
    hovertemplate = "%{text}"
  )
  if (length(value_cols) > 1) {
    for (i in 2:length(value_cols)) {
      dot <- dot %>% plotly::add_trace(
        x = dot_data[[value_cols[i]]],
        type = 'scatter',
        marker = list(opacity = 0.5),
        size = ~n_terms,
        name = df_names[i],
        text = ~paste0(
          'Cluster ', Cluster, '<br>',
          ClusterName, '<br>',
          '-log10(',value_type,') = ', sprintf("%.2f", dot_data[[value_cols[i]]]), '<br>',
          '# terms = ', n_terms
        ),
        hovertemplate = "%{text}"
      )
    }
  }
  # add layout stuff
  dot <- dot%>% plotly::layout(
    title = title,
    margin = list(autoexpand=TRUE, t=50, b=-50),
    xaxis = list(title = "-log10(Padj)"),
    yaxis = list(title = "Cluster")
  )
  return(dot)
}
# cdot <- cluster_dot(cluster_result)
# cdot

#' Term-level Dot Plot for a Specific Cluster
#'
#' Creates a dot plot of individual terms within a specified cluster, showing
#' their significance and number of genes.
#'
#' @param cluster_result A result list returned from \code{\link{cluster}}.
#' @param cluster A single cluster ID (numeric) or term name (character) to plot.
#'   A numeric ID must appear in `cluster_df$Cluster`. A term name must belong to
#'   exactly one cluster; a term shared by several clusters is ambiguous here and is
#'   rejected, naming the clusters it belongs to so one can be passed instead.
#' @param value_type The name of the value column to visualize (e.g., "Padj" or "Pvalue").
#' @param title Optional title for the plot. If NULL, a default title is generated using the representative term.
#'
#' @return A \code{plotly} object representing the dot plot of terms.
#' @examples
#' \donttest{
#' # Load example data
#' cluster_result <- readRDS(system.file("extdata", "cluster_result.rds",
#'                                       package = "richCluster"))
#' tdot <- term_dot(cluster_result, cluster = 1)
#' tdot
#' }
#' @export
term_dot <- function(cluster_result, cluster=1, value_type="Padj", title=NULL) {

  cluster_df <- cluster_result$cluster_df
  df_names <- cluster_result$df_names

  # get cluster # from the term
  if (is.character(cluster)) {
    # DS-13 (SPEC-RC-011): $Cluster on a multi-row subset returns a VECTOR, and
    # the filter() below then recycled element-wise instead of selecting one
    # cluster.  Measured on the shipped example object (2026-08-28): 53 of 298
    # terms sit in 2+ clusters -- "anatomical structure morphogenesis" is in
    # both 25 and 30 -- so this was reached by roughly one term in six.  The
    # function plots a single cluster, so an ambiguous term is refused by name.
    matched <- unique(cluster_df$Cluster[cluster_df$Term == cluster])
    if (length(matched) == 0) {
      stop(sprintf("term \"%s\" is not in any cluster", cluster))
    }
    if (length(matched) > 1) {
      stop(sprintf(
        "term \"%s\" is in %d clusters (%s); pass one of those cluster numbers as `cluster` instead",
        cluster, length(matched), paste(sort(matched), collapse = ", ")))
    }
    cluster <- matched
  } else if (!is.numeric(cluster)) {
    stop("cluster must be numeric (a cluster number) or character (a term name)")
  }

  # A length > 1 id recycles in filter() exactly as the term path used to.
  if (length(cluster) != 1) {
    stop("`cluster` must be a single cluster number or term name")
  }
  # DS-15 (SPEC-RC-011): an unknown id used to draw a silently empty dot plot.
  # The sibling surfaces cluster_network() / cluster_correlation_hmap() already
  # guard this (DS-05, SPEC-RC-010).
  validate_cluster_ids(cluster, cluster_df$Cluster, "cluster")

  # dot data
  dot_data <- cluster_df %>%
    group_by(Cluster) %>%
    filter(Cluster==cluster) %>%
    # summarise(across(starts_with(paste0(value_type, "_")), function(x) mean(x, na.rm=TRUE))) %>%
    mutate(across(starts_with(paste0(value_type, "_")), function(x) -log10(x))) %>%
    mutate(across(starts_with(paste0(value_type, "_")), function(x) ifelse(is.infinite(x), 0, x))) %>%
    mutate(n_genes=sapply(strsplit(GeneID, ','), length))

  # representative term making: minimum of the merged value_type column.
  # DS-06 (SPEC-RC-010): the old dplyr filter compared the argument string to
  # itself and was a no-op. The frame here is already a single cluster.
  representative_term <- unname(get_representative_terms(dot_data, value_type))

  # use it in the default title (if none supplied)
  if (is.null(title)) {
    title <- paste0(value_type, ", ", representative_term, " (cluster ", cluster, ")")
  }

  # get all dot data across different dfs
  value_cols <- grep(paste0("^", value_type, "_"), names(dot_data), value = TRUE)

  # the plotly object
  dot <- plotly::plot_ly(
    dot_data,
    x = dot_data[[value_cols[1]]],
    y = ~Term,
    type = 'scatter',
    mode = 'markers',
    marker = list(opacity = 0.5),
    size = ~n_genes,
    name = df_names[1],
    text = ~paste0(
      Term, '<br>',
      '-log10(',value_type,') = ', sprintf("%.2f", dot_data[[value_cols[1]]]), '<br>',
      '# genes = ', n_genes
    ),
    hovertemplate = "%{text}"
  )
  if (length(value_cols) > 1) {
    for (i in 2:length(value_cols)) {
      dot <- dot %>% plotly::add_trace(
        x = dot_data[[value_cols[i]]],
        type = 'scatter',
        marker = list(opacity = 0.5),
        size = ~n_genes,
        name = df_names[i],
        text = ~paste0(
          Term, '<br>',
          '-log10(',value_type,') = ', sprintf("%.2f", dot_data[[value_cols[i]]]), '<br>',
          '# terms = ', n_genes
        ),
        hovertemplate = "%{text}"
      )
    }
  }
  # add layout stuff
  dot <- dot%>% plotly::layout(
    title = title,
    margin = list(autoexpand=TRUE, t=50, b=-50),
    xaxis = list(title = "-log10(Padj)"),
    yaxis = list(title = "Term")
  )
  return(dot)
}

# tdot <- term_dot(cluster_result, cluster=48)
# tdot
