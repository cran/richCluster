#' @importFrom magrittr %>%
#' @importFrom dplyr across filter group_by mutate pull slice starts_with summarise ungroup
NULL

#' Cluster-level Bar Plot of Enrichment Significance
#'
#' Generates a horizontal bar plot showing average enrichment significance
#' for each cluster, across one or more enrichment datasets.
#'
#' @param cluster_result A result list returned by \code{\link{cluster}}.
#' @param clusters Optional numeric vector of cluster IDs to include. Defaults to all clusters.
#' @param value_type The column name to use for enrichment significance ("Padj" or "Pvalue").
#' @param title Optional plot title. If NULL, a default will be generated.
#'
#' @return A \code{plotly} object representing the bar plot.
#'
#' @examples
#' \donttest{
#' # Load example data
#' cluster_result <- readRDS(system.file("extdata", "cluster_result.rds",
#'                                       package = "richCluster"))
#' cbar <- cluster_bar(cluster_result)
#' cbar
#' }
#' @export
cluster_bar <- function(cluster_result, clusters=NULL, value_type="Padj", title=NULL) {

  cluster_df <- cluster_result$cluster_df
  df_names <- cluster_result$df_names

  # bar data
  bar_data <- cluster_df %>%
    group_by(Cluster) %>%
    summarise(across(starts_with(paste0(value_type, "_")), function(x) mean(x, na.rm=TRUE))) %>%
    mutate(across(starts_with(paste0(value_type, "_")), function(x) -log10(x))) %>%
    mutate(across(starts_with(paste0(value_type, "_")), function(x) ifelse(is.infinite(x), 0, x)))

  # representative term making
  # use the Pvalue/Padj average column
  representative_terms <- cluster_df %>%
    group_by(Cluster) %>%
    filter(value_type==min(value_type, na.rm=TRUE)) %>%
    slice(1) %>%
    ungroup %>%
    pull(Term)

  bar_data$Cluster <- representative_terms

  # generate default title if none supplied
  if (is.null(title)) {
    title <- paste0(value_type, ", all clusters")
  }

  # get all bar data across different dfs
  value_cols <- grep(paste0("^", value_type, "_"), names(bar_data), value = TRUE)

  # the plotly object
  bar <- plotly::plot_ly(
    bar_data,
    x = bar_data[[value_cols[1]]],
    y = ~Cluster,
    type = 'bar',
    name = df_names[1]
  )
  if (length(value_cols) > 1) {
    for (i in 2:length(value_cols)) {
      bar <- bar %>% plotly::add_trace(
        x = bar_data[[value_cols[i]]],
        name = df_names[i]
      )
    }
  }
  # add layout stuff
  bar <- bar%>% plotly::layout(
    title = title,
    margin = list(autoexpand=TRUE, t=50, b=-50),
    xaxis = list(title = "-log10(Padj)"),
    yaxis = list(title = "Term")
  )
  return(bar)
}

# cbar <- cluster_bar(cluster_result)
# cbar

#' Term-level Bar Plot for a Specific Cluster
#'
#' Creates a horizontal bar plot showing enrichment values for individual terms
#' in a selected cluster.
#'
#' @param cluster_result A result list returned by \code{\link{cluster}}.
#' @param cluster Cluster ID (numeric) or term name (character) to visualize.
#' @param value_type The column name to use for enrichment significance ("Padj" or "Pvalue").
#' @param title Optional plot title. If NULL, a default will be generated.
#'
#' @return A \code{plotly} object representing the bar plot.
#'
#' @examples
#' \donttest{
#' # Load example data
#' cluster_result <- readRDS(system.file("extdata", "cluster_result.rds",
#'                                       package = "richCluster"))
#' tbar <- term_bar(cluster_result, cluster = 1)
#' tbar
#' }
#' @export
term_bar <- function(cluster_result, cluster=1, value_type="Padj", title=NULL) {

  cluster_df <- cluster_result$cluster_df
  df_names <- cluster_result$df_names

  # get cluster # from the term
  if (is.character(cluster)) {
    cluster <- cluster_df[cluster_df$Term==cluster, ]$Cluster
  } else if (!is.numeric(cluster)) {
    stop("cluster must be numeric (a cluster number) or character (a term name)")
  }

  # bar data
  bar_data <- cluster_df %>%
    group_by(Cluster) %>%
    filter(Cluster==cluster) %>%
    # summarise(across(starts_with(paste0(value_type, "_")), function(x) mean(x, na.rm=TRUE))) %>%
    mutate(across(starts_with(paste0(value_type, "_")), function(x) -log10(x))) %>%
    mutate(across(starts_with(paste0(value_type, "_")), function(x) ifelse(is.infinite(x), 0, x)))

  # representative term making
  # use the Pvalue/Padj average column
  representative_term <- bar_data %>%
    group_by(Cluster) %>%
    filter(value_type==min(value_type, na.rm=TRUE)) %>%
    slice(1) %>%
    ungroup %>%
    pull(Term)

  # use it in the default title (if none supplied)
  if (is.null(title)) {
    title <- paste0(value_type, ", ", representative_term, " (cluster ", cluster, ")")
  }

  # get all bar data across different dfs
  value_cols <- grep(paste0("^", value_type, "_"), names(bar_data), value = TRUE)

  # the plotly object
  bar <- plotly::plot_ly(
    bar_data,
    x = bar_data[[value_cols[1]]],
    y = ~Term,
    type = 'bar',
    name = df_names[1]
  )
  if (length(value_cols) > 1) {
    for (i in 2:length(value_cols)) {
      bar <- bar %>% plotly::add_trace(
        x = bar_data[[value_cols[i]]],
        name = df_names[i]
      )
    }
  }
  # add layout stuff
  bar <- bar%>% plotly::layout(
    title = title,
    margin = list(autoexpand=TRUE, t=50, b=-50),
    xaxis = list(title = "-log10(Padj)"),
    yaxis = list(title = "Term")
  )
  return(bar)
}

# tbar <- term_bar(cluster_result, cluster=48)
# tbar
