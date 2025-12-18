# NULL placeholder for roxygen namespace declarations
NULL

#' Create a Correlation Heatmap for a Specific Cluster
#'
#' This function generates a correlation heatmap for a specific cluster based on the provided distance matrix.
#'
#' @param final_clusters A dataframe containing the final cluster data.
#' @param distance_matrix A matrix representing the distances between terms.
#' @param cluster_number An integer specifying the cluster number to visualize.
#' @param merged_df A dataframe with all terms used to map term indices to names.
#' @return An interactive heatmaply heatmap.
#' @export
cluster_correlation_hmap <- function(final_clusters, distance_matrix, cluster_number, merged_df) {
  #TODO: Be able to see what genes are involved
  #TODO: Update tooltip

  # Extract and process ClusterIndices
  term_indices <- as.numeric(unlist(strsplit(final_clusters$TermIndices[cluster_number], ", ")))

  # Create an empty matrix for the cluster
  cluster_matrix <- matrix(0, nrow = length(term_indices), ncol = length(term_indices))

  # Populate the cluster matrix with distances from the distance matrix
  for (i in seq_along(term_indices)) {
    term_i <- term_indices[i] + 1  # Adjust for 1-based indexing
    for (j in seq_along(term_indices)) {
      term_j <- term_indices[j] + 1  # Adjust for 1-based indexing
      kappa <- distance_matrix[term_i, term_j]
      if (kappa == -99) {
        kappa <- 1
      }
      cluster_matrix[i, j] <- kappa
    }
  }

  # Get names of all terms from term_indices
  term_names <- merged_df$Term[term_indices + 1]  # Adjust for 1-based indexing

  rownames(cluster_matrix) <- term_names
  colnames(cluster_matrix) <- term_names

  # Create the heatmaply plot
  c <- heatmaply::heatmaply(
    cluster_matrix,
    main = paste0("Correlation Matrix for Cluster ", cluster_number),
    xlab = "Terms",
    ylab = "Terms",
    colors = viridis::viridis(256),
    grid_color = "grey",
    dendrogram = "none",
    margins = c(50, 50, 50, 50),
    limits = c(min(cluster_matrix, na.rm = TRUE), max(cluster_matrix, na.rm = TRUE))
  )
  return(c)
}

# Example usage
# c <- cluster_correlation_hmap(final_clusters, distance_matrix, 7)
# c


#' Create a Network Graph for a Specific Cluster
#'
#' This function generates a network graph for a specific cluster based on the provided distance matrix.
#' The opacity and length of the edges correspond to the given distance_metric (eg, kappa) score similarity 
#' between terms, which is based on shared gene content.
#'
#' @param final_clusters A dataframe containing the final cluster data.
#' @param distance_matrix A matrix representing the distances between terms.
#' @param cluster_number An integer specifying the cluster number to visualize.
#' @param merged_df A dataframe with all terms used to map term indices to names.
#' @return An interactive networkD3 network graph.
#' @export
cluster_network <- function(final_clusters, distance_matrix, cluster_number, merged_df) {
  # view a network graph of all terms in a single cluster
  # opacity + length of edge corresponds to kappa score
  #TODO: double check does higher=shorter+darker?

  # Extract and process ClusterIndices
  term_indices <- as.numeric(unlist(strsplit(final_clusters$TermIndices[cluster_number], ", ")))

  # Create an empty matrix for the cluster
  cluster_matrix <- matrix(0, nrow = length(term_indices), ncol = length(term_indices))

  # Populate the cluster matrix with distances from the distance matrix
  for (i in seq_along(term_indices)) {
    term_i <- term_indices[i] + 1  # Adjust for 1-based indexing
    for (j in seq_along(term_indices)) {
      term_j <- term_indices[j] + 1  # Adjust for 1-based indexing
      kappa <- distance_matrix[term_i, term_j]
      if (kappa == -99) {
        kappa <- 1
      }
      cluster_matrix[i, j] <- kappa
    }
  }

  # Get names of all terms from term_indices
  term_names <- merged_df$Term[term_indices + 1]  # Adjust for 1-based indexing
  rownames(cluster_matrix) <- term_names
  colnames(cluster_matrix) <- term_names

  # Create an igraph object from the cluster matrix
  g <- igraph::graph_from_adjacency_matrix(cluster_matrix, mode = "undirected", weighted = TRUE)

  # Set the vertex attributes
  g_d3 <- networkD3::igraph_to_networkD3(g, term_names)

  # Create character vector of link colors
  link_colors <- sapply(g_d3$links$value, function(value) {
    opacity <- value  # Assuming value is between 0 and 1
    paste0("rgba(0, 0, 0, ", opacity, ")")
  })

  # Create the networkD3 plot
  d3net <-  networkD3::forceNetwork(
    Links = g_d3$links,
    Nodes = g_d3$nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "name",
    Group = "group",
    opacity = 0.9,
    linkWidth = 1.5,
    linkDistance = networkD3::JS("function(d) { return d.value * 50; }"),  # Adjust link distance based on value
    charge = -200,
    fontFamily = "arial",
    zoom = TRUE,
    linkColour = link_colors
    # tooltip = JS("function(d) { return 'Value: ' + d.value; }")  # Show value on hover
  )
  return(d3net)
}

# Example usage
# d <- cluster_network(final_clusters, distance_matrix, 20)
# d
