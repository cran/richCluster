# NULL placeholder for roxygen namespace declarations
NULL

#' Create a Network Graph for the Entire Distance Matrix
#'
#' This function generates a network graph for the entire distance matrix.
#'
#' @param cluster_result Cluster result named list from richCluster::cluster()
#' @return An interactive networkD3 network graph.
#' @export
full_network <- function(cluster_result) {
  distance_matrix = cluster_result$distance_matrix
  # Create an igraph object from the distance matrix
  g <- igraph::graph_from_adjacency_matrix(distance_matrix, mode="undirected", weighted=TRUE)
  term_names <- rownames(distance_matrix)
  # Create the networkD3 plot
  g_d3 <- networkD3::igraph_to_networkD3(g, term_names)
  d3net <- networkD3::forceNetwork(
    Links = g_d3$links,
    Nodes = g_d3$nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "name",
    Group = "group",
    opacity = 0.9,
    linkWidth = 1.5,
    linkDistance =  networkD3::JS("function(d) { return d.value * 50; }"),  # Adjust link distance based on value
    charge = -200,
    fontFamily = "arial",
    zoom = TRUE
  )
  return(d3net)
}

# Example usage
# get only first 100 x 100 matrix of distance matrix
# small_dm <- distance_matrix[1:100, 1:100]
# d_full <- full_network(small_dm)
# d_full
