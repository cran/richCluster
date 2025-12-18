#' @importFrom magrittr %>%
#' @importFrom dplyr filter group_by mutate pull select slice summarise ungroup everything
NULL

# functions to export the cluster_result as csv

#' Export Cluster Result as Dataframe
#'
#' Returns a comprehensive dataframe containing all the different terms in all clusters.
#'
#' @param cluster_result The cluster_result object from cluster()
#'
#' @returns A data.frame view of the clustering
#' @export
export_df <- function(cluster_result) {
  cluster_df <- cluster_result$cluster_df
  df_names <- cluster_result$df_names

  # create representative term column (ClusterName)
  representative_terms <- cluster_df %>%
    group_by(Cluster) %>%
    filter(Pvalue==min(Pvalue, na.rm=TRUE)) %>%
    slice(1) %>%
    ungroup %>%
    pull(Term)

  cluster_df <- cluster_df %>%
    mutate(ClusterName=representative_terms[Cluster]) %>%
    select(Cluster, ClusterName, everything())

  # replace _* with _dfname[*]
  df_colnames <- colnames(cluster_df)
  df_colnames <- sapply(
    df_colnames,
    function(colname) {
      if (!grepl("_\\d+$", colname)) {
        return(colname) # return unchanged if not _*
      }
      index <- sub(".*_(\\d+)$", "\\1", colname)
      # if NA index
      if (is.na(index)) {
        return(colname)
      }
      index <- as.numeric(index)
      if (index > length(df_names)) {
        return(colname)
      } else { # replace!
        sub("_\\d+$", paste0("_", df_names[index]), colname)
      }
    }
  )
  colnames(cluster_df) <- df_colnames
  return(cluster_df)

}
