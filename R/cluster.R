#' @importFrom magrittr %>%
#' @importFrom dplyr across filter group_by mutate n row_number select summarise ungroup
#' @importFrom tidyr separate_rows
NULL

#' Cluster Terms from Enrichment Results
#'
#' This function performs clustering on enrichment results by integrating
#' gene similarity scores and various clustering strategies.
#'
#' @param enrichment_results A list of dataframes, each containing enrichment results.
#'        Each dataframe should include at least the columns 'Term', 'GeneID', and 'Padj'.
#' @param df_names Optional, a character vector of names for the enrichment result dataframes. Must
#'        match the length of `enrichment_results`. Default is `NULL`.
#' @param min_terms Minimum number of terms each final cluster must include
#' @param min_value Upper bound on the significance value: a term is kept when
#'        its `filter_on` column is strictly LESS than `min_value`. Despite the
#'        name this is a maximum, not a minimum; the name is retained from 1.0.2
#'        so that existing calls keep working. Default is `0.1`.
#' @param distance_metric A string specifying the distance metric to use.
#'        Supported options are "kappa", "jaccard", and "dice".
#' @param distance_cutoff A numeric value for the distance cutoff (0 < cutoff <= 1).
#'        The comparison is STRICT: two terms are linked when their similarity
#'        is greater than `distance_cutoff`, not when it equals it.
#' @param linkage_method A string specifying the linkage method to use
#'        (e.g., "average"). Supported options are "single", "complete",
#'        "average", and "ward".
#' @param linkage_cutoff A numeric value between 0 and 1 for the membership cutoff.
#'        The comparison is STRICT: clusters merge when their linkage score is
#'        greater than `linkage_cutoff`, not when it equals it.
#' @param filter_on Name of the column terms are selected on, compared against
#'        `min_value`. Default is `"Padj"`; `"Pvalue"` reproduces 1.0.2
#'        behaviour. When the named column is absent but 'Pvalue' is present,
#'        `cluster()` warns (condition class
#'        `richCluster_filter_on_fallback`) and falls back to 'Pvalue' rather
#'        than erroring.
#' @param gene_delim A single string separating gene identifiers within the
#'        'GeneID' column. Default is `","`, which reproduces 1.0.2 behaviour
#'        exactly. Set it when your gene lists use another separator (for
#'        example `";"`); it is applied both to the per-term gene split and to
#'        the distinct-gene universe count, which must agree.
#'
#' @param verbose Logical; print the C++ core's progress narration to the
#'        console. Default `FALSE` --- the core runs silent. Before 2.0.0 this
#'        narration was unconditional and could not be switched off, and it
#'        scales with the data (one line per merge iteration), so it is now
#'        opt-in.
#'
#' @return A named list of eight elements, in this order:
#' \describe{
#'   \item{distance_matrix}{A numeric \code{n x n} matrix of pairwise SIMILARITY
#'     scores between the \code{n} terms that survived the \code{min_value}
#'     filter. Values lie on \[0, 1\] and LARGER means MORE similar; the diagonal
#'     is 1. For \code{"kappa"}, a negative value is set to 0, so the matrix
#'     holds a non-negative kappa similarity rather than an unmodified Cohen's
#'     kappa. Row and column names are the terms, in \code{merged_df} row order.
#'     Despite the element name it holds similarities, not distances; the name is
#'     retained so that existing code keeps working.}
#'   \item{all_clusters}{A data frame of every merged cluster, BEFORE the
#'     \code{min_terms} filter, with three columns: \code{Cluster} (integer
#'     cluster ID), \code{TermNames} (the cluster's term names, comma-separated)
#'     and \code{TermIndices} (the same terms as ZERO-based row indices into
#'     \code{merged_df}, comma-separated).}
#'   \item{df_list}{The \code{enrichment_results} argument, unmodified.}
#'   \item{merged_df}{The merged and filtered data frame the clustering ran on:
#'     one row per surviving term, carrying \code{Term}, the per-dataset columns
#'     suffixed \code{_1}, \code{_2}, ... , and the pooled \code{GeneID},
#'     \code{Pvalue}, \code{Padj} and \code{DatasetCount} columns. Its row
#'     order is the index basis of \code{TermIndices}.}
#'   \item{cluster_options}{A named list of the parameters clustering ran with:
#'     \code{min_terms}, \code{min_value}, \code{distance_metric},
#'     \code{distance_cutoff}, \code{linkage_method}, \code{linkage_cutoff} and
#'     \code{filter_on}.}
#'   \item{df_names}{A character vector naming the input data frames. Always
#'     present: when \code{df_names} is \code{NULL} or does not match the length
#'     of \code{enrichment_results} it is replaced with \code{"1"},
#'     \code{"2"}, ... .}
#'   \item{final_clusters}{\code{all_clusters} with clusters of fewer than
#'     \code{min_terms} terms dropped and \code{Cluster} renumbered 1..k. TWO
#'     columns only --- \code{Cluster} and \code{TermIndices}; \code{TermNames}
#'     is not carried over.}
#'   \item{cluster_df}{The main output, and the frame the plotting functions
#'     (\code{cluster_hmap()}, \code{cluster_dot()}, \code{cluster_bar()},
#'     \code{cluster_network()}) consume: one row per (cluster, term) pair, being
#'     a \code{Cluster} column followed by that term's row from
#'     \code{merged_df}.}
#' }
#'
#' @section Determinism:
#' Clustering results depend only on the terms and gene sets supplied, not on
#' the order in which rows are given: \code{cluster()} and
#' \code{runRichCluster()} canonicalise row order internally, so re-sorting the
#' input cannot change cluster membership. Note that the merge stage is a greedy
#' agglomeration --- clusters are built by repeatedly merging the best-scoring
#' available pair --- so results reflect that greedy strategy rather than a
#' global optimum, and small changes to \code{linkage_cutoff} can change
#' membership substantially.
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
#' res <- cluster(list(d1, d2), min_terms = 2, distance_metric = "kappa")
#' head(res$cluster_df[, c("Cluster", "Term")])
#' }
#' @export
cluster <- function(enrichment_results, df_names=NULL, min_terms=5, min_value=0.1,
                    distance_metric="kappa", distance_cutoff=0.5,
                    linkage_method="average", linkage_cutoff=0.5,
                    filter_on="Padj", gene_delim=",", verbose=FALSE) {

  if (is.null(df_names) || length(enrichment_results) != length(df_names)) {
    df_names <- as.character(seq_along(enrichment_results))
  }

  validate_inputs(enrichment_results, df_names, distance_metric, distance_cutoff,
                  linkage_method, linkage_cutoff, filter_on)

  # accept a list of dataframes as input
  # call merge_enrichment_results
  # gene_delim must reach the MERGE as well as the C++ split.  Threading it
  # into only one of the two joins gene lists with "," and then splits them
  # with the caller's delimiter, fabricating identifiers present in no input.
  merged_df <- merge_enrichment_results(enrichment_results, gene_delim = gene_delim)

  # SPEC-RC-007 / T1-17.  1.0.2 filtered on the RAW Pvalue, hardcoded, while
  # cluster()'s own documentation named Padj as a required column.  The default
  # is now Padj; "Pvalue" reproduces 1.0.2 exactly.
  #
  # REQ-017-3 / REQ-017-4: a missing requested column WARNS and falls back to
  # Pvalue -- it must never error.  Frames carrying Pvalue but no Padj ran fine
  # at 1.0.2, and erroring on them would break working input for no gain.  The
  # 2.0.0 major bump lifts the semver obligation to keep this, but not the
  # reason for it: the fallback is deliberate, not inherited.
  filter_col <- filter_on
  if (!filter_col %in% names(merged_df)) {
    if ("Pvalue" %in% names(merged_df)) {
      warning(structure(
        class = c("richCluster_filter_on_fallback", "warning", "condition"),
        list(message = sprintf(
               "filter_on = '%s' not found in merged results; falling back to 'Pvalue'.",
               filter_on),
             call = NULL)))
      filter_col <- "Pvalue"
    } else {
      stop(sprintf(
        "Neither filter_on = '%s' nor 'Pvalue' is present in the merged results.",
        filter_on))
    }
  }
  # Base R rather than dplyr, because the column is chosen at runtime and a bare
  # symbol cannot be.  The explicit !is.na() and the rowname reset are what make
  # this byte-equivalent to the 1.0.2 `filter(Pvalue < min_value)` it replaces:
  # dplyr::filter() treats NA as FALSE and renumbers rows, and `[` does neither.
  keep <- !is.na(merged_df[[filter_col]]) & merged_df[[filter_col]] < min_value
  merged_df <- merged_df[keep, , drop = FALSE]
  rownames(merged_df) <- NULL

  
  term_vec <- merged_df$Term
  geneID_vec <- merged_df$GeneID

  # throw error if cluster options are invalid

  cluster_result <- richCluster::runRichCluster(
    term_vec, geneID_vec,
    distance_metric, distance_cutoff,
    linkage_method, linkage_cutoff,
    gene_delim, verbose
  )

  # add the original stuff to the cluster_result
  # (helps visualizations later)
  cluster_options <- list(
    min_terms = min_terms,
    min_value = min_value,
    distance_metric = distance_metric,
    distance_cutoff = distance_cutoff,
    linkage_method = linkage_method,
    linkage_cutoff = linkage_cutoff,
    filter_on = filter_on
  )

  cluster_result$df_list <- enrichment_results
  cluster_result$merged_df <- merged_df
  cluster_result$cluster_options <- cluster_options
  cluster_result$df_names <- df_names

  cluster_result$final_clusters <- filter_clusters(cluster_result$all_clusters, min_terms)
  cluster_result$cluster_df <- make_full_clusterdf(cluster_result$final_clusters, merged_df)

  return(cluster_result)
}


validate_inputs <- function(enrichment_results, df_names=NA_character_,
                            distance_metric="kappa", distance_cutoff=0.5,
                            linkage_method="average", linkage_cutoff=0.5,
                            filter_on="Padj") {
  if (!is.list(enrichment_results)) {
    stop("enrichment_results must be a list of dataframes.")
  }
  if (any(!sapply(enrichment_results, is.data.frame))) {
    stop("Each element of enrichment_results must be a dataframe.")
  }
  if (distance_cutoff <= 0 || distance_cutoff > 1) {
    stop("distance_cutoff must be between 0 and 1.")
  }
  if (linkage_cutoff <= 0 || linkage_cutoff > 1) {
    stop("linkage_cutoff must be between 0 and 1.")
  }
  if (!distance_metric %in% c("kappa", "jaccard", "dice")) {
    stop("Unsupported distance metric. Only 'kappa', 'jaccard', and 'dice' are supported.")
  }
  if (!linkage_method %in% c("single", "complete", "average", "ward")) {
    stop("Unsupported linkage_method. Only 'single', 'complete', 'average', and 'ward' are supported.")
  }
  # REQ-017-8: closed set, matching the whitelist style used just above.
  if (!filter_on %in% c("Padj", "Pvalue")) {
    stop("Unsupported filter_on. Only 'Padj' and 'Pvalue' are supported.")
  }

}

#' Filter Clusters by Number of Terms
#'
#' Filters the full list of clusters by keeping only those with greater
#' than or equal to min_terms # of terms.
#'
#' @param all_clusters A dataframe containing the merged seeds with column named `ClusterIndices`.
#' @param min_terms An integer specifying the minimum number of terms required in a cluster.
#'
#' @return The filtered data frame with clusters filtered to include only those with at least `min_terms` terms.
#'
#' @examples
#' \donttest{
#' cluster_result <- readRDS(system.file("extdata", "cluster_result.rds",
#'                                       package = "richCluster"))
#' # Keep only clusters carrying at least 10 terms.
#' kept <- filter_clusters(cluster_result$all_clusters, min_terms = 10)
#' nrow(cluster_result$all_clusters)
#' nrow(kept)
#' }
#' @export
filter_clusters <- function(all_clusters, min_terms)
{
  filtered_clusters <- all_clusters %>%
    mutate(row_id = row_number()) %>%  # Add a row identifier
    separate_rows(TermIndices, sep = ", ") %>%  # Separate into individual rows
    group_by(row_id, Cluster) %>%  # Group by the original rows
    dplyr::filter(n() >= min_terms) %>%  # Filter groups with at least X terms
    summarise(TermIndices = paste(TermIndices, collapse = ", ")) %>%  # Collapse back to single strings
    ungroup() %>%  # Ungroup to finalize the data frame
    mutate(Cluster = row_number()) %>%  # DS-05 (SPEC-RC-010): canonical ids 1..k
    select(-row_id)  # Remove the temporary row identifier

  return(filtered_clusters)
}


make_full_clusterdf <- function(final_clusters, merged_df) {
  # Initialize an empty data frame to store the results
  full_clusterdf <- data.frame()

  # Loop over each row in final_clusters
  for(i in seq_len(nrow(final_clusters))) {
    row <- final_clusters[i, ]
    TermIndices <- unlist(strsplit(row$TermIndices, ", "))

    # Loop over each term index in TermIndices
    for (termIndex in TermIndices) {
      R_termIndex <- as.integer(termIndex) + 1  # Convert termIndex to integer and adjust for 1-based indexing
      term_row <- merged_df[R_termIndex, ]  # Get the row corresponding to the termIndex

      # Create a new row with the cluster number and term row.
      # DS-05 (SPEC-RC-010): label with the id the final_clusters row carries,
      # not the loop counter, so the agreement is structural not coincidental.
      new_row <- c(Cluster = row$Cluster, term_row)

      # Append the new row to the data frame
      full_clusterdf <- rbind(full_clusterdf, new_row)
    }
  }

  return(full_clusterdf)
}


# Representative term per cluster: the Term carrying the minimum value of the
# merged significance column named by value_type ("Padj"/"Pvalue"). One term
# per cluster, ordered by ascending cluster id (matches group_by(Cluster)
# summarise order at every call site). DS-05/DS-06: SPEC-RC-010.
# C9 (converge ledger, 2026-09-04): a bit-exact tie on value_type is broken by
# the OTHER significance column, then by Term (bytewise), never by row
# position -- row position followed the C++ standard library's hash order and
# flipped the choice between macOS and Linux on the shipped example.
get_representative_terms <- function(cluster_df, value_type) {
  if (!value_type %in% names(cluster_df)) {
    stop(sprintf("value_type '%s' is not a column of cluster_df", value_type))
  }
  vals  <- cluster_df[[value_type]]
  other <- setdiff(c("Padj", "Pvalue"), value_type)
  tie   <- if (length(other) == 1L && other %in% names(cluster_df)) {
    cluster_df[[other]]
  } else {
    rep(NA_real_, nrow(cluster_df))
  }
  terms <- cluster_df$Term
  vapply(split(seq_len(nrow(cluster_df)), cluster_df$Cluster), function(idx) {
    v <- vals[idx]
    if (all(is.na(v))) return(terms[idx[1]])
    o <- order(v, tie[idx], terms[idx], method = "radix", na.last = TRUE)
    terms[idx[o[1]]]
  }, character(1))
}



#' Run clustering in C++ backend
#'
#' @param terms Character vector of term names
#' @param geneIDs Character vector of geneIDs
#' @param distanceMetric e.g. "kappa"
#' @param distanceCutoff numeric between 0 and 1
#' @param linkageMethod e.g. "average"
#' @param linkageCutoff numeric between 0 and 1
#' @param geneDelim single string separating gene identifiers within `geneIDs`.
#'        Default `","`.
#'
#' @param verbose Logical; print the C++ core's progress narration to the
#'        console. Default `FALSE` --- the core runs silent. Before 2.0.0 this
#'        narration was unconditional and could not be switched off, and it
#'        scales with the data (one line per merge iteration), so it is now
#'        opt-in.
#'
#' @return A list of two elements:
#' \describe{
#'   \item{distance_matrix}{A numeric \code{n x n} matrix of pairwise SIMILARITY
#'     scores between the \code{n} input terms. Values lie on \[0, 1\] and LARGER
#'     means MORE similar; the diagonal is 1, the self-similarity every supported
#'     metric agrees on. For \code{"kappa"}, a negative value is set to 0, so the
#'     matrix holds a non-negative kappa similarity rather than an unmodified
#'     Cohen's kappa. Row and column names are \code{terms}, in the caller's
#'     input order. Despite the element name it holds similarities, not
#'     distances; the name is retained so that existing code keeps working.
#'     Linkage converts internally with \code{d = 1 - similarity}.}
#'   \item{all_clusters}{A data frame of every merged cluster, with three
#'     columns: \code{Cluster} (integer cluster ID), \code{TermNames} (the
#'     cluster's term names, comma-separated) and \code{TermIndices} (the same
#'     terms as ZERO-based indices into \code{terms}, comma-separated).}
#' }
#'
#' @section Determinism:
#' Clustering results depend only on the terms and gene sets supplied, not on
#' the order in which rows are given: \code{cluster()} and
#' \code{runRichCluster()} canonicalise row order internally, so re-sorting the
#' input cannot change cluster membership. Note that the merge stage is a greedy
#' agglomeration --- clusters are built by repeatedly merging the best-scoring
#' available pair --- so results reflect that greedy strategy rather than a
#' global optimum, and small changes to \code{linkage_cutoff} can change
#' membership substantially.
#'
#'
#' @details
#' A missing (\code{NA}) or empty gene list is treated as the EMPTY SET: the
#' term is similar to nothing, so it joins no cluster and \code{min_terms}
#' drops it. No row is removed from the caller's data and no count changes.
#' A classed warning (\code{richCluster_empty_gene_list}) names how many terms
#' were affected, so the exclusion is never silent.
#' @examples
#' \donttest{
#' res <- runRichCluster(
#'   terms = c("T1", "T2", "T3"),
#'   geneIDs = c("a,b,c", "b,c,d", "x,y,z"),
#'   distanceMetric = "kappa", distanceCutoff = 0.5,
#'   linkageMethod = "average", linkageCutoff = 0.5)
#' res$distance_matrix
#' }
#' @export
runRichCluster <- function(terms, geneIDs, distanceMetric, distanceCutoff, linkageMethod, linkageCutoff, geneDelim = ",", verbose = FALSE) {
  # OD-7.  NA reaches Rcpp::as and becomes the literal string "NA", so two terms
  # with missing gene lists share that single "gene" and score similarity 1.0 --
  # they cluster together on the strength of both being unknown.  The C++ already
  # does the right thing for an empty string (similarity 0); it was simply never
  # reached for NA.  Normalising at this boundary is the whole fix.
  #
  # The warning fires only where an NA was actually normalised.  cluster() merges
  # first, and the merge has already collapsed all-NA lists to "", so this path
  # sees no NA and does not warn a second time about the same terms.
  geneIDs <- as.character(geneIDs)
  na_hits <- is.na(geneIDs)
  if (any(na_hits)) {
    geneIDs[na_hits] <- ""
    rc_warn_empty_gene_lists(sum(na_hits), "runRichCluster()")
  }
  .Call(`_richCluster_runRichCluster`, terms, geneIDs, distanceMetric, distanceCutoff, linkageMethod, linkageCutoff, geneDelim, verbose)
}
