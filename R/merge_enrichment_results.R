# NULL placeholder for roxygen namespace declarations
NULL

#' Merge List of Enrichment Results
#'
#' This function merges multiple enrichment results ('enrichment_results') into a single dataframe by
#' combining unique GeneID elements across each geneset, and averaging Pvalue / Padj
#' values for each term across all enrichment_results.
#'
#' @param enrichment_results A list of geneset dataframes containing columns c('Term', 'GeneID', 'Pvalue', 'Padj')
#'
#' @return A single merged geneset dataframe with all original columns
#'         suffixed with the index of the geneset, with new columns 'GeneID', 'Pvalue',
#'         'Padj' containing the merged values.
#'
#' @export
merge_enrichment_results <- function(enrichment_results) {
  # Note: Keep track of what index each geneset has in the list

  SEP <- "_" # separator for column suffixes (geneID + '_' + index)

  # Preprocessing: Suffix all non 'Term' columns by their index in the list
  # Allows base::merge by 'Term'
  for (i in seq_along(enrichment_results)) {
    rownames(enrichment_results[[i]]) <- NULL # prevents rownames from causing errors

    # Get renamed columns (suffixed by index)
    colnames(enrichment_results[[i]]) <- format_colnames(colnames(enrichment_results[[i]]))
    all_colnames <- colnames(enrichment_results[[i]])
    nonterm_cols <- all_colnames[all_colnames != 'Term']
    all_colnames[all_colnames != 'Term'] <- paste(nonterm_cols, i, sep=SEP)

    colnames(enrichment_results[[i]]) <- all_colnames
  }

  # Initialize merged_gs with first geneset
  merged_gs <- enrichment_results[[1]]
  if (length(enrichment_results) == 1) {
    return(merged_gs) # Return the single geneset if only one is provided
  }

  # Else, merge the rest of the enrichment_results
  for (i in 2:length(enrichment_results)) {
    merged_gs <- base::merge(merged_gs, enrichment_results[[i]], by='Term', all=TRUE)
  }

  # For each row in merged_gs, combine unique GeneID elements
  geneid_cols <- paste("GeneID", seq_along(enrichment_results), sep=SEP)
  merged_gs$GeneID <- apply(merged_gs[, geneid_cols], 1, function(x) {
    paste(unique(na.omit(x)), collapse = ',')
  })

  # Average the value columns across all enrichment_results
  # Avg Pvalue
  pvalue_cols <- paste("Pvalue", seq_along(enrichment_results), sep=SEP)
  available_pvalue_cols <- pvalue_cols[pvalue_cols %in% colnames(merged_gs)]
  if (length(available_pvalue_cols) > 0) {
    merged_gs$Pvalue <- rowMeans(merged_gs[, available_pvalue_cols], na.rm = TRUE)
  }
  
  # Avg Padj
  padj_cols <- paste("Padj", seq_along(enrichment_results), sep=SEP)
  available_padj_cols <- padj_cols[padj_cols %in% colnames(merged_gs)]
  if (length(available_padj_cols) > 0) {
    merged_gs$Padj <- rowMeans(merged_gs[, available_padj_cols], na.rm = TRUE)
  }
  # Return the merged geneset df
  return(merged_gs)

}

#' Format Column Names for Merging
#'
#' This function maps a vector of column names to standardized names
#' for "GeneID", "Pvalue", and "Padj" based on known variations.
#'
#' @param colnames A character vector of column names to be standardized.
#'
#' @return A character vector of standardized column names.
format_colnames <- function(colnames) {
  # dictionary of common alternative column names
  mappings <- list(
    Term = c("term", "pathway", "keyword", "domain", "description", "title"),
    GeneID = c("geneid", "gene", "gene_symbols", "gene_id"),
    Pvalue = c("pvalue", "pval", "p-value", "value"),
    Padj   = c("padj", "p-adj", "pvalue_adjusted", "pval_adj", "adj_pvalue")
  )
  # look through dictionary to find matching name
  get_good_name <- function(colname) {
    for (good_name in names(mappings)) {
      if (tolower(colname) %in% mappings[[good_name]]) {
        return(good_name)
      }
    }
    return(colname) # return original name if none match
  }
  # apply function to all colnames
  return(sapply(colnames, get_good_name))
}
