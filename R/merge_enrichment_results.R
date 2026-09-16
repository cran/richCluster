# NULL placeholder for roxygen namespace declarations
NULL

#' Merge List of Enrichment Results
#'
#' This function merges multiple enrichment results ('enrichment_results') into a single dataframe by
#' combining unique GeneID elements across each geneset, and averaging Pvalue / Padj
#' values for each term across all enrichment_results.
#'
#' Each gene identifier appears at most once per term, ordered by first
#' appearance across the contributing genesets.
#'
#' @param enrichment_results A list of geneset dataframes containing columns c('Term', 'GeneID', 'Pvalue', 'Padj').
#'        A list of length 1 is supported and returns the same merged columns as a
#'        longer list, row-ordered by 'Term'.
#'
#' @param gene_delim A single non-empty string separating gene identifiers
#'        within the 'GeneID' column. Must match the delimiter the caller's
#'        data actually uses: the merge splits on it and re-joins on it, so a
#'        mismatch fabricates identifiers that appear in no input. Defaults to
#'        ",", which reproduces the pre-parameter behaviour exactly.
#'
#' @return A single merged geneset dataframe with all original columns
#'         suffixed with the index of the geneset, with new columns 'GeneID', 'Pvalue',
#'         'Padj' containing the merged values, and a trailing integer column
#'         'DatasetCount' giving the number of input datasets that reported each term.
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
#' d1 <- utils::read.delim(system.file("extdata", "HF36wk_vs_HF12wk.txt",
#'                                      package = "richCluster"))
#' d2 <- utils::read.delim(system.file("extdata", "HF36wk_vs_WT12wk.txt",
#'                                      package = "richCluster"))
#' merged <- merge_enrichment_results(list(d1, d2))
#' head(merged[, c("Term", "Pvalue", "Padj", "DatasetCount")])
#' }
#' @export
merge_enrichment_results <- function(enrichment_results, gene_delim = ",") {

  # The merge splits and re-joins gene lists.  Under an empty delimiter
  # strsplit() returns one element PER CHARACTER rather than erroring, so an
  # unusable value corrupts silently instead of failing -- the same shape the
  # C++ entry guard blocks at src/RichCluster.cpp.  Refuse it here too, before
  # any work happens.
  if (!is.character(gene_delim) || length(gene_delim) != 1L ||
      is.na(gene_delim) || !nzchar(gene_delim)) {
    stop("gene_delim must be a single non-empty string", call. = FALSE)
  }
  # Note: Keep track of what index each geneset has in the list

  SEP <- "_" # separator for column suffixes (geneID + '_' + index)

  # Preprocessing: Suffix all non 'Term' columns by their index in the list
  # Allows base::merge by 'Term'
  for (i in seq_along(enrichment_results)) {
    rownames(enrichment_results[[i]]) <- NULL # prevents rownames from causing errors

    # Presence marker: added before the rename so it is suffixed with everything
    # else, survives the merge as NA wherever the term is absent, and so records
    # "this dataset reported this term" without sniffing for non-NA values.
    enrichment_results[[i]]$.rc_contrib <- TRUE

    # Get renamed columns (suffixed by index)
    colnames(enrichment_results[[i]]) <- format_colnames(colnames(enrichment_results[[i]]))
    all_colnames <- colnames(enrichment_results[[i]])
    nonterm_cols <- all_colnames[all_colnames != 'Term']
    all_colnames[all_colnames != 'Term'] <- paste(nonterm_cols, i, sep=SEP)

    colnames(enrichment_results[[i]]) <- all_colnames
  }

  # C6: duplicate 'Term' values within an input frame cross-join.
  #
  # base::merge(by = 'Term', all = TRUE) is a SQL-style join, so a term carried
  # k times in one frame and j times in another emits k * j merged rows -- and
  # each of those rows pastes together a gene set that appears in NO input.
  # Those fabricated terms then enter the distance matrix, get clustered and get
  # counted, with nothing on screen to say so.
  #
  # The join is not itself wrong; its silence is. So this WARNS rather than
  # erroring or de-duplicating:
  #   - erroring would reject input the rest of the package is deliberately
  #     built to handle -- SPEC-RC-011's canonical sort key is (fold(Term),
  #     Term, GeneID) precisely because Term alone is not a total order once
  #     duplicates exist;
  #   - de-duplicating would silently discard the caller's rows and would also
  #     change the plotting layer, which indexes the distance matrix by term
  #     NAME and already returns only the first match.
  # Bundled data carries 0 duplicated terms, so this guard moves no baseline.
  dup_counts <- vapply(
    enrichment_results,
    function(df) if ("Term" %in% names(df)) sum(duplicated(df[["Term"]])) else 0L,
    integer(1))
  if (any(dup_counts > 0L)) {
    offenders <- which(dup_counts > 0L)
    # Rows base::merge(all = TRUE) will emit: per term, the product of its
    # per-frame occurrence counts, absent frames contributing a factor of 1.
    per_frame <- lapply(enrichment_results, function(df)
      if ("Term" %in% names(df)) table(df[["Term"]]) else table(character(0)))
    all_terms <- unique(unlist(lapply(per_frame, names), use.names = FALSE))
    projected <- sum(vapply(all_terms, function(tm)
      prod(vapply(per_frame, function(tb)
        if (tm %in% names(tb)) as.integer(tb[[tm]]) else 1L, integer(1))),
      numeric(1)))
    dup_terms <- unique(unlist(lapply(enrichment_results[offenders], function(df)
      df[["Term"]][duplicated(df[["Term"]])]), use.names = FALSE))
    # Base indexing rather than utils::head(): this file imports nothing from
    # utils, and one convenience call is not worth an import or a check NOTE.
    sample_terms <- dup_terms[seq_len(min(3L, length(dup_terms)))]
    warning(structure(
      class = c("richCluster_duplicate_terms", "warning", "condition"),
      list(message = sprintf(
             paste0("Dataset(s) %s carry duplicated 'Term' values (%s duplicate row(s); ",
                    "e.g. %s). Merging joins on 'Term', so these cross-join: the merged ",
                    "frame will have %d rows for %d distinct terms, and the extra rows ",
                    "carry gene sets that appear in no input. De-duplicate 'Term' per ",
                    "dataset if that is not what you want."),
             paste(offenders, collapse = ", "),
             paste(dup_counts[offenders], collapse = ", "),
             paste(sprintf("'%s'", sample_terms), collapse = ", "),
             as.integer(projected), length(all_terms)),
           call = NULL)))
  }

  # Initialize merged_gs with first geneset
  merged_gs <- enrichment_results[[1]]

  # Merge the rest of the enrichment_results.
  # The length > 1 guard is required, not decorative: `2:1` iterates c(2, 1) in
  # R, so an unguarded loop would merge a one-element list against itself.
  if (length(enrichment_results) > 1) {
    for (i in 2:length(enrichment_results)) {
      merged_gs <- base::merge(merged_gs, enrichment_results[[i]], by='Term', all=TRUE)
    }
  }

  # C12: order the rows canonically, on BOTH paths.
  #
  # base::merge(sort = TRUE) orders its result by 'Term' under LC_COLLATE, and a
  # single dataset never reaches merge at all, so the row order used to depend on
  # the session's locale.  Every index this package exports is a position in this
  # frame -- the distance matrix, TermIndices, and the seed walk david_cluster()
  # follows -- so a locale change moved results: the same data returned 532, 533
  # or 535 DAVID clusters depending only on the collation (measured 2026-09-05 on
  # Linux, macOS and Windows).  cluster() was unaffected because C9 made
  # everything downstream of the seed step canonical; david_cluster() never had
  # that treatment.
  #
  # method = "radix" sorts in the C locale by definition, so this order is the
  # same in every locale and on every platform.  The sort is stable, so rows
  # sharing a Term keep the order the join gave them.
  merged_gs <- merged_gs[order(merged_gs$Term, method = "radix"), , drop = FALSE]
  rownames(merged_gs) <- NULL

  # For each row in merged_gs, combine unique GeneID elements.
  # Same available_* shape the Pvalue / Padj blocks use, plus an explicit guard:
  # indexing the full geneid_cols unconditionally yielded 'undefined columns
  # selected', which names nothing the caller can act on.
  geneid_cols <- paste("GeneID", seq_along(enrichment_results), sep=SEP)
  available_geneid_cols <- geneid_cols[geneid_cols %in% colnames(merged_gs)]
  if (length(available_geneid_cols) == 0L) {
    stop(structure(
      class = c("richCluster_missing_geneid", "error", "condition"),
      list(message = sprintf(
             "No GeneID column found. Expected one of: %s. Columns present: %s.",
             paste(geneid_cols, collapse = ", "),
             paste(colnames(merged_gs), collapse = ", ")),
           call = NULL)))
  }
  # Dedup at the gene level, not the string level: unique() over whole
  # delimited strings left "A,B,C" + "B,C,D" as "A,B,C,B,C,D". Splitting
  # first, then unique(), keeps first appearance -- columns are scanned in
  # ascending index order and tokens left to right, which is the required order.
  merged_gs$GeneID <- apply(merged_gs[, available_geneid_cols, drop = FALSE], 1, function(x) {
    toks <- unlist(strsplit(stats::na.omit(x), gene_delim, fixed = TRUE), use.names = FALSE)
    paste(unique(toks[nzchar(toks)]), collapse = gene_delim)
  })


  # OD-7.  na.omit() above already drops absent columns, so a term whose gene
  # list is missing in EVERY contributing dataset lands here as "" -- the empty
  # set, which is the decided semantics.  Announce it: a term that silently
  # cannot cluster is the shape a user cannot debug.
  rc_warn_empty_gene_lists(sum(!nzchar(merged_gs$GeneID)), "merge_enrichment_results()")

  # Average the value columns across all enrichment_results
  # Avg Pvalue
  pvalue_cols <- paste("Pvalue", seq_along(enrichment_results), sep=SEP)
  available_pvalue_cols <- pvalue_cols[pvalue_cols %in% colnames(merged_gs)]
  # The length-1 case is split out: merged_gs[, <one col>] drops to a vector and
  # rowMeans() then fails with "'x' must be an array of at least two dimensions".
  # The > 1 branch is byte-identical to 1.0.2, so the demo-pair arithmetic is
  # unchanged.
  if (length(available_pvalue_cols) == 1L) {
    merged_gs$Pvalue <- merged_gs[[available_pvalue_cols]]
  } else if (length(available_pvalue_cols) > 1L) {
    merged_gs$Pvalue <- rowMeans(merged_gs[, available_pvalue_cols], na.rm = TRUE)
  }
  
  # Avg Padj
  padj_cols <- paste("Padj", seq_along(enrichment_results), sep=SEP)
  available_padj_cols <- padj_cols[padj_cols %in% colnames(merged_gs)]
  if (length(available_padj_cols) == 1L) {
    merged_gs$Padj <- merged_gs[[available_padj_cols]]
  } else if (length(available_padj_cols) > 1L) {
    merged_gs$Padj <- rowMeans(merged_gs[, available_padj_cols], na.rm = TRUE)
  }

  # Count how many input datasets reported each term, then drop the markers.
  # Assigned after Padj so 'DatasetCount' is the last column: every pre-existing
  # column keeps its positional index, which is what makes the addition safe for
  # 1.0.2 callers that index by position.
  contrib_cols <- paste(".rc_contrib", seq_along(enrichment_results), sep=SEP)
  merged_gs$DatasetCount <- as.integer(
    rowSums(!is.na(merged_gs[, contrib_cols, drop = FALSE]))
  )
  merged_gs <- merged_gs[, setdiff(colnames(merged_gs), contrib_cols), drop = FALSE]

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
    Pvalue = c("pvalue", "pval", "p-value"),
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
  mapped <- sapply(colnames, get_good_name)

  # Leftmost wins: a later column mapping onto an already-claimed standardised
  # name keeps its original name, and the collision is reported. Without this,
  # two source columns can collapse onto one name and base::merge then mangles
  # the duplicate.
  for (target in names(mappings)) {
    hits <- which(mapped == target)
    if (length(hits) > 1L) {
      warning(structure(
        class = c("richCluster_colname_collision", "warning", "condition"),
        list(message = sprintf(
               "Columns %s all standardise to '%s'; keeping '%s' and leaving the rest unchanged.",
               paste(sprintf("'%s'", colnames[hits]), collapse = ", "), target, colnames[hits[1]]),
             call = NULL)))
      mapped[hits[-1L]] <- colnames[hits[-1L]]
    }
  }

  return(mapped)
}
