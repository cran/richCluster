# NULL placeholder for roxygen namespace declarations
NULL

# OD-7 (author decision, 2026-08-30): a MISSING gene list is an EMPTY SET.
#
# An empty set is similar to nothing, so a term carrying one cannot join a
# cluster and min_terms drops it naturally -- no row is removed from the
# caller's data and no count changes.  The one thing that must not happen is
# for that to be silent, hence this warning.
#
# It is classed rather than bare so a caller can catch it selectively, matching
# richCluster_duplicate_terms / richCluster_missing_geneid / the rest.
rc_warn_empty_gene_lists <- function(n, where) {
  n <- as.integer(n)
  if (is.na(n) || n <= 0L) return(invisible(FALSE))
  warning(structure(
    class = c("richCluster_empty_gene_list", "warning", "condition"),
    list(message = sprintf(
           paste0("%d term(s) have an empty gene list in %s. An empty gene list ",
                  "is treated as the empty set: such terms are similar to nothing ",
                  "and will not join any cluster."),
           n, where),
         call = NULL)))
  invisible(TRUE)
}
