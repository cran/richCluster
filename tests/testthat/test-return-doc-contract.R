# Guards the @return blocks of cluster() and runRichCluster() against silent
# drift from what the functions actually return.
#
# The documented blocks had drifted badly before this file existed:
# runRichCluster()'s \value named a `linkage_tree` element that has never
# existed anywhere in src/, and cluster()'s named a `clusters` element that
# does not exist while omitting `all_clusters`, `final_clusters` and
# `cluster_df` -- the last of which is the frame every plotting function
# consumes, so a user reading ?cluster could not discover the package's main
# output.  Prose cannot be checked by the compiler; this file checks it.
#
# The documented set is read from the INSTALLED .Rd (tools::Rd_db), not from a
# constant duplicated here, so the assertion binds the shipped documentation --
# the thing a user actually reads -- rather than a copy of it.


# --- Pull the \item tag names out of a \describe block inside \value --------
#
# Rd is parsed, not grepped: an Rd_db element is a nested list whose nodes
# carry an "Rd_tag" attribute.  \value holds a \describe, whose \item nodes
# each hold list(name, body); we want the names, in document order.
rd_value_item_names <- function(rd_file) {
  db <- tools::Rd_db("richCluster")
  if (!rd_file %in% names(db)) {
    stop("man/", rd_file, " is not in the installed Rd db.", call. = FALSE)
  }
  rd <- db[[rd_file]]

  tag_of <- function(x) {
    a <- attr(x, "Rd_tag")
    if (is.null(a)) NA_character_ else a
  }

  value <- rd[vapply(rd, tag_of, character(1)) == "\\value"]
  if (length(value) != 1L) {
    stop(rd_file, " has ", length(value), " \\value sections; expected 1.",
         call. = FALSE)
  }
  value <- value[[1L]]

  describe <- value[vapply(value, tag_of, character(1)) == "\\describe"]
  if (length(describe) != 1L) {
    stop(rd_file, " \\value has ", length(describe), " \\describe blocks; ",
         "expected 1.  The element-by-element \\describe form is what makes ",
         "the documented names machine-checkable -- keep it.", call. = FALSE)
  }
  describe <- describe[[1L]]

  items <- describe[vapply(describe, tag_of, character(1)) == "\\item"]
  unname(vapply(items,
                function(it) paste(unlist(it[[1L]]), collapse = ""),
                character(1)))
}


# --- A fixture small enough to run in-line ---------------------------------
#
# cluster() needs only Term / GeneID / Pvalue / Padj; merge_enrichment_results()
# supplies the rest.  min_terms = 2 keeps filter_clusters() from emptying
# final_clusters, so cluster_df is a real frame rather than a degenerate one.
rc_doc_contract_frames <- function() {
  terms <- sprintf("T%02d", 1:8)
  genes <- c("g1,g2,g3", "g1,g2,g3", "g1,g2,g4", "g1,g2,g3",
             "g5,g6,g7", "g5,g6,g7", "g5,g6,g8", "g5,g6,g7")
  df <- data.frame(Term = terms, GeneID = genes,
                   Pvalue = rep(1e-6, length(terms)),
                   Padj   = rep(1e-6, length(terms)),
                   stringsAsFactors = FALSE)
  list(df, df)
}


test_that("runRichCluster() \\value documents exactly the elements returned", {
  documented <- rd_value_item_names("runRichCluster.Rd")

  utils::capture.output(
    actual <- runRichCluster(c("T1", "T2", "T3", "T4"),
                             c("a,b,c", "a,b,d", "x,y,z", "a,b,c"),
                             "kappa", 0.5, "average", 0.5, ",")
  )

  expect_identical(documented, names(actual))
})


test_that("runRichCluster()'s documented all_clusters columns are the real ones", {
  # The three-column shape is the half of the drift a names()-only check misses:
  # the block used to name Cluster and TermIndices and silently drop TermNames.
  utils::capture.output(
    actual <- runRichCluster(c("T1", "T2", "T3", "T4"),
                             c("a,b,c", "a,b,d", "x,y,z", "a,b,c"),
                             "kappa", 0.5, "average", 0.5, ",")
  )
  expect_identical(names(actual$all_clusters),
                   c("Cluster", "TermNames", "TermIndices"))
})


test_that("cluster() \\value documents exactly the elements returned", {
  documented <- rd_value_item_names("cluster.Rd")

  utils::capture.output(
    actual <- cluster(rc_doc_contract_frames(),
                      df_names = c("a", "b"), min_terms = 2)
  )

  expect_identical(documented, names(actual))
})


test_that("cluster() documents cluster_df, the frame the plotting functions take", {
  # Called out separately because this is the omission with the real user cost:
  # cluster_df is what cluster_hmap() / cluster_dot() / cluster_bar() consume,
  # and it was absent from ?cluster entirely.
  expect_true("cluster_df" %in% rd_value_item_names("cluster.Rd"))
})
