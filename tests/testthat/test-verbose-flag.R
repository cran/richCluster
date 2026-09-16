# tests/testthat/test-verbose-flag.R
# The C++ core narrates its progress to stdout.  Until the `verbose` flag, that
# narration was UNCONDITIONAL: 14 Rcpp::Rcout sites across src/RichCluster.cpp
# (11) and src/DavidClustering.cpp (3), with no way to switch them off.
#
# WHY THIS IS A RELEASE BLOCKER, not a style preference.  CRAN bounced 1.0.2
# over a single print() (recorded in cran-comments.md), and this is the same
# objection at fourteen times the scale.  Measured before the fix:
# suppressMessages() still leaked 23 lines, because Rcout writes to stdout while
# suppressMessages() only silences the message() condition on stderr.  The
# volume also scales with the data -- one line per merge iteration -- so a large
# run buries the caller's own output.
#
# The package itself already conceded the defect in two places rather than fix
# it: test-shipped-numeric-regression.R wraps every call in capture.output(),
# and test-cpp-entry-guards.R sinks to nullfile().  Both workarounds are left in
# place; they are harmless, and this file is the guard that makes them
# unnecessary rather than mandatory.
#
# Inputs are inline, so this file has no fixture dependency and ships in the
# tarball.  Do NOT rename to test-src-*: that prefix is tarball-excluded.

# Three terms; the first two share two of three genes, the third shares none.
# The actual partition is irrelevant here -- this file asserts on OUTPUT, never
# on cluster membership, so it cannot go stale when the numerics move.
rc_v_terms <- c("T1", "T2", "T3")
rc_v_genes <- c("g1,g2,g3", "g1,g2,g4", "g5,g6,g7")

rc_v_frame <- function(seed) {
  data.frame(
    Term   = rc_v_terms,
    GeneID = rc_v_genes,
    Pvalue = c(0.001, 0.002, 0.003) * seed,
    Padj   = c(0.010, 0.020, 0.030) * seed,
    stringsAsFactors = FALSE
  )
}

# capture.output() takes R's stdout connection, which is exactly where
# Rcpp::Rcout writes -- so a silent run returns character(0) and a narrating one
# returns its lines.  force() keeps the promise from escaping the capture.
#
# invisible() is LOAD-BEARING, not decoration: capture.output() auto-prints a
# VISIBLE result, and every function under test returns a list.  Without it the
# captured text is that list's print method -- measured: 12 lines opening
# `$distance_matrix` -- and the silence assertions fail against output the C++
# never wrote.  invisible() marks the value non-printing so what is captured is
# only what the run actually emitted.
rc_v_stdout <- function(expr) {
  utils::capture.output(invisible(force(expr)))
}

test_that("runRichCluster() is silent by default", {
  out <- rc_v_stdout(
    runRichCluster(rc_v_terms, rc_v_genes, "kappa", 0.5, "average", 0.5)
  )
  # Measured before the fix: 14 lines, opening with "Starting richCluster...".
  expect_identical(out, character(0))
})

test_that("runRichCluster(verbose = TRUE) narrates", {
  out <- rc_v_stdout(
    runRichCluster(rc_v_terms, rc_v_genes, "kappa", 0.5, "average", 0.5,
                   verbose = TRUE)
  )
  expect_gt(length(out), 0)
  # The opening banner is the one line every run emits regardless of data, so
  # it is the only content this file pins.  Asserting on the merge-iteration
  # lines would couple an output test to the clustering numerics.
  expect_true(any(grepl("Starting richCluster", out, fixed = TRUE)))
})

test_that("cluster() is silent by default and threads verbose through", {
  input <- list(rc_v_frame(1), rc_v_frame(2))

  quiet <- rc_v_stdout(
    cluster(input, min_terms = 1, min_value = 1)
  )
  expect_identical(quiet, character(0))

  loud <- rc_v_stdout(
    cluster(input, min_terms = 1, min_value = 1, verbose = TRUE)
  )
  expect_gt(length(loud), 0)
  expect_true(any(grepl("Starting richCluster", loud, fixed = TRUE)))
})

test_that("runDavidClustering() is silent by default", {
  out <- rc_v_stdout(
    richCluster:::runDavidClustering(rc_v_terms, rc_v_genes, 0.5, 1, 1, 0.5)
  )
  # Measured before the fix: 3 lines -- "Calculating kappa scores...",
  # "Finding initial seeds...", "Merging seeds...".
  expect_identical(out, character(0))
})

test_that("runDavidClustering(verbose = TRUE) narrates", {
  out <- rc_v_stdout(
    richCluster:::runDavidClustering(rc_v_terms, rc_v_genes, 0.5, 1, 1, 0.5, verbose = TRUE)
  )
  expect_gt(length(out), 0)
  expect_true(any(grepl("Calculating kappa scores", out, fixed = TRUE)))
})

test_that("runDavidClusteringWithKappa() honours verbose too", {
  # The second DAVID entry point runs the same pipeline; a flag threaded to one
  # and not the other would leave a silent path and a loud one for the same
  # computation.
  quiet <- rc_v_stdout(
    richCluster:::runDavidClusteringWithKappa(rc_v_terms, rc_v_genes, 0.5, 1, 1, 0.5)
  )
  expect_identical(quiet, character(0))

  loud <- rc_v_stdout(
    richCluster:::runDavidClusteringWithKappa(rc_v_terms, rc_v_genes, 0.5, 1, 1, 0.5,
                                verbose = TRUE)
  )
  expect_gt(length(loud), 0)
})

test_that("david_cluster() is silent by default and threads verbose through", {
  input <- list(rc_v_frame(1), rc_v_frame(2))

  quiet <- rc_v_stdout(
    david_cluster(input, similarity_threshold = 0.5,
                  initial_group_membership = 1, final_group_membership = 1,
                  multiple_linkage_threshold = 0.5)
  )
  expect_identical(quiet, character(0))

  loud <- rc_v_stdout(
    david_cluster(input, similarity_threshold = 0.5,
                  initial_group_membership = 1, final_group_membership = 1,
                  multiple_linkage_threshold = 0.5, verbose = TRUE)
  )
  expect_gt(length(loud), 0)
})

test_that("no ungated Rcout survives in the C++ sources", {
  # A source-level guard, so a future edit cannot reintroduce an unconditional
  # print without turning this file red.  src/ is absent from the .Rcheck tree,
  # so this skips under R CMD check and runs during development.
  src <- file.path("..", "..", "src")
  skip_if_not(dir.exists(src), "src/ not present (installed tests)")

  files <- list.files(src, pattern = "\\.cpp$", full.names = TRUE)
  # RcppExports.cpp carries the generated `Rcpp::Rcout = Rcpp::Rcpp_cout_get()`
  # declaration, which is a definition rather than a print site.
  files <- files[basename(files) != "RcppExports.cpp"]

  offenders <- character(0)
  for (f in files) {
    lines <- readLines(f, warn = FALSE)
    hits <- grep("Rcpp::Rcout", lines, fixed = TRUE)
    for (h in hits) {
      # Accept a site only when `verbose` guards it, either on the same line
      # (`if (verbose) Rcpp::Rcout << ...`) or on the line above.
      window <- paste(lines[max(1, h - 1):h], collapse = " ")
      if (!grepl("verbose", window, fixed = TRUE)) {
        offenders <- c(offenders, sprintf("%s:%d", basename(f), h))
      }
    }
  }
  expect_identical(offenders, character(0))
})
