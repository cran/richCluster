# tests/testthat/helper-oracle.R
# ===========================================================================
# T3-02 (WAVE 0) -- AN INDEPENDENT PURE-R COHEN'S KAPPA ORACLE.
#
# The single most load-bearing artifact in the plan.  It is the only
# verification of T0-01 that depends on neither a hand-copied constant, nor a
# defective self-oracle, nor an undefined fixture.
#
# ---------------------------------------------------------------------------
# HOW THIS FILE WAS WRITTEN, AND WHY IT MATTERS
# ---------------------------------------------------------------------------
# EVERY LINE BELOW IS DERIVED FROM THE MATHEMATICAL DEFINITION OF COHEN'S
# KAPPA OVER A 2x2 AGREEMENT TABLE.  Nothing was ported, transcribed or read
# from src/DistanceMetric.cpp or src/DavidClustering.cpp.  An oracle derived
# from the code it exists to check is worthless -- it reproduces that code's
# defects and reports agreement.  This is the one place in the release where
# independence of derivation matters more than agreement, and it is free.
#
# The plan specifies the CONTRACT and deliberately says nothing about
# vectorisation, tcrossprod, loops or data structures (critique M7, m5).  The
# shapes chosen here are this team's own; the other team's oracle is expected
# to look different and to agree to 1e-12 (V1).
#
# ---------------------------------------------------------------------------
# THE CONTRACT (PLAN_V2_T3.md T3-02)
# ---------------------------------------------------------------------------
# Build a terms x genes BINARY MEMBERSHIP MATRIX over the fixture.  For each
# term pair (A, B) form the 2x2 agreement table with
#
#     N = the number of DISTINCT GENES IN THE DATASET
#
# -- not |A u B|, not the number of gene tokens, not the number of terms.  Then
#
#     both    = |A n B|
#     neither = N - |A| - |B| + |A n B|
#     oab     = (both + neither) / N
#     aab     = (|A||B| + (N - |A|)(N - |B|)) / N^2
#     kappa   = (oab - aab) / (1 - aab),      kappa = 1 when aab == 1
#
# ---------------------------------------------------------------------------
# OD-5 IS SETTLED, AND IT CHANGES THIS CONTRACT.  READ THIS BEFORE EDITING.
# ---------------------------------------------------------------------------
# OD-5 `OD-DISJOINT-KAPPA` was settled by the author on 2026-08-09:
#
#     FLOOR ALL NEGATIVE KAPPA VALUES AT 0.  The exported similarity scale
#     is [0, 1].
#
# This is a THIRD option, not the remove/retain binary the plan's OD-5 bullet
# offered, and BOTH of those are superseded.  T3-02's contract as printed in
# PLAN_V2_T3.md encodes the *removal* branch (no clamp at all) and is
# SUPERSEDED ON THAT POINT.  The clamp is not optional and it is not confined
# to the disjoint case:
#
#     kappa_reported = max(kappa, 0)      for every pair, disjoint or not
#
# On the 580-term anchor fixture this clamps 26159 upper-triangle pairs (mean
# -0.013692, most negative -0.107613) and moves the upper-triangle sum from
# 14239.673391 to 14597.851964.  An oracle WITHOUT the clamp disagrees with
# corrected C++ on all 26159 of them and T0-01 V1's 1e-12 gate then FAILS ON
# CORRECT CODE -- which is precisely the failure mode OD-5 was gated to
# prevent.
#
# The clamp is a SEPARATE, LAST step, applied to the finished chance-corrected
# value.  It is deliberately NOT folded into the closed form:
#   * `clamp = FALSE` recovers the raw chance-corrected kappa, which is what
#     the UNCLAMPED INDEPENDENCE CHECK against DavidClustering needs (V2), and
#   * keeping the two separable is what lets a reader see that the clamp is a
#     reporting-scale decision and not part of the mathematics.
#
# ---------------------------------------------------------------------------
# THE DIAGONAL
# ---------------------------------------------------------------------------
# This oracle is mathematically correct and therefore reports kappa(A, A) == 1
# on the diagonal.  The DAVID kappa matrix exposed by T3-04a carries 0 there
# (calculateKappaScores() loops `for (j = i + 1; ...)` and never writes the
# diagonal), and the stock DistanceMetric matrix carries the -99
# SAME_TERM_DISTANCE sentinel.  EVERY comparison against either C++ matrix is
# therefore over the UPPER TRIANGLE ONLY.  That is a property of those two
# matrices, not a defect of this oracle, and it must not be "fixed" here.
# ===========================================================================


# --- The closed form, on counts --------------------------------------------

#' Cohen's kappa for two gene sets, from their cardinalities.
#'
#' The whole of the mathematics lives here.  Everything else in this file is
#' bookkeeping that produces the four counts.
#'
#' Vectorised over `n_a`, `n_b`, `n_common` and `n_universe` by ordinary R
#' recycling, so it serves both the scalar hand cases (V3), the 200000-draw
#' invariant sweep (V4) and the full 580 x 580 matrix (V1, V2) unchanged.
#'
#' @param n_a,n_b sizes of the two gene sets, |A| and |B|.
#' @param n_common size of their intersection, |A n B|.
#' @param n_universe N -- the number of DISTINCT GENES IN THE DATASET.  This is
#'   the parameter every naive implementation gets wrong: it is a property of
#'   the whole input, not of the pair.
#' @param clamp apply the settled OD-5 floor.  TRUE is the shipping contract;
#'   FALSE returns the raw chance-corrected value and exists for the unclamped
#'   independence check against DavidClustering (V2) and for exhibiting what
#'   the clamp actually changed.
#' @return kappa, in [0, 1] when `clamp = TRUE` and in [-1, 1] otherwise.
rc_oracle_kappa <- function(n_a, n_b, n_common, n_universe, clamp = TRUE) {

  N <- n_universe

  if (any(N <= 0)) {
    stop("rc_oracle_kappa(): n_universe must be positive -- kappa is undefined ",
         "over an empty gene universe.", call. = FALSE)
  }
  if (any(n_a < 0 | n_b < 0 | n_common < 0)) {
    stop("rc_oracle_kappa(): set sizes must be non-negative.", call. = FALSE)
  }
  if (any(n_a > N | n_b > N)) {
    stop("rc_oracle_kappa(): |A| and |B| cannot exceed the gene universe N.",
         call. = FALSE)
  }
  if (any(n_common > pmin(n_a, n_b))) {
    stop("rc_oracle_kappa(): |A n B| cannot exceed min(|A|, |B|).", call. = FALSE)
  }
  if (any(n_a + n_b - n_common > N)) {
    stop("rc_oracle_kappa(): |A u B| cannot exceed the gene universe N.",
         call. = FALSE)
  }

  # The 2x2 agreement table.  Rows: gene in A / not in A.  Columns: gene in B /
  # not in B.  Every one of the N genes falls in exactly one cell.
  #
  #                  in B          not in B
  #   in A           both          n_a - both
  #   not in A       n_b - both    neither
  #
  both    <- n_common
  neither <- N - n_a - n_b + n_common          # N - |A u B|

  # Observed agreement: the two diagonal cells, as a proportion of N.
  oab <- (both + neither) / N

  # Chance agreement: the two marginals multiply independently.  P(both in) is
  # (|A|/N)(|B|/N) and P(both out) is ((N-|A|)/N)((N-|B|)/N); written over the
  # common denominator N^2 so the numerator stays an exact integer.
  aab <- (n_a * n_b + (N - n_a) * (N - n_b)) / (N * N)

  kappa <- (oab - aab) / (1 - aab)

  # aab == 1 happens only in the degenerate cases |A| = |B| = N (both terms
  # span the whole universe) and |A| = |B| = 0.  Both are the maximally
  # similar case -- the two sets agree on every one of the N genes -- so kappa
  # is 1 there.  The quotient above is 0/0 and must be replaced, not rescued.
  # (`>= 1` rather than `== 1` only as a rounding guard: aab <= 1 always.)
  degenerate <- aab >= 1
  kappa[degenerate] <- 1

  # --- The settled OD-5 clamp.  Last, separate, and applied to ALL negatives.
  if (isTRUE(clamp)) kappa <- pmax(kappa, 0)

  kappa
}


#' The Dice coefficient, 2|A n B| / (|A| + |B|).
#'
#' Present only so V4 can assert the kappa -> Dice limit as N -> inf, which is
#' a property of the oracle's own arithmetic.  Independent of the package's
#' `distance_metric = "dice"` (D-02), which this file does not touch.
rc_oracle_dice <- function(n_a, n_b, n_common) {
  denom <- n_a + n_b
  out <- 2 * n_common / denom
  out[denom == 0] <- 1          # two empty sets are identical
  out
}


# --- The terms x genes binary membership matrix ----------------------------

#' Split a vector of delimited gene-ID strings into character vectors of
#' DISTINCT gene names.
#'
#' `unique()` is load-bearing, not defensive tidying.  The anchor fixture DOES
#' repeat gene names inside a single term's string (measured: the largest term
#' carries 2328 tokens over a 1819-gene universe), and a membership MATRIX is
#' binary by definition -- a gene is in a term or it is not.  Counting tokens
#' instead of members would give |A| > N and put the 2x2 table outside its own
#' universe.
rc_oracle_gene_sets <- function(gene_ids, delim = ",") {
  if (anyNA(gene_ids)) {
    stop("rc_oracle_gene_sets(): NA gene-ID strings have no membership ",
         "interpretation; resolve the NA policy before calling the oracle.",
         call. = FALSE)
  }
  lapply(strsplit(as.character(gene_ids), delim, fixed = TRUE), unique)
}


#' The gene universe: every distinct gene name appearing anywhere in the input.
#'
#' This is the N of the contract.  Sorted so that the membership matrix's
#' column order is a deterministic function of the input alone.
rc_oracle_universe <- function(gene_ids, delim = ",") {
  sort(unique(unlist(rc_oracle_gene_sets(gene_ids, delim), use.names = FALSE)))
}


#' The terms x genes binary membership matrix.
#'
#' @return a numeric matrix of 0s and 1s, one row per term (in input order) and
#'   one column per distinct gene (in sorted order).  Numeric rather than
#'   logical so that the pairwise intersection counts below are an ordinary
#'   matrix product.
rc_oracle_membership <- function(gene_ids, delim = ",") {
  sets  <- rc_oracle_gene_sets(gene_ids, delim)
  genes <- sort(unique(unlist(sets, use.names = FALSE)))
  M <- matrix(0, nrow = length(sets), ncol = length(genes),
              dimnames = list(NULL, genes))
  for (i in seq_along(sets)) M[i, match(sets[[i]], genes)] <- 1
  M
}


#' The full pairwise kappa matrix over a vector of gene-ID strings.
#'
#' @param gene_ids character vector of delimited gene-ID strings, one per term.
#' @param clamp the settled OD-5 floor; see the file header.
#' @param delim the gene delimiter.
#' @return an n_terms x n_terms numeric matrix, symmetric, with 1 on the
#'   diagonal (kappa(A, A) == 1).  Un-dimnamed, indexed in the order of
#'   `gene_ids`, so it lines up positionally with the C++ matrices.
#'
#' The diagonal is 1 BY THE MATHEMATICS.  The DAVID matrix carries 0 there and
#' the stock DistanceMetric matrix carries -99; comparisons against either are
#' over the upper triangle only.  See the file header.
rc_oracle_kappa_matrix <- function(gene_ids, clamp = TRUE, delim = ",") {

  M <- rc_oracle_membership(gene_ids, delim)
  n <- nrow(M)
  N <- ncol(M)

  # Pairwise intersection counts.  M[i, ] and M[j, ] are 0/1 indicator rows, so
  # their inner product counts exactly the genes present in both terms; the
  # matrix product does all n^2 of them at once.  Every partial sum is a small
  # integer, so this is exact in double precision.
  common <- tcrossprod(M)

  sizes <- rowSums(M)                       # |A| for each term
  n_a <- matrix(sizes, nrow = n, ncol = n)              # varies down rows
  n_b <- matrix(sizes, nrow = n, ncol = n, byrow = TRUE) # varies across columns

  K <- rc_oracle_kappa(n_a, n_b, common, N, clamp = clamp)
  dim(K) <- c(n, n)
  dimnames(K) <- NULL
  K
}


# --- V1 cross-team exchange ------------------------------------------------
#
# V1 requires the TWO TEAMS' oracles to agree to 1e-12 on the fixture BEFORE
# either team edits any package code.  A disagreement found now costs nothing;
# found after five files have moved it costs the release.
#
# The exchanged object is deliberately NOT written through rc_emit(): the
# T3-01(d) registry is a registry of PACKAGE ENTRY POINTS (runRichCluster(),
# david_cluster(), merge_enrichment_results(), ...) and rc_artifact_registry()
# errors on any tag it does not define -- correctly, since an unregistered tag
# silently becoming a baseline is the defect that registry exists to prevent.
# The oracle is not a package entry point and inventing a registry family for
# it would be an edit to T3-01's file.  It is instead written self-describing
# under its own filename, carrying the same information the registry carries:
# what produced it, on what input, under which OD-5 branch.

#' Path of this team's exchanged oracle matrix.
rc_oracle_crossteam_file <- function() {
  file.path(rc_fixture_dir(), "oracle-kappa-580.rds")
}

#' Path at which the OTHER team's oracle matrix is expected.
#'
#' Overridable with RC_CROSSTEAM_ORACLE so the file can be dropped anywhere.
rc_oracle_counterpart_file <- function() {
  env <- Sys.getenv("RC_CROSSTEAM_ORACLE", unset = "")
  if (nzchar(env)) return(env)
  file.path(rc_fixture_dir(), "oracle-kappa-580-crossteam.rds")
}

#' Build the self-describing exchange record for a kappa matrix.
#'
#' Carries its own provenance so that a matrix which turns up in the wrong
#' place cannot be silently mistaken for a baseline: the input it was computed
#' from is identified by term count, universe size and a checksum of the frozen
#' fixture's own gene strings, and the OD-5 branch is recorded explicitly
#' because an unclamped and a clamped oracle are different objects.
rc_oracle_record <- function(kappa, gene_ids, clamp) {
  list(
    produced_by  = "T3-02 pure-R Cohen's kappa oracle (rc_oracle_kappa_matrix)",
    input        = "T3-01 frozen anchor fixture, fx-input.rds $GeneID",
    n_terms      = length(gene_ids),
    n_universe   = length(rc_oracle_universe(gene_ids)),
    input_digest = sum(utf8ToInt(paste(gene_ids, collapse = "|"))),
    # SPEC-RC-005, finding A4: name the expression the digest was taken over, so
    # a value that does not reproduce can be diagnosed rather than merely
    # disputed.  A checksum without its input expression is not a checksum.
    input_digest_expr = "sum(utf8ToInt(paste(gene_ids, collapse = \"|\")))",
    od5_branch   = if (isTRUE(clamp)) "clamp ALL negatives at 0 (SETTLED)"
                   else "unclamped raw chance-corrected kappa",
    clamp        = clamp,
    diagonal     = "1 (kappa(A, A) == 1); compare upper triangle only",
    upper_tri_sum = sum(kappa[upper.tri(kappa)]),
    kappa        = kappa,
    created_at   = Sys.time(),
    r_version    = R.version.string
  )
}

#' Write this team's oracle matrix for exchange with the other team.
#'
#' Refuses to overwrite a differing record unless `overwrite = TRUE`, for the
#' same reason rc_freeze_fixture() does: silently re-emitting after the oracle
#' has been edited destroys the evidence that the two teams ever agreed.
rc_oracle_emit_crossteam <- function(gene_ids, clamp = TRUE, overwrite = FALSE) {
  dir <- rc_fixture_dir()
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  path <- rc_oracle_crossteam_file()
  rec  <- rc_oracle_record(rc_oracle_kappa_matrix(gene_ids, clamp = clamp),
                           gene_ids, clamp)
  if (file.exists(path) && !isTRUE(overwrite)) {
    old <- readRDS(path)
    if (!identical(old$kappa, rec$kappa)) {
      stop(sprintf(paste0(
        "Refusing to overwrite '%s': the recomputed oracle differs from the ",
        "emitted one.  Either the oracle changed or the fixture did; find out ",
        "which before overwriting."), path), call. = FALSE)
    }
    return(invisible(path))
  }
  saveRDS(rec, path)
  invisible(path)
}

# --- SPEC-RC-004: generate the exchange record when it is not stored -------
#
# This generator calls NO richCluster function -- bare or qualified -- and NO
# compiled entry point.  That independence is the whole value of the oracle: an
# oracle derived from the code it checks reproduces that code's defects and
# reports agreement.  AC-008 enforces this by intersecting codetools::findGlobals
# against the package namespace, so a BARE call is caught too, and by a textual
# arm that also catches string dispatch.
#
# rc_assert_invariant() is defined in helper-fixture.R and is a TEST HELPER, not
# a package function, so naming it here does not violate REQ-RC004-005.

#' Generate the oracle exchange record from the anchor fixture.
rc_oracle_generate <- function(clamp = TRUE) {
  fx <- rc_fixture_frozen()
  K  <- rc_oracle_kappa_matrix(fx$GeneID, clamp = clamp)

  rc_assert_invariant("O1 dim(K)", dim(K), c(580L, 580L))
  ut <- K[upper.tri(K)]
  if (isTRUE(clamp)) {
    rc_assert_invariant("O2 clamped upper-triangle sum", sum(ut), 14597.85196436438, tolerance = 1e-9)
    rc_assert_invariant("O5 clamped range [0, 1]", sum(ut < 0 | ut > 1), 0L)
    rc_assert_invariant("O6 diagonal", unique(diag(K)), 1)
  } else {
    rc_assert_invariant("O3 raw upper-triangle sum", sum(ut), 14239.67339060683, tolerance = 1e-9)
    rc_assert_invariant("O4 pairs the OD-5 clamp moves", sum(ut < 0), 26159L)
  }

  rc_oracle_record(K, fx$GeneID, clamp)
}

#' Load the oracle exchange record: stored file when present, else generate.
#' REQ-RC004-004 / REQ-RC004-012.
rc_oracle_crossteam <- function() {
  path <- rc_oracle_crossteam_file()
  if (file.exists(path)) return(readRDS(path))
  if (is.null(.rc_memo$oracle)) .rc_memo$oracle <- rc_oracle_generate(clamp = TRUE)
  .rc_memo$oracle
}
