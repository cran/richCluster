# tests/testthat/test-t002-denominator.R
# ===========================================================================
# SPEC-RC-007 / T0-02 -- the kappa denominator is computed in double.
#
# src/DavidClustering.cpp:93 casts BOTH numerator products to double but left
# the denominator as `totalGeneCount * totalGeneCount` -- an int * int product.
# int32 max is 2,147,483,647, so the product overflows at N > 46,340: signed
# integer overflow, which is undefined behaviour, and `aab` becomes garbage.
#
# The two RED assertions below were both measured against the unfixed build
# before the change (2026-08-25):
#   N = 66,200 -> every pairwise kappa is 1.004032, i.e. OUT of [-1, 1]  (V1)
#   N = 48,200 -> every pairwise kappa is 0.9166301, in range but WRONG   (V2)
# V2 exists because the range check alone is not sufficient: between the
# overflow threshold and roughly N = 65,000 the wrapped denominator still
# yields an in-range value, so only an oracle comparison detects the defect.
#
# At the frozen fixture's N = 1819 the product is 3,308,761 -- four orders of
# magnitude below the ceiling -- so the cast provably changes no number there.
# That no-op claim is asserted by V3 and by AC-002-1 (the v102-david identity).
# ===========================================================================

# A synthetic term set whose distinct-gene universe is `n_terms * block +
# shared_n`.  Every term shares one common block, so kappa is non-degenerate.
rc_t002_synth <- function(n_terms, block, shared_n) {
  shared <- sprintf("S%06d", seq_len(shared_n))
  terms  <- sprintf("TERM_%02d", seq_len(n_terms))
  genes  <- vapply(seq_len(n_terms), function(i) {
    own <- sprintf("G%06d", ((i - 1) * block + 1):(i * block))
    paste(c(own, shared), collapse = ",")
  }, character(1))
  list(Term = terms, GeneID = genes,
       N = length(unique(unlist(strsplit(genes, ",")))))
}

rc_t002_kappa <- function(x) {
  out <- NULL
  utils::capture.output(
    out <- richCluster:::runDavidClusteringWithKappa(
      x$Term, x$GeneID, 0.5, 3, 3, 0.5)$kappa_matrix)
  out
}


test_that("V1: kappa stays within [-1, 1] when N exceeds the int32 threshold", {
  x <- rc_t002_synth(n_terms = 12, block = 5500, shared_n = 200)
  # The premise of the test: if this is not above 46,340 the generator is
  # broken, and a pass would be meaningless.  Fix the generator, not the test.
  expect_gt(x$N, 46340L)

  km <- rc_t002_kappa(x)
  ut <- km[upper.tri(km)]

  expect_false(anyNA(ut))
  expect_true(all(is.finite(ut)))
  expect_true(all(ut >= -1 & ut <= 1))
})


test_that("V2: kappa matches the double-precision oracle above the threshold", {
  # N = 48,200: past the overflow threshold, but the wrapped denominator still
  # produces an in-range value, so V1 alone cannot see the defect here.
  x <- rc_t002_synth(n_terms = 12, block = 4000, shared_n = 200)
  expect_gt(x$N, 46340L)

  km <- rc_t002_kappa(x)

  # Every pair has the same shape: |A| = |B| = block + shared, and the only
  # genes in common are the shared block.
  # Doubles throughout: the oracle's own N * N term overflows R's integer type
  # at exactly the same threshold the C++ defect does.
  n_a <- 4000 + 200
  expected <- rc_oracle_kappa(n_a = n_a, n_b = n_a, n_common = 200,
                              n_universe = as.numeric(x$N), clamp = FALSE)
  expect_equal(unique(as.numeric(round(km[upper.tri(km)], 12))),
               round(expected, 12))
})


test_that("V3: below the threshold the cast is a no-op", {
  # N = 4,200; the product is 17,640,000 -- two orders below the ceiling.  This
  # assertion is green both before and after the change, by construction: a
  # double represents every integer up to 2^53 exactly and the operand order is
  # unchanged.  It is here so a regression at ordinary scale is visible.
  x <- rc_t002_synth(n_terms = 8, block = 500, shared_n = 200)
  expect_lt(x$N, 46340L)

  km <- rc_t002_kappa(x)
  n_a <- 500 + 200
  expected <- rc_oracle_kappa(n_a = n_a, n_b = n_a, n_common = 200,
                              n_universe = as.numeric(x$N), clamp = FALSE)
  expect_equal(unique(as.numeric(round(km[upper.tri(km)], 12))),
               round(expected, 12))
})
