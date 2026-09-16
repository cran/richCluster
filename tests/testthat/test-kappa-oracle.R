# tests/testthat/test-kappa-oracle.R
# ===========================================================================
# T3-02 (WAVE 0) -- THE INDEPENDENT PURE-R COHEN'S KAPPA ORACLE.
#
# The oracle itself is helper-oracle.R; this file is its four verifications.
# It is the only verification of T0-01 that depends on neither a hand-copied
# constant, nor a defective self-oracle, nor an undefined fixture, and it is
# written BEFORE any numeric package code moves.
#
# INDEPENDENCE OF DERIVATION.  helper-oracle.R was written from the
# mathematical definition of Cohen's kappa over a 2x2 agreement table.  No
# line of src/DistanceMetric.cpp or src/DavidClustering.cpp was ported,
# transcribed or consulted.  An oracle derived from the code it exists to
# check reproduces that code's defects and reports agreement -- which is
# exactly what the NULL-ACTION CHECK at the bottom of this file is designed to
# expose.
#
# OD-5 IS SETTLED AND SUPERSEDES T3-02'S PRINTED CONTRACT.  The ruling is
# FLOOR ALL NEGATIVE KAPPA AT 0 -- a third option, not the remove/retain
# binary the plan's OD-5 bullet offered.  PLAN_V2_T3.md's T3-02 contract
# encodes the *removal* branch (no clamp) and is superseded on that point.
# The oracle applies the clamp by default; `clamp = FALSE` recovers the raw
# chance-corrected value, which is what V2's independence check needs.
#
# THE DIAGONAL.  The oracle correctly reports kappa(A, A) == 1.  The DAVID
# matrix carries 0 there (calculateKappaScores() never writes the diagonal)
# and the stock DistanceMetric matrix carries the -99 sentinel, so EVERY
# comparison below is over the UPPER TRIANGLE ONLY.  That is a property of
# those matrices, not of this oracle (T3-04a V2 owns the reasoning).
#
# MEASURED CONTEXT, recorded as comments and NOT asserted as constants.
# PLAN_V2_SPINE.md C2's remedy is that no fixture-derived constant is
# asserted; every assertion below is either a closed form derived in place
# from its own inputs, or structural.  On the 580-term anchor fixture:
#   upper-tri sum, raw (no clamp)      14239.673391   min -0.107613
#   upper-tri sum, clamped -- THE RULING   14597.851964   min 0, 0 negatives
#   pairs the clamp moves                  26159 of 167910
#   max |raw oracle - DAVID matrix|        0 exactly
#   max |oracle - stock DistanceMetric|    6.955802
#
# These tests read tests/testthat/fixtures/, which is .Rbuildignore'd, so they
# skip on an installed tarball and run from the source tree.
# ===========================================================================

# --- Shared, memoised runs -------------------------------------------------
#
# MOVED to helper-t3-02.R by SPEC-RC-004 M4, so that the v102-dependent blocks
# relocated to test-src-oracle-v102.R can reach them too.  Helpers are sourced
# before every test file, so rc_t3_02_*() resolves in both configurations.

# --- Contract preconditions ------------------------------------------------

test_that("the oracle's universe is the dataset's distinct gene count", {
  rc_skip_if_no_fixture()

  # N is "the number of distinct genes in the DATASET" -- not |A u B|, not the
  # token count, not the term count.  T3-01 V1 owns these two numbers; they are
  # re-derived here through the oracle's own universe function because getting
  # N wrong is the single most likely way for two oracles to disagree.
  genes <- rc_t3_02_genes()
  M <- rc_oracle_membership(genes)

  expect_identical(nrow(M), length(genes))
  expect_identical(ncol(M), length(rc_oracle_universe(genes)))
  expect_identical(ncol(M), length(unique(unlist(strsplit(genes, ",", fixed = TRUE)))))

  # Binary by construction, and |A| is a MEMBER count, not a token count.
  #
  # SPEC-RC-004 RE-DERIVATION (2026-08-26).  The third assertion below used to
  # read expect_gt(max tokens per term, ncol(M)) -- the fixture repeated gene
  # names inside single term strings, so the largest term carried MORE tokens
  # (measured then: > 1819) than the universe had genes.  SPEC-RC-006's
  # token-level deduplication inside merge_enrichment_results() removed those
  # repeats, so the property is now false BY DESIGN: measured on the re-frozen
  # post-dedup fixture, max tokens per term = 1639 <= ncol(M) = 1819.
  #
  # The old assertion survived only because the stored fx-input.rds was still
  # the PRE-dedup artifact; REQ-RC004-020's mandated re-freeze surfaced it.
  # This is the fence working, not a regression -- and it is NOT a silent
  # re-baseline: the member/token distinction is re-anchored below onto what is
  # still true, namely that no term can claim more members than the universe
  # holds and that tokens and members now coincide because dedup made them.
  expect_true(all(M %in% c(0, 1)))
  expect_true(max(rowSums(M)) <= ncol(M))
  expect_true(max(lengths(strsplit(genes, ",", fixed = TRUE))) <= ncol(M))
  expect_identical(unname(rowSums(M)),
                   as.numeric(lengths(strsplit(genes, ",", fixed = TRUE))))
})


# --- V1 -- the two teams' oracles agree with each other to 1e-12 -----------
#
# "before either team edits any package code.  A disagreement here means one
# oracle is wrong, and finding that out now costs nothing; finding it out after
# five files have moved costs the release."
#
# The exchanged object is written self-describing under its own filename
# rather than through rc_emit(): the T3-01(d) registry is a registry of
# PACKAGE ENTRY POINTS and errors on unregistered tags, and the oracle is not
# a package entry point.  See helper-oracle.R for the reasoning.

test_that("V1: this team's oracle matrix is emitted, self-describing, for exchange", {
  rc_skip_if_no_fixture()

  path <- rc_oracle_emit_crossteam(rc_t3_02_genes(), clamp = TRUE)
  expect_true(file.exists(path))

  rec <- readRDS(path)
  expect_identical(rec$kappa, rc_t3_02_oracle())          # round trips
  expect_identical(rec$clamp, TRUE)                       # the settled branch
  expect_identical(rec$n_terms, length(rc_t3_02_genes()))
  expect_identical(rec$n_universe, ncol(rc_oracle_membership(rc_t3_02_genes())))

  # Re-emitting identical content is a no-op; re-emitting DIFFERENT content
  # errors rather than silently destroying the evidence that the teams agreed.
  expect_silent(rc_oracle_emit_crossteam(rc_t3_02_genes(), clamp = TRUE))
})

test_that("V1: the two teams' oracles agree to 1e-12 on the fixture", {
  rc_skip_if_no_fixture()

  # SPEC-RC-004 REQ-RC004-021 -- DECLARED SKIP EXCEPTION (site S3).
  # Ruled a permanent, legitimate absence by SPEC-RC-004 spec.md section 5:
  # this needs the OTHER team's independently computed record, which has never
  # been supplied and does not exist in this repository.  Check-time generation
  # cannot manufacture a second team's independent computation, so this is a
  # cross-team exchange gate rather than a fixture-availability gate.
  counterpart <- rc_oracle_counterpart_file()
  skip_if_not(file.exists(counterpart), sprintf(paste0(
    "the other team's oracle has not been supplied.  V1 needs their record at ",
    "'%s' (or set RC_CROSSTEAM_ORACLE).  This team's half is emitted at '%s'."),
    counterpart, rc_oracle_crossteam_file()))

  theirs <- readRDS(counterpart)
  ours   <- rc_t3_02_oracle()

  # Same input, same OD-5 branch -- otherwise the comparison is meaningless and
  # the disagreement it would report is about the branch, not about either
  # oracle.  An unclamped and a clamped oracle differ on 26159 fixture pairs
  # by construction; that is the OD-5 trap and it must ERROR, not report.
  expect_identical(theirs$n_terms,    length(rc_t3_02_genes()))
  expect_identical(theirs$n_universe, ncol(rc_oracle_membership(rc_t3_02_genes())))
  if (!identical(theirs$clamp, TRUE)) {
    stop("OD-5 BRANCH MISMATCH: the supplied counterpart oracle was computed ",
         "with clamp = ", format(theirs$clamp), ".  OD-5 is settled on FLOOR ",
         "ALL NEGATIVES AT 0; an unclamped oracle disagrees with corrected ",
         "C++ on 26159 fixture pairs and fails T0-01 V1 on correct code.",
         call. = FALSE)
  }

  expect_identical(dim(theirs$kappa), dim(ours))
  expect_lt(max(abs(theirs$kappa - ours)), 1e-12)
})


# --- V2 -- the oracle agrees with the UNTOUCHED shipped COMPUTATION --------
#
# DavidClustering::calculateKappaScores() already counts the intersection
# correctly with a manual membership loop, so it is a valid kappa and predates
# every edit in this plan.  Agreement is required to 1e-12 on the T3-01 DAVID
# fixture (frozen 580-term vectors, 0.5 / 3 / 3 / 0.5).
#
# RUNNABILITY PRECONDITION -- T3-04a.  As shipped this V could not be run at
# all: kappaMatrix and calculateKappaScores() were private, run() returned only
# list(clusters = ...), and runDavidClustering is not in NAMESPACE.  T3-04a
# exposes the matrix additively and proves the exposure output-neutral.  The
# claim supported here is therefore "agrees with the untouched COMPUTATION,
# observed behind an additive accessor" -- NOT "agrees with a literally
# untouched tree".  IF THE AGREEMENT FAILS, THE ORACLE IS WRONG, NOT THE
# PACKAGE.
#
# T0-02's denominator cast: required by the plan to be present before this gate
# is judged.  At the fixture's N = 1819 the product is 3.3M, four orders of
# magnitude below the int32 ceiling, so the cast provably changes no number
# here (T0-02 V2 asserts exactly that by bit-identity) -- it makes the
# agreement a theorem instead of a coincidence.  Nothing below depends on
# whether it has landed yet.
#
# CONDITIONAL ON OD-5, RESTATED FOR THE RULING.  PLAN_V2_T3.md says this V
# "holds as written under OD-5 = remove; under retain the oracle and
# DavidClustering legitimately differ on disjoint pairs".  Under the SETTLED
# ruling -- floor ALL negatives -- they legitimately differ on ALL negative
# pairs, not merely the disjoint ones.  The V therefore splits in two, and
# both halves are asserted:
#   (a) the oracle's RAW core agrees with DavidClustering everywhere; and
#   (b) the CLAMPED shipping oracle differs from it on EXACTLY the negative
#       set and agrees elsewhere.
# (a) is the real independence check -- it is the half that cannot pass if the
# oracle has reproduced the intersection defect.

test_that("V2(a): the RAW oracle agrees with DavidClustering to 1e-12", {
  rc_skip_if_no_fixture()

  km  <- rc_t3_02_david_matrix()
  raw <- rc_t3_02_oracle_raw()
  expect_identical(dim(raw), dim(km))

  ut <- upper.tri(km)          # diagonal excluded: DAVID's is 0, the oracle's 1
  expect_lt(max(abs(raw[ut] - km[ut])), 1e-12)
})

test_that("V2(b): the CLAMPED oracle differs from DavidClustering on exactly the negatives", {
  rc_skip_if_no_fixture()

  km  <- rc_t3_02_david_matrix()
  ok  <- rc_t3_02_oracle()
  ut  <- upper.tri(km)

  negative <- km[ut] < 0
  # There ARE negatives; otherwise the OD-5 clamp is untestable on this fixture
  # and this test proves nothing.  Structural, not the measured 26159.
  expect_gt(sum(negative), 0L)

  # Off the clamped set the two are the same computation.
  expect_lt(max(abs(ok[ut][!negative] - km[ut][!negative])), 1e-12)
  # On it the oracle reports exactly 0 -- the settled OD-5 contract.
  expect_true(all(ok[ut][negative] == 0))
  # And the difference set is EXACTLY the negative set: no third population.
  expect_identical(which(abs(ok[ut] - km[ut]) > 1e-12), which(negative))
})

test_that("V2: the clamped oracle is on the settled [0, 1] scale", {
  rc_skip_if_no_fixture()

  ok  <- rc_t3_02_oracle()
  raw <- rc_t3_02_oracle_raw()

  expect_identical(sum(ok < 0), 0L)               # the ruling: no negatives
  expect_true(all(ok >= 0 & ok <= 1))
  expect_gt(sum(raw < 0), 0L)                     # the clamp did real work
  # The clamp moves the total UP, and only up.
  expect_true(all(ok >= raw))
  expect_gt(sum(ok[upper.tri(ok)]), sum(raw[upper.tri(raw)]))
})


# --- V3 -- closed forms ----------------------------------------------------
#
# "The oracle reproduces the hand cases of T0-01 V2 to < 1e-15."
#
# Every constant below is DERIVED IN PLACE from its own inputs (PLAN_V2_SPINE
# C2's remedy), so the assertion cannot be satisfied by a transcription error
# that happens to match.  Asserted with a tolerance, never with identical() or
# == : the exact double is one ulp from the naively written decimal, and
# revision 1's "assert exact equality" failed on correct code (critique m2a).

test_that("V3: the T0-01 V2 hand case, N=100 |A|=20 |B|=30 |AnB|=10", {
  N <- 100; a <- 20; b <- 30; common <- 10

  # Derived in place, from the definition, with no reference to the oracle:
  #   both = 10, neither = 100 - 20 - 30 + 10 = 60  ->  oab = 70/100
  #   aab  = (20*30 + 80*70) / 100^2 = 6200/10000
  #   kappa = (0.70 - 0.62) / (1 - 0.62) = 0.08/0.38 = 4/19
  expected <- ((common + (N - a - b + common)) / N -
                 (a * b + (N - a) * (N - b)) / N^2) /
              (1 - (a * b + (N - a) * (N - b)) / N^2)
  expect_lt(abs(expected - 4 / 19), 1e-15)              # the derivation checks out

  # Positive, so the OD-5 clamp is a no-op here and both branches agree.
  expect_lt(abs(rc_oracle_kappa(a, b, common, N, clamp = FALSE) - expected), 1e-15)
  expect_lt(abs(rc_oracle_kappa(a, b, common, N, clamp = TRUE)  - expected), 1e-15)

  # The plan quotes 0.21052631578947381 for this case.  That decimal is the
  # compiled build's double; it and 4/19 differ by less than an ulp of the
  # tolerance, which is why the plan mandates < 1e-15 rather than identical().
  expect_lt(abs(rc_oracle_kappa(a, b, common, N, clamp = TRUE) -
                  0.21052631578947381), 1e-15)
})

test_that("V3: self-identity, disjointness and the degenerate universe", {

  # kappa(A, A) == 1, exactly.  oab = 1, so the quotient is (1-aab)/(1-aab).
  expect_identical(rc_oracle_kappa(37, 37, 37, 500, clamp = FALSE), 1)
  expect_identical(rc_oracle_kappa(37, 37, 37, 500, clamp = TRUE),  1)

  # Disjoint sets, |A| + |B| < N.  Derived in place:
  #   N = 100, |A| = 20, |B| = 30, |AnB| = 0
  #   both = 0, neither = 50   -> oab = 0.50
  #   aab  = (600 + 5600)/10000 = 0.62
  #   raw kappa = (0.50 - 0.62)/0.38 = -0.12/0.38 = -6/19  < 0
  N <- 100; a <- 20; b <- 30
  raw_disjoint <- ((0 + (N - a - b)) / N - (a * b + (N - a) * (N - b)) / N^2) /
                  (1 - (a * b + (N - a) * (N - b)) / N^2)
  expect_lt(abs(raw_disjoint - (-6 / 19)), 1e-15)
  expect_lt(raw_disjoint, 0)

  # The raw core reproduces it; the SETTLED OD-5 clamp reports 0 instead.
  # PLAN_V2_T3.md's contract has no clamp and PLAN_V2_T0.md T0-01 V3's
  # "k < 0 for disjoint sets" clause encodes the removal branch; both are
  # superseded by the ruling.  Asserting either of them here would fail on
  # correct code.
  expect_lt(abs(rc_oracle_kappa(a, b, 0, N, clamp = FALSE) - raw_disjoint), 1e-15)
  expect_identical(rc_oracle_kappa(a, b, 0, N, clamp = TRUE), 0)

  # The degenerate aab == 1 branch: both terms span the whole universe, which
  # is the MAXIMALLY SIMILAR case, so kappa is 1 and not 0.  (T0-01 change (3)
  # makes DistanceMetric agree; DavidClustering already returns 1.0.)
  expect_identical(rc_oracle_kappa(500, 500, 500, 500, clamp = FALSE), 1)
  expect_identical(rc_oracle_kappa(0, 0, 0, 500, clamp = FALSE), 1)
})

test_that("V3: the oracle refuses impossible count tuples", {
  # An oracle that silently accepts |AnB| > min(|A|,|B|), or |A| > N, would
  # return a finite number for an input that has no 2x2 table, and a caller
  # could then compare it to something.
  expect_error(rc_oracle_kappa(20, 30, 25, 100), "cannot exceed min")
  expect_error(rc_oracle_kappa(120, 30, 10, 100), "cannot exceed the gene universe")
  expect_error(rc_oracle_kappa(60, 60, 0, 100), "|A u B| cannot exceed", fixed = TRUE)
  expect_error(rc_oracle_kappa(1, 1, 1, 0), "must be positive")
})


# --- V4 -- invariants over 200000 random draws -----------------------------
#
# Properties of the ORACLE's own arithmetic, independent of the package.
# Seeded, so a failure is reproducible.

test_that("V4: range, symmetry and self-similarity over 200000 random draws", {
  set.seed(20260810)
  n <- 200000L

  N <- sample(2:5000, n, replace = TRUE)
  a <- as.integer(round(runif(n) * N))
  b <- as.integer(round(runif(n) * N))
  lo <- pmax(0, a + b - N)                  # |AnB| >= |A|+|B|-N
  hi <- pmin(a, b)                          # |AnB| <= min(|A|,|B|)
  common <- lo + as.integer(round(runif(n) * (hi - lo)))

  raw <- rc_oracle_kappa(a, b, common, N, clamp = FALSE)
  cl  <- rc_oracle_kappa(a, b, common, N, clamp = TRUE)

  # Range.  The raw scale is [-1, 1]; the SETTLED exported scale is [0, 1].
  expect_true(all(raw >= -1 & raw <= 1))
  expect_true(all(cl >= 0 & cl <= 1))
  expect_true(all(is.finite(raw)))
  # Both bounds are actually reached, so the range assertion is not vacuous.
  expect_lt(min(raw), 0)
  expect_gt(max(raw), 0.99)

  # Symmetry: kappa(A, B) == kappa(B, A), exactly.  Every term of both oab and
  # aab is symmetric under the swap and IEEE multiplication is commutative, so
  # this is bit-equality and not a tolerance.
  expect_identical(raw, rc_oracle_kappa(b, a, common, N, clamp = FALSE))
  expect_identical(cl,  rc_oracle_kappa(b, a, common, N, clamp = TRUE))

  # Self-similarity: kappa(A, A) == 1 for every draw, including the degenerate
  # |A| = N case that T0-01 V3 calls out explicitly.
  expect_true(all(rc_oracle_kappa(a, a, a, N, clamp = FALSE) == 1))
  expect_true(all(rc_oracle_kappa(N, N, N, N, clamp = FALSE) == 1))

  # The clamp is exactly max(raw, 0) and nothing else.
  expect_identical(cl, pmax(raw, 0))
})

test_that("V4: kappa -> Dice as N -> inf", {
  set.seed(20260810)
  n <- 200000L
  N <- 1e6                                   # T0-01 V3's stated regime

  a <- sample(1:100, n, replace = TRUE)      # |A|, |B| << N is the limit's
  b <- sample(1:100, n, replace = TRUE)      # premise, not a convenience
  common <- as.integer(round(runif(n) * pmin(a, b)))

  k <- rc_oracle_kappa(a, b, common, N, clamp = FALSE)
  d <- rc_oracle_dice(a, b, common)
  expect_lt(max(abs(k - d)), 1e-4)

  # The limit is a limit: the gap shrinks like 1/N.  Asserting convergence
  # rather than a single threshold is what makes this a property and not a
  # tuned constant.
  gaps <- vapply(c(1e3, 1e4, 1e5, 1e6, 1e7), function(NN)
    abs(rc_oracle_kappa(40, 60, 20, NN, clamp = FALSE) - rc_oracle_dice(40, 60, 20)),
    numeric(1))
  expect_true(all(diff(gaps) < 0))
  expect_lt(gaps[5], gaps[1] / 1000)
})



