# tests/testthat/test-david-kappa-matrix.R
# ===========================================================================
# T3-04a (WAVE 0) -- ADDITIVE, OUTPUT-NEUTRAL EXPOSURE OF THE DAVID KAPPA
# MATRIX.  This is the FIRST package-code edit in the project, and its entire
# acceptance argument is bit-identity against the snapshot T3-03 froze.
#
# THE GAP IT CLOSES.  DavidClustering::calculateKappaScores() already counts
# the gene intersection correctly with a manual membership loop, so it is the
# in-repo reference implementation the plan's single most load-bearing
# verification (the 1e-12 oracle-agreement gate: T3-02 V2, T0-01 V4) is judged
# against.  As shipped, NONE of those gates was executable: `kappaMatrix` and
# `calculateKappaScores()` were both private members of the C++ class, `run()`
# returned only `list(clusters = ...)`, and `runDavidClustering` is not in
# NAMESPACE -- the matrix was not observable from R by any means.
#
# WHAT MOVED: only observability.  `run()`, `calculateKappaScores()`,
# `runDavidClustering()` and `david_cluster()` are untouched; a public const
# accessor on the C++ class and a second, additional entry point
# (`runDavidClusteringWithKappa()`) were added alongside them.  The claim this
# licenses is "agrees with the untouched COMPUTATION, observed behind an
# additive accessor" -- NOT "agrees with a literally untouched tree".
#
# V1 IS TWO ASSERTIONS AGAINST TWO DISTINCT TAGS, ONE PER DAVID ENTRY POINT,
# and they are not collapsible: runDavidClustering() takes the frozen 580-term
# Term/GeneID VECTORS at 0.5/3/3/0.5, while david_cluster()'s first formal
# `enrichment_results` is a LIST OF DATA FRAMES and it applies no p-value
# filter (3403 terms).  They are different objects at different entry points;
# citing one while comparing against the other's tag is exactly the mismatch
# T3-01 V6 requires rc_expect_identical() to ERROR on, and the third test
# below asserts that error rather than trusting the prose.  In particular
# "david_cluster() on the DAVID fixture" is NOT A RUNNABLE CALL and must never
# be written.
#
# WHAT IS NOT ASSERTED HERE.  V4's second half -- "agrees with the pure-R
# oracle to 1e-12" -- needs the oracle, and the oracle is T3-02, which
# DEPENDS ON THIS ITEM (the edge runs T3-04a -> T3-02, one way only).  It is
# therefore not runnable at this item's landing by construction of the
# dependency graph, and its shipped home is T3-02 V2, which this item is what
# makes executable at all.  This file asserts V4's runnable half: the exposed
# matrix must DISAGREE with DistanceMetric::getKappa -- a build in which the
# two agree has reproduced the std::set_intersection defect in both places and
# every downstream gate is worthless.
#
# These tests read tests/testthat/fixtures/, which is .Rbuildignore'd, so they
# skip on an installed tarball and run from the source tree.
#
# 2026-08-17 (SPEC-RC-009): T3-04a's V1 neutrality argument was discharged at
# its landing on the stock computation.  DS-03's fix deliberately moves the
# merge stage, so V1/V3 are re-anchored to the v110-* family (the corrected
# DAVID path); the deliberate divergence from v102-david is asserted
# explicitly below.
# ===========================================================================

# --- Shared, memoised runs -------------------------------------------------
#
# MOVED to helper-t3-04a.R by SPEC-RC-004 M4, so that the v102-dependent blocks
# relocated to test-src-david-v102.R can reach them too.  Helpers are sourced
# before every test file, so rc_t3_04a_*() resolves in both configurations.

# --- V1 -- OUTPUT NEUTRALITY (primary) -------------------------------------
#
# The item is breaking_change: NO, and V1 PROVES that rather than asserting
# it.  An implementation that fails (a) has moved the C++ return shape; one
# that fails (b) has moved the exported wrapper's output.  Either way it is
# not additive.

test_that("V1(a): runDavidClustering() on the DAVID fixture is identical to v110-david", {
  rc_skip_if_no_fixture()
  rc_skip_if_no_artifact("v110-david")

  out <- rc_t3_04a_david()
  rc_expect_identical("v110-david", out, "runDavidClustering()")
  # The names() clause: no existing element removed, renamed or reordered.
  expect_identical(names(out), names(rc_artifact_object("v110-david")))
})

# (b) re-pointed from v102-dc to w1-dc by SPEC-RC-006 (author ruling 2026-08-24).
# v102-dc stores the WHOLE result object, merged_df included, and SPEC-RC-006
# deliberately reshapes that frame: +DatasetCount, and gene tokens deduplicated
# on 1723 of 3403 rows.  So (b) began failing the moment RC-006 landed.  It was
# measured component by component first, and the clustering itself is unmoved --
# clusters and final_clusters are both identical() across the change; only
# merged_df and cluster_df differ, exactly where RC-006 intends.  SPEC-RC-006
# section F.3 predicted this and its M8 emitted w1-dc as the post-RC-006 anchor,
# but no milestone re-pointed this assertion at it.  That is what this edit does.
# The neutrality contract is unchanged in kind: (b) still proves the exported
# wrapper's output has not moved, now against the current-era anchor.
# SPEC-RC-009 re-anchor: (b) was re-pointed v102-dc -> w1-dc by SPEC-RC-006 and
# is re-pointed w1-dc -> v110-dc here.  DS-03 deliberately moves the merge
# stage, so the david_cluster() return moves with it and w1-dc ceases to be the
# current-era anchor.  The neutrality contract is unchanged in kind: (b) still
# proves the exported wrapper's output has not moved, now against the anchor
# frozen immediately after this SPEC's fix.
test_that("V1(b): david_cluster(list(d1, d2)) at defaults is identical to v110-dc", {
  skip_on_cran()   # full two-dataset DAVID run, ~28 s locally; CRAN caps a check at 10 min
  rc_skip_if_no_fixture()
  rc_skip_if_no_artifact("v110-dc")

  dc <- rc_t3_04a_dc()
  rc_expect_identical("v110-dc", dc, "david_cluster()")
  expect_identical(names(dc), names(rc_artifact_object("v110-dc")))
})




# --- V2 -- the matrix is reachable and well formed -------------------------
#
# NULL-ACTION: this V fails by NON-EXISTENCE on the stock build.  Measured
# there: zero namespace objects matching "kappa", two registered .Call entry
# points (runDavidClustering, runRichCluster), and names(runDavidClustering())
# == "clusters".  There was no R-visible path to the matrix at all, which is
# precisely the defect.

test_that("V2: the exposure is reachable from R", {
  expect_true(is.function(richCluster:::runDavidClusteringWithKappa))
  expect_true(is.function(richCluster:::david_kappa_matrix))
  expect_true("_richCluster_runDavidClusteringWithKappa" %in%
                names(getDLLRegisteredRoutines("richCluster")$.Call))
})

test_that("V2: the exposed matrix is a well-formed n_terms x n_terms kappa matrix", {
  rc_skip_if_no_fixture()

  dfx <- rc_david_fixture()
  res <- rc_t3_04a_david_with_kappa()
  km  <- res$kappa_matrix

  # Additive on the C++ return: `clusters` first, `kappa_matrix` appended.
  expect_identical(names(res), c("clusters", "kappa_matrix"))

  n <- length(dfx$Term)
  expect_true(is.matrix(km))
  expect_identical(typeof(km), "double")
  expect_identical(dim(km), c(n, n))

  # Symmetric (calculateKappaScores() writes [i][j] and [j][i] the same double,
  # so this is exact, not approximate) and finite everywhere.
  expect_true(identical(km, t(km)))
  expect_true(all(is.finite(km)))

  # THE DIAGONAL IS 0, NOT 1.  calculateKappaScores() loops
  # `for (j = i + 1; ...)` and never writes it, so it keeps the constructor's
  # 0.0 fill (DavidClustering.cpp: kappaMatrix.resize(n_terms,
  # std::vector<double>(n_terms, 0.0))).  EXCLUDE THE DIAGONAL from every
  # oracle comparison built on this matrix (T3-02 V2, T0-01 V4): a correct
  # pure-R Cohen's-kappa oracle gives kappa(A, A) == 1 there, so an inclusive
  # comparison would fail on n_terms entries for a reason that has nothing to
  # do with the intersection defect.  Do NOT "fix" the diagonal here -- this
  # item is additive and output-neutral, and the diagonal is never read by the
  # clustering.
  expect_identical(unique(diag(km)), 0)
})

test_that("V2: the thin R accessor returns the same matrix as the C++ entry point", {
  # A tiny synthetic input: this asserts the wrapper's plumbing, not a number,
  # so it needs neither the fixture nor a 7-second run.  Both sides are
  # computed inside this test run, so it is exempt from the tag rule.
  terms <- c("t1", "t2", "t3", "t4")
  genes <- c("a,b,c,d", "a,b,c,e", "a,b,f,g", "x,y,z")
  direct <- rc_t3_04a_quiet(
    richCluster:::runDavidClusteringWithKappa(terms, genes, 0.5, 3, 3, 0.5))
  viaR <- rc_t3_04a_quiet(
    richCluster:::david_kappa_matrix(terms, genes))
  expect_identical(viaR, direct$kappa_matrix)
  expect_identical(dim(viaR), c(4L, 4L))
})


# --- V3 -- it is the same matrix the clustering used -----------------------
#
# Both sides of every comparison below are computed inside this test run, so
# they are unambiguous by construction and exempt from the tag rule.

test_that("V3: the exposed matrix is the one the clustering consumed", {
  rc_skip_if_no_fixture()
  rc_skip_if_no_artifact("v110-david")

  dfx  <- rc_david_fixture()
  base <- rc_t3_04a_david_with_kappa()

  # Same pipeline: the clustering returned ALONGSIDE the exposed matrix is
  # bit-identical to what the untouched entry point returns, and to the
  # clustering frozen under v110-david (re-anchored by SPEC-RC-009: DS-03
  # deliberately moves the merge stage, so v102-david is no longer current).
  expect_identical(base$clusters, rc_t3_04a_david()$clusters)
  expect_identical(base$clusters, rc_artifact_object("v110-david")$clusters)

  # Perturb ONE gene string and re-run.  The perturbed term is the first term
  # the frozen clustering places in a cluster, so a change to it must be
  # visible in BOTH outputs.  If the accessor were reading a stale or
  # separately computed copy, one of the two would not move.
  frozen <- rc_artifact_object("v110-david")$clusters
  idx0   <- as.integer(strsplit(frozen$TermIndices[1], ", ", fixed = TRUE)[[1]])[1]
  i      <- idx0 + 1L                       # C++ term indices are 0-based
  genes  <- dfx$GeneID
  genes[i] <- "RC_T3_04A_PERTURBATION_SENTINEL"

  pert <- rc_t3_04a_quiet(richCluster:::runDavidClusteringWithKappa(
    dfx$Term, genes,
    dfx$similarity_threshold, dfx$initial_group_membership,
    dfx$final_group_membership, dfx$multiple_linkage_threshold))

  expect_false(identical(pert$kappa_matrix, base$kappa_matrix))
  expect_false(identical(pert$clusters, base$clusters))
  # and specifically the perturbed term's own row moved
  expect_false(identical(pert$kappa_matrix[i, ], base$kappa_matrix[i, ]))
})


