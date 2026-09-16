# tests/testthat/test-cpp-entry-guards.R
# D1 / D2 / D3 -- guards at the runRichCluster() entry point.
#
# runRichCluster() is in NAMESPACE, so every check cluster() performs in
# validate_inputs() is skippable by calling it directly; D1 and D3 are both
# reachable from cluster() as well.  All inputs are inline, so this file has no
# fixture dependency and ships in the tarball.  Do NOT rename to test-src-*:
# that prefix is tarball-excluded.

# Two terms sharing no genes: kappa(A, B) = 0 (measured).  Any cutoff in the
# documented domain leaves them apart, so a merge here means the cutoff was
# never checked.
rc_guard_terms <- c("A", "B")
rc_guard_genes <- c("g1,g2", "g3,g4")

rc_guard_run <- function(...) {
  sink(nullfile())
  on.exit(sink(), add = TRUE)
  runRichCluster(rc_guard_terms, rc_guard_genes, ...)
}

test_that("D3: out-of-domain cutoffs error instead of returning a wrong answer", {
  msg_d <- "distance_cutoff must be between 0 and 1."
  msg_l <- "linkage_cutoff must be between 0 and 1."

  # Measured before the guard: dc = lc = -1 returned the single cluster
  # "B, A" -- kappa 0 > -1 admits every pair, so two terms sharing no genes
  # were merged.
  expect_error(rc_guard_run("kappa", -1, "average", -1), msg_d, fixed = TRUE)

  # Measured before the guard: each of these returned 2 singleton clusters,
  # because every comparison against NaN is false.  A naive range check
  # (x <= 0 || x > 1) still lets them through; !(x > 0) is what rejects them.
  expect_error(rc_guard_run("kappa", NA_real_, "average", 0.5), msg_d, fixed = TRUE)
  expect_error(rc_guard_run("kappa", NaN,      "average", 0.5), msg_d, fixed = TRUE)
  expect_error(rc_guard_run("kappa", Inf,      "average", 0.5), msg_d, fixed = TRUE)
  expect_error(rc_guard_run("kappa", 0.5, "average", NA_real_), msg_l, fixed = TRUE)
  expect_error(rc_guard_run("kappa", 0.5, "average", NaN),      msg_l, fixed = TRUE)
  expect_error(rc_guard_run("kappa", 0.5, "average", Inf),      msg_l, fixed = TRUE)

  # The endpoints of the documented domain, 0 < cutoff <= 1.
  expect_error(rc_guard_run("kappa", 0,   "average", 0.5), msg_d, fixed = TRUE)
  expect_error(rc_guard_run("kappa", 1.5, "average", 0.5), msg_d, fixed = TRUE)
  expect_error(rc_guard_run("kappa", 0.5, "average", 0),   msg_l, fixed = TRUE)
  expect_error(rc_guard_run("kappa", 0.5, "average", 1.5), msg_l, fixed = TRUE)
})

test_that("D3: the guard admits the whole documented domain", {
  # cutoff = 1 is legal (the rule is 0 < cutoff <= 1, not < 1), and the guard
  # must not narrow what cluster()'s validate_inputs() accepts.
  for (cut in c(1e-9, 0.5, 1)) {
    res <- rc_guard_run("kappa", cut, "average", cut)
    expect_s3_class(res$all_clusters, "data.frame")
  }
  # The disjoint pair still does not merge at the default cutoffs.
  res <- rc_guard_run("kappa", 0.5, "average", 0.5)
  expect_equal(nrow(res$all_clusters), 2)
})

test_that("D1: an empty gene_delim errors instead of hanging forever", {
  # Before the guard this did not return: both tokenisers advance with
  # start = end + delimiter.length(), and find("", start) returns start with
  # length() 0.  Measured against the shipped function itself --
  # StringUtils::splitStringToUnorderedSet("g1,g2", "") had not returned after
  # 8 s -- and through runRichCluster(), where the R session had to be killed
  # (SIGKILL after 25 s, exit 137).  The loop allocates nothing, so it spins
  # silently: no memory growth, no output, no way out but kill -9.
  expect_error(rc_guard_run("kappa", 0.5, "average", 0.5, ""),
               "gene_delim must be a non-empty string", fixed = TRUE)
})

test_that("D1: the documented top-level API no longer hangs on gene_delim = ''", {
  d1 <- data.frame(Term = c("T1", "T2"), GeneID = c("a,b,c", "a,b,d"),
                   Pvalue = 1e-5, Padj = 1e-4, stringsAsFactors = FALSE)
  d2 <- data.frame(Term = c("T1", "T2"), GeneID = c("a,b,c", "a,b,e"),
                   Pvalue = 1e-5, Padj = 1e-4, stringsAsFactors = FALSE)
  # Only that it errors, and that the message names the offending parameter --
  # whether the rejection lands in R or in C++ is not this file's business.
  sink(nullfile())
  on.exit({ while (sink.number() > 0) sink() }, add = TRUE)
  expect_error(cluster(list(d1, d2), gene_delim = ""), "gene_delim")
})

test_that("D2: a long C++ run answers R's interrupt checks", {
  # R checks elapsed-time limits at the same points it checks for Ctrl-C, so a
  # limit that fires mid-run is a direct test of Rcpp::checkUserInterrupt().
  # Measured before the checks were added: this input ran 12.9 s to completion
  # under a 1 s limit, and four real SIGINTs over 8 s did not stop a 500-term
  # run either.  Deliberately not skip_on_cran()'d: the test costs ~2 s once the
  # checks are in (it stops at the 1 s limit), and gating it would leave the
  # only regression guard for interruptibility switched off by default.  The
  # timing margin is wide -- 13.6 s unguarded against a 5 s ceiling.
  n <- 1600
  set.seed(42)
  genes <- vapply(seq_len(n),
                  function(i) paste(sample(paste0("g", 1:200), 40), collapse = ","),
                  character(1))
  terms <- paste0("T", seq_len(n))

  on.exit({
    setTimeLimit(cpu = Inf, elapsed = Inf)
    while (sink.number() > 0) sink()
  }, add = TRUE)

  started <- Sys.time()
  outcome <- tryCatch({
    setTimeLimit(elapsed = 1, transient = TRUE)
    sink(nullfile())
    runRichCluster(terms, genes, "kappa", 0.1, "average", 0.1)
    while (sink.number() > 0) sink()
    "completed"
  },
  interrupt = function(e) "stopped",
  error = function(e) paste0("error: ", conditionMessage(e)))
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  setTimeLimit(cpu = Inf, elapsed = Inf)
  while (sink.number() > 0) sink()

  # The run must not have finished, and it must have given up promptly -- the
  # unguarded run of this same input takes ~13 s, so a stop inside 5 s can only
  # come from an interrupt check inside the C++ loops.
  expect_false(identical(outcome, "completed"))
  expect_lt(elapsed, 5)
  # runRichCluster()'s blanket catch (...) must not have swallowed the
  # interrupt: that would report it as an unknown C++ exception and leave the
  # run just as unstoppable.
  expect_false(grepl("Unknown C++ exception", outcome, fixed = TRUE))
})
