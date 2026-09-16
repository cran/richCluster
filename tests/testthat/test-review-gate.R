# tests/testthat/test-review-gate.R
# ===========================================================================
# SPEC-RC-005 -- THE WAVE REJECTION GATE.
#
# THE PROBLEM THIS EXISTS FOR
# ---------------------------
# Two reviewers returned NEEDS REVISION on Wave 0.  No revision pass ran, and
# later stages proceeded on top of unrevised work, because the wave workflow had
# no mechanism by which a rejection could halt anything.
#
# THE PROCEDURE
# -------------
#   1. A reviewer returning NEEDS REVISION adds ONE ROW PER FINDING to
#      review-ledger.csv, with status = OPEN and resolved_by empty.
#   2. NO LATER WAVE MAY BE STARTED while any row reads status = OPEN.
#   3. Closing a row requires setting status = CLOSED and naming, in
#      resolved_by, the SPEC that closed it.
#
# WHAT THIS TEST ENFORCES
# -----------------------
# Step 2, mechanically: while any row is OPEN the suite is RED, and a red suite
# halts every stage gated on a green suite.  That is real enforcement.
#
# WHAT THIS TEST CANNOT ENFORCE  (stated plainly, not papered over)
# -----------------------------------------------------------------
# It CANNOT DETECT A REVIEW THAT WAS NEVER RECORDED.  If a reviewer returns
# NEEDS REVISION and nobody adds the rows, this ledger is silent and this test
# is green.  Recording the verdict is a HUMAN OBLIGATION that no test can
# discharge.  Do not read a green gate as "the work was reviewed"; read it only
# as "no recorded finding is open".
#
# The ledger lives under tests/testthat/ -- NOT under tools/, which
# .Rbuildignore excludes -- so the gate ships in the tarball and runs in the
# same configuration CRAN runs.  A gate absent from the shipped build would
# repeat reviewer finding B1 inside its own fix.
# ===========================================================================

rc_review_ledger_path <- function() testthat::test_path("review-ledger.csv")

rc_review_ledger <- function() {
  utils::read.csv(rc_review_ledger_path(), stringsAsFactors = FALSE,
                  colClasses = "character")
}

test_that("the review ledger is present and well formed", {
  # A gate that vanishes when its input is missing is not a gate: FAIL, never
  # skip.  (Same reasoning as reviewer finding B1.)
  expect_true(file.exists(rc_review_ledger_path()),
              info = paste0("review-ledger.csv is missing at ",
                            rc_review_ledger_path(),
                            " -- the wave rejection gate cannot run."))

  led <- rc_review_ledger()
  expect_identical(names(led),
                   c("wave", "finding_id", "reviewer", "severity",
                     "verdict", "status", "resolved_by", "summary"))
  expect_gt(nrow(led), 0)
  expect_true(all(led$status  %in% c("OPEN", "CLOSED")),
              info = paste("bad status values:",
                           paste(setdiff(led$status, c("OPEN", "CLOSED")),
                                 collapse = ", ")))
  expect_true(all(led$verdict %in% c("ACCEPTED", "NEEDS REVISION")),
              info = paste("bad verdict values:",
                           paste(setdiff(led$verdict,
                                         c("ACCEPTED", "NEEDS REVISION")),
                                 collapse = ", ")))
  # A closed finding must name the SPEC that closed it.
  closed <- led[led$status == "CLOSED", , drop = FALSE]
  expect_true(all(nzchar(closed$resolved_by)),
              info = paste("CLOSED with no resolved_by:",
                           paste(closed$finding_id[!nzchar(closed$resolved_by)],
                                 collapse = ", ")))
})

test_that("THE GATE: no review finding is left open", {
  expect_true(file.exists(rc_review_ledger_path()))
  led  <- rc_review_ledger()
  open <- led[led$status == "OPEN", , drop = FALSE]

  expect_identical(
    nrow(open), 0L,
    info = paste0(
      "WAVE HALTED -- ", nrow(open), " review finding(s) still OPEN.\n",
      "No later wave may be started until every one is CLOSED.\n\n",
      paste(sprintf("  wave %s  %s  [%s]  %s",
                    open$wave, open$finding_id, open$severity, open$summary),
            collapse = "\n")))
})

test_that("every Wave-0 reviewer finding is recorded", {
  led <- rc_review_ledger()
  w0  <- led[led$wave == "0", , drop = FALSE]
  expect_setequal(w0$finding_id,
                  c("A1", "A2", "A3", "A4", "A5", "A6", "A7",
                    "B1", "B2", "B3", "B4"))
})
