# tests/testthat/test-news-disclosure.R
# ===========================================================================
# SPEC-RC-005 -- reviewer finding A6.
#
# The NEWS.md disclosure obligations were asserted OUTSIDE the package, so
# nothing enforced them.  They are enforced here instead.
#
# These assertions match on CONTENT, never on line numbers, byte offsets or
# bullet counts, so SPEC-RC-001's insertions into NEWS.md cannot break them.
#
# There is deliberately NO conditional-absence path anywhere in this file.  An
# unreachable NEWS.md is precisely the condition the file exists to detect, so
# it must be reported as a failure (REQ-RC-077, same reasoning as finding B1);
# passing over it silently would reproduce the defect being closed.
#
# Wording re-ruled by the author on 2026-09-11: NEWS.md states what 2.0.0 now
# does, in a neutral register, rather than that earlier results were wrong.
# The same seven disclosures are still required; only the phrases matched
# changed.  Whitespace is matched with \s+ so re-wrapping NEWS.md cannot
# break a check.
#
# Softened further on 2026-09-15 per collaborator review: "results will
# differ" overstated the claim (some analyses may be unaffected depending on
# options and data), so NEWS.md now says "results may differ".  The pattern
# matches either wording so a future re-strengthening is not itself a
# regression.
# ===========================================================================

test_that("A6: NEWS.md is installed and reachable from the suite", {
  # A failure here means the disclosure cannot be checked at all.
  path <- system.file("NEWS.md", package = "richCluster")
  expect_true(nzchar(path),
              info = "NEWS.md is not installed; the disclosure cannot be verified")
  expect_true(file.exists(path))
})

test_that("A6: NEWS.md carries every required 1.0.2 disclosure", {
  path <- system.file("NEWS.md", package = "richCluster")
  expect_true(nzchar(path))
  news <- paste(readLines(path, warn = FALSE), collapse = "\n")

  required <- list(
    "results differ from 1.0.2 and earlier" =
      "results\\s+(?:will|may)\\s+differ\\s+from\\s+those\\s+of\\s+richCluster\\s+1\\.0\\.2",
    "re-running analyses is recommended" =
      "recommend\\s+re-running",
    "single and complete were exchanged" =
      "(?s)\"single\".{0,120}\"complete\".{0,160}exchanged",
    "jaccard now returns the Jaccard index" =
      "now\\s+returns\\s+the\\s+Jaccard\\s+index",
    "kappa is on a [0, 1] scale" =
      "\\[0,\\s+1\\]\\s+scale",
    "saved cluster_result objects should be recomputed" =
      "(?s)cluster_result.{0,80}should\\s+be\\s+recomputed",
    "there is no compatibility mode" =
      "no\\s+compatibility\\s+mode"
  )

  for (nm in names(required)) {
    expect_true(
      grepl(required[[nm]], news, ignore.case = TRUE, perl = TRUE),
      info = paste0("NEWS.md does not disclose: ", nm,
                    "\n  (searched for /", required[[nm]], "/)"))
  }
})
