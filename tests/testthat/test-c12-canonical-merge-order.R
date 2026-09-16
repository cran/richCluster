# C12 -- the merged row order must not depend on the session's collation.
#
# base::merge(by = "Term") sorts its result under LC_COLLATE, and the
# single-dataset branch used a bare order(), so merge_enrichment_results()
# returned a different row order in different locales.  cluster() survived that
# because C9 made everything downstream of the seed step canonical, but
# david_cluster() walks terms greedily in the order it is given and returned a
# different clustering per locale on the package's own example data -- 532, 533
# or 535 clusters depending only on LC_COLLATE (measured 2026-09-05 across
# Linux, macOS and Windows; see the converge ledger).
#
# The contract these tests hold is byte order (radix), which is defined without
# reference to any locale.  Note that is.unsorted() is itself collation-aware
# and therefore cannot express this contract.

c12_frame <- function(terms) {
  n <- seq_along(terms)
  data.frame(
    Term       = terms,
    Annot      = "GO",
    Annotated  = 10L,
    Significant = 3L,
    Pvalue     = n / 1000,
    Padj       = n / 500,
    GeneID     = paste0("g", n),
    stringsAsFactors = FALSE
  )
}

# Case and the hyphen-versus-space position are precisely what a collation
# reorders: under C.UTF-8 these sort case-insensitively with punctuation
# folded, under C they sort by byte.
c12_terms <- c("ADP metabolic process", "actin binding",
               "actin filament-based movement", "actin-mediated cell contraction",
               "Actin polymerization", "actin cytoskeleton")

test_that("C12: a single dataset merges to canonical byte order", {
  m <- merge_enrichment_results(list(c12_frame(c12_terms)))
  expect_identical(m$Term, sort(c12_terms, method = "radix"))
})

test_that("C12: two datasets merge to canonical byte order", {
  m <- merge_enrichment_results(list(c12_frame(c12_terms[1:4]),
                                     c12_frame(c12_terms[3:6])))
  expect_identical(m$Term, sort(unique(c12_terms), method = "radix"))
})

test_that("C12: the merged order is identical under every available collation", {
  old <- Sys.getlocale("LC_COLLATE")
  on.exit(suppressWarnings(Sys.setlocale("LC_COLLATE", old)), add = TRUE)

  seen <- list()
  for (loc in unique(c("C", "C.UTF-8", "en_US.UTF-8", old))) {
    if (!nzchar(suppressWarnings(Sys.setlocale("LC_COLLATE", loc)))) next
    seen[[loc]] <- merge_enrichment_results(list(c12_frame(c12_terms)))$Term
  }
  skip_if(length(seen) < 2, "only one collation is available on this host")
  for (nm in names(seen)) {
    expect_identical(seen[[nm]], seen[[1]], info = paste("collation:", nm))
  }
})
