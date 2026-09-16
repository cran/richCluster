# tests/testthat/test-cpp-units.R
# ===========================================================================
# T3-04 (WAVE 0) -- DRIVER FOR THE STANDALONE C++ UNIT HARNESS.
#
# The harness itself is tools/cpp-units/rc_cpp_units.cpp.  This file compiles
# it against the package's own src/*.cpp with R's own toolchain, runs it, and
# reports PASS/FAIL PER ASSERTION (T3-04 V1) -- one testthat expectation per
# RCUNIT line, so a failure names the exact case rather than an exit code.
#
# ***  THESE ASSERTIONS ARE EXPECTED TO BE RED ON THE STOCK BUILD.  ***
#
# That is the specification, not a defect.  T3-04's NULL-ACTION CHECK:
#
#     "V2 run against the stock build must FAIL on the kappa and jaccard hand
#      cases (returning the intersection-undercounted kappa and 0.2 instead of
#      0.25) and must PASS on the trivially correct average() case.  A harness
#      that passes everything on stock code is not testing the right entry
#      points."
#
# Wave 0 exists because fix_conflict FC-9 is binding: no test in richCluster
# 1.0.2 can distinguish before from after for any numeric correction in this
# plan.  A test that CAN distinguish is therefore red until its owning item
# lands.  Each RCUNIT line carries that owning item, and the failure message
# repeats it.  Do not "fix" the harness to make these green; land the item.
#
#   K1   T0-01 V2          kappa closed form         -> green when T0-01 lands
#   J1   T0-05 V1          jaccard closed form       -> green when T0-05 lands
#   J2   T0-05 V1          dice closed form          -> green when T0-05 lands
#   L1   T3-04 null-action average(), disjoint pair  -> GREEN ON STOCK ALREADY
#   L2a  T0-04 V2          single()  >= average()    -> green when T0-04 lands
#   L2b  T0-04 V2          average() >= complete()   -> green when T0-04 lands
#   L3   T0-04 V4          negative linkage survives -> green when T0-04 lands
#   L4a  T0-04 V5          empty cross-product       -> green when T0-04 lands
#   L4b  T0-04 V5          no NaN                    -> GREEN ON STOCK ALREADY
#   L4c  T0-04 V5          no linkage == 100         -> GREEN ON STOCK ALREADY
#   W1   T4-01             ward singleton anchor     -> GREEN ON STOCK ALREADY
#   W2   T4-01             ward 2x2 closed form      -> green when T4-01 lands
#   W3   T4-01             ward sepSq clamp          -> green when T4-01 lands
#   W4   T4-01             ward empty cross-product  -> GREEN ON STOCK ALREADY
#   W5   T4-01             ward negative similarity  -> GREEN ON STOCK ALREADY
#   W6   T4-01             ward size penalty (value) -> green when T4-01 lands
#   W6b  T4-01             ward size penalty (ineq)  -> green when T4-01 lands
#   K2   N2                kappa of an identical pair -> green when N2 lands
#   K3   N2                kappa == jaccard == dice   -> green when N2 lands
#   K4   N2                identical pair, wider N    -> GREEN ON STOCK ALREADY
#   W7   N1                ward on an OVERLAPPING pair -> green when N1 lands
#   W8   N1                ward on a NESTED pair       -> green when N1 lands
#   W9   N1                ward vs the ESS definition  -> green when N1 lands
#   W9b  N1                W9 covered both branches    -> GREEN ON STOCK ALREADY
#   S1   T3-04 exposure    splitStringToVector       -> GREEN ON STOCK ALREADY
#   S2   T3-04 exposure    splitStringToUnorderedSet -> GREEN ON STOCK ALREADY
#   S3   T3-04 exposure    countUniqueElements       -> GREEN ON STOCK ALREADY
#
# T3-04 V3 (a -fsanitize=address,undefined build) is opt-in via
# RC_CPP_UNITS_SANITIZE=1 rather than unconditional: requiring a
# sanitizer-capable toolchain of every check host is not portable, and the
# sanitized run is an implementer/CI verification.  With the flag set, the
# harness is rebuilt with the sanitizers and the run is additionally asserted
# to emit no sanitizer diagnostic.
# ===========================================================================


# --- Locating the sources --------------------------------------------------
#
# tools/ is in the tarball but is NOT installed, so the harness is reachable
# from the SOURCE tree only.  Two layouts matter:
#
#   source tree / devtools::test()   <pkg>/tests/testthat   -> ../..
#   R CMD check on a tarball         <pkg>.Rcheck/tests/testthat
#                                    -> ../../00_pkg_src/<pkg>
#
# RC_CPP_SRC overrides both (point it at the package root).  The override is
# VALIDATED like any other candidate: a wrong RC_CPP_SRC must produce a clear
# skip, not a confusing compiler error a hundred lines later.
rc_cpp_root_ok <- function(cand) {
  file.exists(file.path(cand, "src", "DistanceMetric.cpp")) &&
    file.exists(file.path(cand, "tools", "cpp-units", "rc_cpp_units.cpp"))
}

rc_cpp_pkg_root <- function() {
  env <- Sys.getenv("RC_CPP_SRC", unset = "")
  if (nzchar(env)) {
    if (rc_cpp_root_ok(env)) return(normalizePath(env, mustWork = FALSE))
    # SPEC-RC-005, disposition row 2: the operator DECLARED an intent to run the
    # harness and named a wrong root.  `fatal` routes this to the fail field --
    # skipping here would silently discard that declaration.
    return(structure(NA_character_, fatal = TRUE, reason = sprintf(
      paste0("RC_CPP_SRC is set to '%s', but that directory has no ",
             "src/DistanceMetric.cpp and tools/cpp-units/rc_cpp_units.cpp. ",
             "Point it at the package root, or unset it."), env)))
  }
  base <- tryCatch(testthat::test_path("."), error = function(e) ".")
  candidates <- c(
    file.path(base, "..", ".."),
    Sys.glob(file.path(base, "..", "..", "00_pkg_src", "*")),
    file.path(base, "..", "..", "..")
  )
  for (cand in candidates) {
    if (rc_cpp_root_ok(cand)) return(normalizePath(cand, mustWork = FALSE))
  }
  structure(NA_character_, reason = paste0(
    "package sources not found: tools/ is not installed, so this test needs ",
    "the source tree (or set RC_CPP_SRC to the package root)"))
}

#' One value from `R CMD config`, "" when unset or unavailable.
rc_cpp_rconfig <- function(var) {
  out <- suppressWarnings(tryCatch(
    system2(file.path(R.home("bin"), "R"), c("CMD", "config", var),
            stdout = TRUE, stderr = FALSE),
    error = function(e) character(0)))
  if (length(out) == 0L) "" else paste(out, collapse = " ")
}

#' Split a config string into tokens, dropping empties.
rc_cpp_tokens <- function(x) {
  x <- trimws(x)
  if (!nzchar(x)) return(character(0))
  t <- strsplit(x, "[[:space:]]+")[[1]]
  t[nzchar(t)]
}


# --- Build + run -----------------------------------------------------------

#' Compile and run the harness; return its RCUNIT lines parsed into a
#' data.frame, plus the raw output and the exit status.
#'
#' Uses R's own CXX17 so the harness is built by the same toolchain that built
#' the package, and compiles the package's REAL src/*.cpp -- the harness must
#' exercise the shipped code, never a copy of it.
rc_cpp_units_run <- function(sanitize = FALSE) {
  root <- rc_cpp_pkg_root()
  if (is.na(root)) {
    reason <- attr(root, "reason")
    if (isTRUE(attr(root, "fatal"))) {
      # Row 2: the operator named a root explicitly and named it wrong.
      return(list(status = NA_integer_, out = character(0), df = NULL,
                  skip = NA_character_, fail = reason))
    }
    # Row 1: tools/ is genuinely unreachable from an installed package.
    return(list(status = NA_integer_, out = character(0), df = NULL,
                skip = reason, fail = NA_character_))
  }

  # R reports CXX17 as a bare compiler, but CXX as compiler-plus-standard-flag
  # ("x86_64-linux-gnu-g++ -std=gnu++20").  system2() takes a single command,
  # so the first token is the command and any remainder becomes leading args.
  cc  <- rc_cpp_tokens(rc_cpp_rconfig("CXX17"))
  std <- rc_cpp_rconfig("CXX17STD")
  flg <- rc_cpp_rconfig("CXX17FLAGS")
  if (!length(cc)) {
    cc  <- rc_cpp_tokens(rc_cpp_rconfig("CXX"))
    flg <- rc_cpp_rconfig("CXXFLAGS")
    std <- ""
  }
  if (!length(cc)) {
    # Row 3: the toolchain cannot be interrogated at all.
    return(list(status = NA_integer_, out = character(0), df = NULL,
                skip = "no C++ compiler reported by R CMD config",
                fail = NA_character_))
  }
  cxx <- cc[1L]
  cc_args <- if (length(cc) > 1L) cc[-1L] else character(0)

  dir <- file.path(tempdir(), paste0("rc-cpp-units-", if (sanitize) "asan" else "plain"))
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  bin <- file.path(dir, "rc_cpp_units")

  srcs <- file.path(root, "src",
                    c("DistanceMetric.cpp", "LinkageMethod.cpp", "StringUtils.cpp"))
  main <- file.path(root, "tools", "cpp-units", "rc_cpp_units.cpp")

  extra <- if (sanitize) {
    c("-fsanitize=address,undefined", "-fno-omit-frame-pointer", "-O1", "-g")
  } else {
    character(0)
  }

  args <- c(cc_args,
            rc_cpp_tokens(std), rc_cpp_tokens(flg),
            extra,
            paste0("-I", shQuote(file.path(root, "src"))),
            "-o", shQuote(bin),
            shQuote(main), shQuote(srcs))

  build <- suppressWarnings(system2(cxx, args, stdout = TRUE, stderr = TRUE))
  build_status <- attr(build, "status")
  if (!is.null(build_status) && build_status != 0L) {
    # Row 4 -- FINDING B1.  The instrument is PRESENT and BROKEN.  Skipping here
    # is what let the instrument silently disappear; it is a failure.
    return(list(status = NA_integer_, out = build, df = NULL,
                skip = NA_character_,
                fail = paste0("the C++ unit harness FAILED TO COMPILE.  The ",
                              "instrument is present but broken -- this is a ",
                              "failure, not an absence.\n",
                              paste(build, collapse = "\n"))))
  }

  out <- suppressWarnings(system2(bin, character(0), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  status <- if (is.null(status)) 0L else as.integer(status)

  lines <- grep("^RCUNIT\t", out, value = TRUE)
  df <- NULL
  if (length(lines)) {
    parts <- strsplit(sub("^RCUNIT\t", "", lines), "\t", fixed = TRUE)
    keep <- vapply(parts, length, integer(1)) == 7L
    parts <- parts[keep]
    df <- data.frame(
      id       = vapply(parts, `[`, character(1), 1L),
      status   = vapply(parts, `[`, character(1), 2L),
      observed = vapply(parts, `[`, character(1), 3L),
      expected = vapply(parts, `[`, character(1), 4L),
      tol      = vapply(parts, `[`, character(1), 5L),
      item     = vapply(parts, `[`, character(1), 6L),
      desc     = vapply(parts, `[`, character(1), 7L),
      stringsAsFactors = FALSE)
  }
  list(status = status, out = out, df = df,
       skip = NA_character_, fail = NA_character_)
}

# The harness is compiled once per process and reused by both test blocks.
rc_cpp_units_cached <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) cache <<- rc_cpp_units_run(sanitize = FALSE)
    cache
  }
})


# --- V1 / V2 -- per-assertion reporting ------------------------------------

test_that("V1/V2: the C++ unit harness builds, runs, and reports per assertion", {
  res <- rc_cpp_units_cached()
  # FINDING B1: a present-but-broken instrument FAILS.  This assertion must come
  # BEFORE the skip -- a skip evaluated first would exit the block silently.
  expect_true(is.na(res$fail), info = res$fail)
  testthat::skip_if(!is.na(res$skip), res$skip)

  # V1 -- it ran at all, and every declared case reported a line.
  expect_false(is.null(res$df))
  expect_true(nrow(res$df) > 0)
  expect_setequal(res$df$id,
                  c("K1", "K2", "K3", "K4", "J1", "J2", "L1", "L2a", "L2b", "L3",
                    "L4a", "L4b", "L4c",
                    "W1", "W2", "W3", "W4", "W5", "W6", "W6b",
                    "W7", "W8", "W9", "W9b",
                    "S1", "S2", "S3",
                    "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "L5"))

  # V2 -- one expectation per assertion.  A failure names the case, the
  # observed and expected values, and the item that turns it green.
  for (i in seq_len(nrow(res$df))) {
    r <- res$df[i, ]
    testthat::expect_equal(
      r$status, "PASS",
      info = sprintf(
        paste0("[%s] %s\n  owning item : %s\n  observed    : %s\n",
               "  expected    : %s  (tol %s)\n",
               "  This case is RED until its owning item lands -- see the ",
               "header of tests/testthat/test-cpp-units.R."),
        r$id, r$desc, r$item, r$observed, r$expected, r$tol))
  }
})


# --- V3 -- sanitizer build (opt-in) ----------------------------------------

test_that("V3: the harness builds and runs clean under -fsanitize=address,undefined", {
  testthat::skip_if_not(
    identical(Sys.getenv("RC_CPP_UNITS_SANITIZE"), "1"),
    paste0("T3-04 V3 needs a sanitizer-capable toolchain, which not every ",
           "check host has.  Set RC_CPP_UNITS_SANITIZE=1 to run it."))

  res <- rc_cpp_units_run(sanitize = TRUE)
  # FINDING B1 (REQ-RC-052): the operator opted in; a broken sanitizer build is
  # a failure, not an absence.  Assert before the skip.
  expect_true(is.na(res$fail), info = res$fail)
  testthat::skip_if(!is.na(res$skip), res$skip)

  # "without diagnostics ON THE PASSING CASES": the run must reach its summary
  # line and emit no sanitizer diagnostic.  Cases that legitimately FAIL on
  # this build are still failures -- they are not sanitizer findings.
  expect_true(any(grepl("^RCUNIT-SUMMARY\t", res$out)))
  diagnostics <- grep("(runtime error:|AddressSanitizer|LeakSanitizer|UndefinedBehaviorSanitizer)",
                      res$out, value = TRUE)
  expect_equal(diagnostics, character(0))

  # The sanitized build must agree with the plain build case for case; a
  # divergence would mean the plain build's result depends on UB.
  plain <- rc_cpp_units_cached()
  expect_true(is.na(plain$fail), info = plain$fail)
  testthat::skip_if(!is.na(plain$skip), plain$skip)
  expect_equal(res$df[order(res$df$id), c("id", "status")],
               plain$df[order(plain$df$id), c("id", "status")],
               ignore_attr = TRUE)
})
