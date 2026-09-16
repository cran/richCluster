# tests/testthat/helper-characterization.R
# ===========================================================================
# T3-03 (WAVE 0) -- THE v1.0.2 CHARACTERIZATION SNAPSHOT.
#
# ***  THE VALUES THIS FILE FREEZES ARE KNOWN WRONG.  THEY ARE NOT A TARGET.  ***
#
# This file is the one-shot generator for the `v102-*` artifact family: what
# richCluster 1.0.2 ACTUALLY computed, captured on the STOCK build before any
# numeric code was edited.  It exists for exactly one purpose -- to prove which
# later changes are NO-OPS.  An item that claims to be NON-BREAKING proves it by
# asserting bit-identity against a frozen registry tag; an item marked
# non-breaking without such an assertion has not been verified.
#
# Nothing in the package may ever be "fixed" to match these numbers.  When a
# later item legitimately MOVES them it does not re-freeze this family; it
# asserts against a different, later tag (t0-<nn>-*, pre-<item>-*, w1-*, w2-*).
#
# STRUCTURE.  rc_v102_recompute() is the single definition of "how the object
# under tag <t> is produced".  Both the emitter (rc_freeze_v102_snapshot) and
# the V1 reproduction test call it, so emission and verification cannot drift
# apart -- a snapshot whose generator is not the thing the test re-runs is not
# a characterization of anything.
#
# THIS GENERATOR ONLY REPRODUCES ON THE STOCK BUILD, BY DESIGN.  From Wave 1
# onward the corrected code returns different objects; that is the whole point.
# The V1 reproduction test in test-characterization.R is therefore gated behind
# RC_V102_STOCK_BUILD=1 and is a Wave-0 gate, not a standing regression test.
#
# Measurement discipline (T3-01): run ONE linkage method per R process.
# rc_freeze_v102_snapshot() accepts a `tags` subset precisely so a caller can
# drive it one tag per process.  Output suppression is done in exactly ONE
# place (rc_v102_recompute, quiet = TRUE) and is never nested: repeated .Call
# plus NESTED sink() in a single process is what produced the spurious
# "non-determinism" the revision-1 critic had to retract.
#
# Depends on helper-artifacts.R (rc_emit, rc_artifact_registry, rc_fixture_dir)
# and helper-fixture.R (rc_fixture_frozen).  testthat sources every helper-*.R
# before any test, so call-time lookup always resolves.
# ===========================================================================


# --- The header every frozen file carries (T3-03 V3) -----------------------

#' The characterization header stamped into every `v102-*` artifact.
#'
#' T3-03 V3: the snapshot is a CHARACTERIZATION, not a target, and every file
#' must say so on its face.  rc_emit_v102() writes this string into the stored
#' record as its first element, so `readRDS(<file>)[[1]]` is the warning.
RC_V102_CHARACTERIZATION_HEADER <- paste0(
  "CHARACTERIZATION SNAPSHOT OF richCluster 1.0.2 -- THESE VALUES ARE KNOWN WRONG.\n",
  "This artifact records what the STOCK 1.0.2 build ACTUALLY computed, not what it\n",
  "should have computed.  It is NOT a target: no code may ever be changed to match\n",
  "it, and it is not evidence that any number in it is correct.  The kappa\n",
  "statistic itself was computed incorrectly in 1.0.2 (std::set_intersection over\n",
  "std::unordered_set undercounts the gene intersection), so on the 580-term anchor\n",
  "fixture the exported matrix carries 141 upper-triangle entries outside the\n",
  "mathematically possible range [-1, 1] and sums to -3611.890440.\n",
  "IT EXISTS FOR EXACTLY ONE PURPOSE: to prove which later changes are NO-OPS.  An\n",
  "item marked non-breaking without a bit-identical assertion against a registered\n",
  "tag has not been verified.  When a later item legitimately MOVES these numbers\n",
  "it does not re-freeze this family -- it asserts against a different, later tag\n",
  "(t0-<nn>-*, pre-<item>-*, w1-*, w2-*).\n",
  "Frozen by T3-03 (Wave 0), on the stock build, before any numeric code was edited.")


# --- What the STOCK build can be asked for ---------------------------------
#
# The stock DistanceMetric::computeDistance implements exactly two metrics and
# throws "unsupported distance metric" on anything else, so "every metric x
# linkage combination" is 2 x 4 = 8 on this build.  `dice` is a registry-valid
# metric (RC_METRICS) but is ADDED by D-02 / T0-05: no v102-*-dice-* tag can
# exist, because the build this family characterises cannot compute one.
# LinkageMethod::computeLinkage dispatches all four linkages (ward is an alias
# of average at 1.0.2 -- that is itself a characterised defect, not a target).
RC_V102_METRICS  <- c("kappa", "jaccard")
RC_V102_LINKAGES <- c("single", "complete", "average", "ward")

#' Every registry tag T3-03 freezes.
#'
#' 8 `v102-dm-*` + 8 `v102-clusters-*` + `v102-merge` + `v102-david` +
#' `v102-dc` = 19.  `v102-cl-nopadj` is deliberately absent: it is a cluster()
#' snapshot on frames that do not exist in Wave 0, and T3-01(d) assigns it to
#' T1-17, not to this item.
rc_v102_tags <- function() {
  grid <- expand.grid(linkage = RC_V102_LINKAGES, metric = RC_V102_METRICS,
                      stringsAsFactors = FALSE)
  c(sprintf("v102-dm-%s-%s", grid$metric, grid$linkage),
    sprintf("v102-clusters-%s-%s", grid$metric, grid$linkage),
    "v102-merge", "v102-david", "v102-dc")
}

#' The demo pair, read exactly as RC_ARG_DEMO_PAIR describes it.
rc_v102_demo_pair <- function() {
  list(
    read.delim(system.file("extdata", "HF36wk_vs_HF12wk.txt", package = "richCluster")),
    read.delim(system.file("extdata", "HF36wk_vs_WT12wk.txt", package = "richCluster"))
  )
}


# --- The single definition of how each tag's object is produced ------------

#' Recompute the object a `v102-*` tag names, from its registry argument set.
#'
#' This is the ONLY place that says how a tag's object is produced.  The
#' emitter and the V1 reproduction test both call it.
#'
#' @param tag a `v102-*` registry tag (see rc_v102_tags()).
#' @param quiet suppress the C++ layer's progress chatter.  Implemented with a
#'   single, never-nested utils::capture.output().
#' @return the object the tag names -- NOT a record; rc_emit() wraps it.
rc_v102_recompute <- function(tag, quiet = TRUE) {
  reg <- rc_artifact_registry(tag)
  if (!identical(reg$family, "v102")) {
    stop(sprintf("rc_v102_recompute() only produces the v102-* family; got '%s'", tag),
         call. = FALSE)
  }

  run <- switch(
    reg$shape,

    # runRichCluster(), frozen 580-term vectors, <metric>, 0.5, <linkage>, 0.5
    "dm" = ,
    "clusters" = function() {
      # Driven by the registry's OWN argument set, not by a second parse of the
      # tag, so the object cannot be produced with arguments the tag does not
      # claim.  rc_shape_spec() fills distanceMetric / linkageMethod with the
      # concrete values parsed out of the tag; the two cutoffs are literals.
      fx <- rc_fixture_frozen()
      res <- richCluster::runRichCluster(fx$Term, fx$GeneID,
                                         reg$args$distanceMetric,
                                         as.numeric(reg$args$distanceCutoff),
                                         reg$args$linkageMethod,
                                         as.numeric(reg$args$linkageCutoff))
      if (identical(reg$shape, "dm")) res$distance_matrix else res$all_clusters
    },

    # runDavidClustering(), frozen 580-term vectors, 0.5 / 3 / 3 / 0.5.
    # Not in NAMESPACE at 1.0.2, hence the ::: -- T3-04a is what makes the
    # DAVID internals observable; this entry point itself is reachable.
    "david" = function() {
      dfx <- rc_david_fixture()
      richCluster:::runDavidClustering(dfx$Term, dfx$GeneID,
                                       dfx$similarity_threshold,
                                       dfx$initial_group_membership,
                                       dfx$final_group_membership,
                                       dfx$multiple_linkage_threshold)
    },

    # david_cluster(), list(d1, d2), documented defaults.  A DIFFERENT object
    # at a DIFFERENT entry point from `david`: the wrapper merges internally
    # and applies NO p-value filter, so its term set is the full 3403.
    "dc" = function() {
      richCluster::david_cluster(rc_v102_demo_pair())
    },

    "merge" = function() {
      richCluster::merge_enrichment_results(rc_v102_demo_pair())
    },

    stop(sprintf("T3-03 does not freeze shape '%s' (tag '%s')", reg$shape, tag),
         call. = FALSE)
  )

  if (isTRUE(quiet)) {
    out <- NULL
    utils::capture.output(out <- run())   # one sink level, never nested
    out
  } else {
    run()
  }
}


# --- Emission --------------------------------------------------------------

#' Emit one `v102-*` artifact and stamp the characterization header on it.
#'
#' Delegates the canonical record (entry point, component, build, argument set,
#' object) to T3-01(d)'s rc_emit(); this only prepends the T3-03 V3 header.
#' The header is an ADDITIONAL field -- record$object, record$entry_point and
#' record$args are untouched, so rc_expect_identical() is unaffected.
rc_emit_v102 <- function(tag, object, entry_point, overwrite = FALSE) {
  path <- rc_emit(tag, object, entry_point, overwrite = overwrite)
  record <- readRDS(path)
  if (!identical(record$rc_characterization_header, RC_V102_CHARACTERIZATION_HEADER)) {
    record <- c(list(rc_characterization_header = RC_V102_CHARACTERIZATION_HEADER),
                record[setdiff(names(record), "rc_characterization_header")])
    saveRDS(record, path)
  }
  invisible(path)
}

#' Freeze the v1.0.2 characterization snapshot.  WAVE 0, STOCK BUILD, run once.
#'
#' Refuses to run against a build whose numeric layer has moved: it is checked
#' by rc_emit()'s own no-silent-rebaseline guard, which errors rather than
#' overwriting an artifact whose object differs.
#'
#' @param tags subset of rc_v102_tags() to emit.  Pass one tag per R process to
#'   honour the T3-01 measurement discipline.
rc_freeze_v102_snapshot <- function(tags = rc_v102_tags(), quiet = TRUE,
                                    overwrite = FALSE) {
  unknown <- setdiff(tags, rc_v102_tags())
  if (length(unknown)) {
    stop("not part of the T3-03 snapshot: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  paths <- character(0)
  for (tag in tags) {
    reg <- rc_artifact_registry(tag)
    obj <- rc_v102_recompute(tag, quiet = quiet)
    paths <- c(paths, rc_emit_v102(tag, obj, reg$entry_point, overwrite = overwrite))
  }
  invisible(paths)
}


# --- The shipped demonstration RDS, captured as-is -------------------------
#
# inst/extdata/cluster_result.rds is a SHIPPED DATA FILE, not the return value
# of an entry point, so it has no T3-01(d) registry tag and deliberately gets
# none: rc_artifact_registry() rejects this name, which is what stops anyone
# citing it as if it were an artifact identity.  It is snapshotted because
# Wave 5 (FC-3 / D-09) regenerates it, and the 1.0.2 file is the only evidence
# of what the package demonstrated to users while it was on CRAN.

rc_v102_shipped_path <- function() {
  file.path(rc_fixture_dir(), "v102-shipped-cluster-result.rds")
}

#' Freeze the shipped inst/extdata/cluster_result.rds as-is.  WAVE 0.
#'
#' Stores the header, the file's md5 (which pins the BYTES) and the
#' deserialised object (which pins the CONTENT).  The 1.0.2 bytes themselves
#' remain recoverable from git; this makes the change detectable from a test.
rc_freeze_v102_shipped_cluster_result <- function(overwrite = FALSE) {
  src <- system.file("extdata", "cluster_result.rds", package = "richCluster")
  if (!nzchar(src)) stop("inst/extdata/cluster_result.rds not found", call. = FALSE)
  dir <- rc_fixture_dir()
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  path <- rc_v102_shipped_path()
  record <- list(
    rc_characterization_header = RC_V102_CHARACTERIZATION_HEADER,
    rc_artifact_version = 1L,
    tag         = "v102-shipped-cluster-result",
    entry_point = NA_character_,
    component   = paste0("the file inst/extdata/cluster_result.rds as shipped with ",
                         "richCluster 1.0.2 -- a stored data file, not the return ",
                         "value of an entry point, and therefore deliberately NOT a ",
                         "T3-01(d) registry tag"),
    build       = paste0("pre-fix numeric behaviour: src/ unmodified ",
                         "(see the record's package_version field for ",
                         "the DESCRIPTION version at emission time)"),
    args        = list(source = "system.file(\"extdata\", \"cluster_result.rds\", package = \"richCluster\")"),
    note        = paste0("Wave 5 (FC-3 / D-09) regenerates this file.  md5 pins the ",
                         "shipped bytes; object pins the shipped content."),
    md5         = unname(tools::md5sum(src)),
    # SPEC-RC-005, finding A4: a checksum without the expression that produced
    # it is not a checksum -- a reader cannot tell WHAT was digested, and a
    # reported md5 that does not reproduce cannot be diagnosed.
    md5_input   = paste0("tools::md5sum(system.file(\"extdata\", ",
                         "\"cluster_result.rds\", package = \"richCluster\"))"),
    file_size   = unname(file.info(src)$size),
    object      = readRDS(src),
    created_at  = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    r_version   = R.version.string,
    package_version = tryCatch(as.character(utils::packageVersion("richCluster")),
                               error = function(e) NA_character_)
  )
  if (file.exists(path) && !isTRUE(overwrite)) {
    prev <- readRDS(path)
    if (!identical(prev$object, record$object) || !identical(prev$md5, record$md5)) {
      stop(sprintf(paste0("Refusing to overwrite '%s': the shipped ",
                          "inst/extdata/cluster_result.rds has changed since it was ",
                          "frozen.  Pass overwrite = TRUE only if re-baselining the ",
                          "1.0.2 demonstration surface is genuinely intended."), path),
           call. = FALSE)
    }
    return(invisible(path))
  }
  saveRDS(record, path)
  invisible(path)
}


# --- Integrity constants (T3-03 V2) ----------------------------------------
#
# The recorded stock observations.  V2 asserts the FROZEN artifacts carry
# these, which is what confirms the snapshot was taken against the right build:
# a snapshot taken accidentally against a partially fixed build shows 0
# out-of-range entries instead of 141.
#
# All are measured on the DOCUMENTED-DEFAULT path (kappa / average) over the
# UPPER TRIANGLE of the exported 580 x 580 matrix, except where stated.  The
# whole matrix carries 862 out-of-range entries = 2 * 141 (symmetric) + 580
# (the -99 diagonal); the upper-triangle count is the one T3-03 V2 quotes.
RC_V102_INTEGRITY <- list(
  n_terms                 = 580L,
  n_genes                 = 1819L,
  dm_upper_sum            = -3611.890440,
  dm_pairs_ge_cutoff      = 98L,     # upper-triangle entries >= 0.5
  dm_outside_range        = 141L,    # upper-triangle entries outside [-1, 1]
  dm_min                  = -6.637792,
  dm_max                  = 1.0,
  dm_diagonal             = -99,     # richCluster::SAME_TERM_DISTANCE
  raw_clusters            = 505L,    # nrow(all_clusters), kappa / average
  david_clusters          = 90L,     # runDavidClustering(), 0.5 / 3 / 3 / 0.5
  david_terms_covered     = 523L,
  dc_clusters             = 557L,    # david_cluster(list(d1, d2)), defaults
  dc_terms_covered        = 2922L,
  merge_rows              = 3403L
)

#' Distinct term indices covered by a `clusters` data frame's TermIndices.
rc_terms_covered <- function(clusters_df) {
  length(unique(unlist(strsplit(clusters_df$TermIndices, ", ", fixed = TRUE))))
}
