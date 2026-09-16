# tests/testthat/helper-artifacts.R
# ===========================================================================
# T3-01(d) (WAVE 0) -- THE ARTIFACT TAG REGISTRY AND THE EMISSION / COMPARISON
# HELPER.  Both teams import this file verbatim.
#
# "Artifact identity" is the acceptance mechanism for every assertion in
# Waves 2-5.  EVERY TAG NAMES ITS ENTRY POINT AND ITS FULL ARGUMENT SET,
# because a tag that does not pin the entry point is not an identity at all:
# a cluster()-level distance_matrix on the demo pair is 3403 x 3403 while a
# runRichCluster()-level one on the frozen fixture is 580 x 580, and asserting
# identical() across the two is unsatisfiable by construction.
#
# NEVER CITE AN ARTIFACT BY WAVE NAME ("the Wave-2 artifact") -- cite the tag.
#
# Public surface (SIGNATURES ARE FROZEN -- see the note above rc_emit()):
#   rc_emit(tag, object, entry_point, overwrite = FALSE)   write fixtures/<tag>.rds
#   rc_expect_identical(tag, object, entry_point, args = NULL)  compare against it
#   rc_artifact_registry(tag)                      resolve a tag -> its record
#   rc_artifact_families()                         the registry table itself
#   rc_artifact_object(tag) / rc_artifact_record(tag) / rc_artifact_exists(tag)
#
# entry_point is REQUIRED on both rc_emit() and rc_expect_identical().  T3-01 V6
# requires rc_expect_identical() to ERROR (not report a mismatch) when the stored
# entry point differs from the caller's, and a caller-declared entry point is the
# only thing the helper can compare "the caller's" against.  A defaulted or
# omitted entry_point would make the guard bypassable by omission, which is the
# opposite of "impossible to write by accident".
# ===========================================================================


# --- Where artifacts live --------------------------------------------------

#' Directory holding the frozen fixture and every emitted artifact.
#'
#' Resolution order: RC_FIXTURE_DIR env var, then testthat::test_path("fixtures"),
#' then "fixtures" relative to the working directory (which is tests/testthat
#' during a testthat run).
rc_fixture_dir <- function() {
  env <- Sys.getenv("RC_FIXTURE_DIR", unset = "")
  if (nzchar(env)) return(env)
  if (requireNamespace("testthat", quietly = TRUE)) {
    p <- tryCatch(testthat::test_path("fixtures"), error = function(e) NULL)
    if (!is.null(p) && dir.exists(p)) return(p)
  }
  "fixtures"
}

#' File backing a given artifact tag.
rc_artifact_path <- function(tag) {
  file.path(rc_fixture_dir(), paste0(rc_check_tag_string(tag), ".rds"))
}

#' TRUE when the artifact for this tag has been emitted.
rc_artifact_exists <- function(tag) file.exists(rc_artifact_path(tag))

#' TRUE when this tag's object can be GENERATED at check time.
#'
#' SPEC-RC-004, author ruling 2026-08-26.  Only the v110-* family qualifies.
#' v110 records what the CURRENT, corrected build produces, so it is
#' regenerable by construction -- which is exactly what the frozen v102-*
#' historical family is not (REQ-RC004-011, spec.md section 4.1: regenerating
#' a 1.0.2 characterization from corrected code turns every assertion against
#' it into `current == current`).  Generating v110 keeps the three blocks
#' SPEC-RC-009 re-anchored to that family inside the SHIPPED suite, so they run
#' under R CMD check instead of skipping.
rc_artifact_generatable <- function(tag) {
  isTRUE(tryCatch(identical(rc_artifact_registry(tag)$family, "v110"),
                  error = function(e) FALSE))
}

#' Generate the record for a generatable tag, memoised per R process.
#'
#' Errors -- never skips -- on any tag outside the generatable set, and on any
#' generator failure (REQ-RC004-007 / REQ-RC004-015).  The record is built from
#' the registry exactly as rc_emit() builds it, so rc_expect_identical()'s
#' entry-point and registry-drift guards bind on a generated record precisely as
#' they do on a stored one.  It is flagged `generated = TRUE` so a reader can
#' tell a check-time record from a frozen one.
rc_artifact_generate <- function(tag) {
  tag <- rc_check_tag_string(tag)
  if (!rc_artifact_generatable(tag)) {
    stop(sprintf(paste0(
      "SPEC-RC-004: artifact '%s' is NOT generatable at check time.  Only the ",
      "v110-* family is -- it records what the current build produces.  A ",
      "v102-* tag records the 1.0.2 build and CANNOT be regenerated from ",
      "corrected code (REQ-RC004-011): doing so would turn every assertion ",
      "against it into `current == current`."), tag), call. = FALSE)
  }
  if (is.null(.rc_artifact_memo[[tag]])) {
    reg <- rc_artifact_registry(tag)
    .rc_artifact_memo[[tag]] <- list(
      rc_artifact_version = 1L,
      tag             = tag,
      entry_point     = reg$entry_point,
      component       = reg$component,
      build           = reg$build,
      args            = reg$args,
      note            = reg$note,
      object          = rc_v110_recompute(tag),
      generated       = TRUE,
      created_at      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      r_version       = R.version.string,
      package_version = tryCatch(as.character(utils::packageVersion("richCluster")),
                                 error = function(e) NA_character_)
    )
  }
  .rc_artifact_memo[[tag]]
}

.rc_artifact_memo <- new.env(parent = emptyenv())

#' Skip the calling test when this artifact has not been emitted.
#'
# SPEC-RC-004 REQ-RC004-021 -- DECLARED SKIP EXCEPTION.
# This stays a skip, and that is deliberate.  After SPEC-RC-004 M4, every
# caller of this function that names a v102-* tag lives in a test-src-*.R
# file, and .Rbuildignore excludes those from the tarball -- so in the shipped
# configuration no reachable caller can produce a skip here (verified:
# acceptance.md AC-001 enumerates the shipped skips and this is not among
# them).  The remaining shipped callers name v110-* tags, which are
# GENERATABLE, and the guard returns without skipping for those.  In the
# SOURCE tree it still guards a genuine absence: a v102-* artifact the
# developer has not emitted.  Converting it to an error would fail the suite
# for a developer with an incomplete fixture corpus, which is a different
# defect, not a fix.
rc_skip_if_no_artifact <- function(tag) {
  if (rc_artifact_generatable(tag)) return(invisible(TRUE))
  testthat::skip_if_not(rc_artifact_exists(tag),
                        sprintf("artifact '%s' has not been emitted", tag))
}


# --- Controlled vocabularies ----------------------------------------------

# The five entry points any registry tag can name.
RC_ENTRY_POINTS <- c(
  "runRichCluster()",
  "runDavidClustering()",
  "david_cluster()",
  "cluster()",
  "merge_enrichment_results()"
)

RC_METRICS  <- c("kappa", "jaccard", "dice")            # dice is added by D-02 / T0-05
RC_LINKAGES <- c("single", "complete", "average", "ward")

# Argument-value descriptors reused across families.  These describe the DATA
# arguments; the scalar arguments are written out literally per family below.
RC_ARG_FROZEN_TERMS <- paste0(
  "rc_fixture_frozen()$Term -- the frozen 580-term vector from ",
  "tests/testthat/fixtures/fx-input.rds")
RC_ARG_FROZEN_GENES <- paste0(
  "rc_fixture_frozen()$GeneID -- the frozen 580 GeneID strings from ",
  "tests/testthat/fixtures/fx-input.rds (N = 1819 distinct genes)")
RC_ARG_DEMO_PAIR <- paste0(
  "list(d1, d2) where ",
  "d1 = read.delim(system.file(\"extdata\", \"HF36wk_vs_HF12wk.txt\", package = \"richCluster\")) and ",
  "d2 = read.delim(system.file(\"extdata\", \"HF36wk_vs_WT12wk.txt\", package = \"richCluster\"))")


# --- Entry-point normalisation --------------------------------------------

#' Canonicalise an entry-point string.
#'
#' Accepts "runDavidClustering", "runDavidClustering()" and
#' "richCluster:::runDavidClustering()"; errors on anything outside
#' RC_ENTRY_POINTS so a typo cannot silently become a new entry point.
rc_normalise_entry_point <- function(entry_point) {
  if (!is.character(entry_point) || length(entry_point) != 1L || is.na(entry_point)) {
    stop("entry_point must be a single non-NA character string, one of: ",
         paste(RC_ENTRY_POINTS, collapse = ", "), call. = FALSE)
  }
  x <- trimws(entry_point)
  x <- sub("^[A-Za-z.][A-Za-z0-9._]*:::?", "", x)   # drop a pkg:: / pkg::: prefix
  x <- sub("\\(\\s*\\)$", "", x)                    # drop a trailing ()
  x <- paste0(trimws(x), "()")
  if (!(x %in% RC_ENTRY_POINTS)) {
    stop(sprintf("unknown entry point '%s'.  Known entry points: %s",
                 entry_point, paste(RC_ENTRY_POINTS, collapse = ", ")), call. = FALSE)
  }
  x
}


# --- The registry ----------------------------------------------------------
#
# Six SHAPES.  Each shape fixes an entry point, the component of that entry
# point's return value the tag names, and the argument set:
#
#   dm-<metric>-<linkage>       runRichCluster()          $distance_matrix
#   clusters-<metric>-<linkage> runRichCluster()          $all_clusters
#   cl-dm-<metric>-<linkage>    cluster()                 $distance_matrix
#   david                       runDavidClustering()      whole return
#   dc                          david_cluster()           whole return
#   merge                       merge_enrichment_results() whole return
#
# Six FAMILIES.  Each family fixes the build state the artifact was taken on
# and which shapes it admits:
#
#   v102-*         pre-fix numeric behaviour    dm, clusters, david, dc, merge
#                  (src/ unmodified; the record's derived package_version field
#                   carries the DESCRIPTION version -- SPEC-RC-005 finding B3)
#                  (plus the one-off v102-cl-nopadj)
#   w1-*           end of Wave 1 (after T2-05)  all six shapes
#   w2-*           end of Wave 3 (after T1-15)  all six shapes
#   t0-<nn>-*      immediately after T0-<nn>    dm, clusters, cl-dm, dc, merge
#   pre-<item>-*   immediately before <item>    all six shapes
#   v110-*         corrected DAVID path         david, dc
#                  (immediately after SPEC-RC-009's DS-03/DS-04 fix; the
#                   RE-FREEZE POINT for the DAVID path -- v102-david/v102-dc
#                   remain the frozen historical record and current output
#                   must now DIFFER from them)
#
# Both intra-wave families (t0-<nn>-*, pre-<item>-*) are emitted ON DEMAND, NOT
# EXHAUSTIVELY: emit exactly the shapes the next item's verifications assert
# against, so no tag is defined-but-never-emitted -- the defect that made w2-*
# unusable in the previous revision.
#
# The t0-<nn>-* family exists because Wave 1 is a chain: w1-* is not emitted
# until Wave 2, i.e. AFTER every Wave-1 item including T2-05, so a Wave-1 item
# cannot assert against a w1- tag -- the tag does not exist when its test runs.
# The pre-<item>-* family exists because a no-op proof is a before/after
# identity, and for an item landing mid-wave the "before" state is frequently
# not one of the frozen wave tags (T1-07 moves the exported diagonal and T1-10
# moves merged_df$GeneID, both BEFORE T1-16 lands, so neither w1-dm-* nor
# v102-dc is the correct baseline for the hoist).
#
# The w1-cl-* family exists because T0-06 V1/V2 assert a cluster()-LEVEL
# identity while the anchor fixture is a runRichCluster()-level object.  They
# are different objects and each needs its own tag.  Do NOT resolve that by
# re-anchoring the fixture at cluster() level -- that reintroduces C1.
#
# The david_cluster() entry point has its own family (v102-dc / w1-dc / w2-dc),
# separate from the runDavidClustering() family (v102-david / w1-david /
# w2-david).  They are different objects, and an assertion naming one while
# comparing against the other's tag is exactly the entry-point mismatch that
# rc_expect_identical() errors on.

#' The registry table: one row per (family, shape) pair.
#'
#' <metric> and <linkage> are placeholders; rc_artifact_registry() substitutes
#' the concrete values parsed out of a concrete tag.
rc_artifact_families <- function() {
  fam <- list(
    v102 = list(build  = paste0("pre-fix numeric behaviour: src/ unmodified ",
                                "(see the record's package_version field for ",
                                "the DESCRIPTION version at emission time)"),
                shapes = c("dm", "clusters", "david", "dc", "merge")),
    w1   = list(build  = "end of Wave 1 (after T2-05)",
                shapes = c("dm", "clusters", "cl-dm", "david", "dc", "merge")),
    w2   = list(build  = "end of Wave 3 (emitted in one pass immediately after T1-15)",
                shapes = c("dm", "clusters", "cl-dm", "david", "dc", "merge")),
    t0   = list(build  = "immediately after item T0-<nn> lands (Wave 1 chain)",
                shapes = c("dm", "clusters", "cl-dm", "dc", "merge")),
    pre  = list(build  = "immediately before item <item> lands",
                shapes = c("dm", "clusters", "cl-dm", "david", "dc", "merge")),
    v110 = list(build  = "corrected DAVID path, canonical merge order (SPEC-RC-009 DS-03/DS-04 plus C12)",
                shapes = c("david", "dc"))
  )
  prefix <- c(v102 = "v102", w1 = "w1", w2 = "w2", t0 = "t0-<nn>", pre = "pre-<item>",
              v110 = "v110")
  rows <- do.call(rbind, lapply(names(fam), function(f) {
    do.call(rbind, lapply(fam[[f]]$shapes, function(s) {
      sh <- rc_shape_spec(s, metric = "<metric>", linkage = "<linkage>")
      data.frame(family      = f,
                 tag_pattern = paste0(prefix[[f]], "-", sh$shape_pattern),
                 shape       = s,
                 entry_point = sh$entry_point,
                 component   = sh$component,
                 build       = fam[[f]]$build,
                 arguments   = rc_args_to_string(sh$args),
                 stringsAsFactors = FALSE)
    }))
  }))
  # The one non-parameterised tag: a cluster()-level snapshot on frames that do
  # not exist in Wave 0, so T3-03 cannot snapshot it -- T1-17 emits it.
  sh <- rc_shape_spec("cl-nopadj")
  rows <- rbind(rows,
                data.frame(family      = "v102",
                           tag_pattern = "v102-cl-nopadj",
                           shape       = "cl-nopadj",
                           entry_point = sh$entry_point,
                           component   = sh$component,
                           build       = paste0("pre-fix numeric behaviour: src/ unmodified ",
                                                "(see the record's package_version field for ",
                                                "the DESCRIPTION version at emission time)",
                                                "; EMITTED BY T1-17, not by T3-03"),
                           arguments   = rc_args_to_string(sh$args),
                           stringsAsFactors = FALSE))
  rows[order(rows$family, rows$shape), ]
}

#' Entry point, component and argument set for one shape.
rc_shape_spec <- function(shape, metric = NULL, linkage = NULL) {
  switch(shape,
    "dm" = list(
      shape_pattern = "dm-<metric>-<linkage>",
      entry_point = "runRichCluster()",
      component   = "$distance_matrix (580 x 580 on the frozen fixture)",
      args = list(terms          = RC_ARG_FROZEN_TERMS,
                  geneIDs        = RC_ARG_FROZEN_GENES,
                  distanceMetric = metric,
                  distanceCutoff = "0.5",
                  linkageMethod  = linkage,
                  linkageCutoff  = "0.5")),
    "clusters" = list(
      shape_pattern = "clusters-<metric>-<linkage>",
      entry_point = "runRichCluster()",
      component   = "$all_clusters",
      args = list(terms          = RC_ARG_FROZEN_TERMS,
                  geneIDs        = RC_ARG_FROZEN_GENES,
                  distanceMetric = metric,
                  distanceCutoff = "0.5",
                  linkageMethod  = linkage,
                  linkageCutoff  = "0.5")),
    "cl-dm" = list(
      shape_pattern = "cl-dm-<metric>-<linkage>",
      entry_point = "cluster()",
      component   = paste0("$distance_matrix (3403 x 3403 on the demo pair -- a ",
                           "DIFFERENT object from the runRichCluster()-level dm-* shape)"),
      args = list(enrichment_results = RC_ARG_DEMO_PAIR,
                  df_names           = "unset (NULL)",
                  min_terms          = "5 (documented default)",
                  min_value          = "0.1 (documented default)",
                  distance_metric    = metric,
                  distance_cutoff    = "0.5",
                  linkage_method     = linkage,
                  linkage_cutoff     = "0.5",
                  filter_on          = "unset (argument does not exist at 1.0.2; added by T1-17)",
                  gene_delim         = "unset (argument does not exist at 1.0.2; added by T0-06)")),
    "david" = list(
      shape_pattern = "david",
      entry_point = "runDavidClustering()",
      component   = "the whole runDavidClustering() return value",
      args = list(terms                    = RC_ARG_FROZEN_TERMS,
                  geneIDs                  = RC_ARG_FROZEN_GENES,
                  similarityThreshold      = "0.5",
                  initialGroupMembership   = "3",
                  finalGroupMembership     = "3",
                  multipleLinkageThreshold = "0.5")),
    "dc" = list(
      shape_pattern = "dc",
      entry_point = "david_cluster()",
      component   = paste0("the whole david_cluster() return value -- the exported wrapper ",
                           "merges internally and applies NO p-value filter, so its term set ",
                           "is the full 3403 and its return additionally carries merged_df, ",
                           "df_list, df_names, cluster_options and cluster_df"),
      args = list(enrichment_results         = RC_ARG_DEMO_PAIR,
                  df_names                   = "unset (NULL)",
                  similarity_threshold       = "0.5",
                  initial_group_membership   = "3",
                  final_group_membership     = "3",
                  multiple_linkage_threshold = "0.5")),
    "merge" = list(
      shape_pattern = "merge",
      entry_point = "merge_enrichment_results()",
      component   = "the whole merge_enrichment_results() return value (3403 rows, Term-sorted)",
      args = list(enrichment_results = RC_ARG_DEMO_PAIR)),
    "cl-nopadj" = list(
      shape_pattern = "cl-nopadj",
      entry_point = "cluster()",
      component   = "the whole cluster() return value",
      args = list(enrichment_results = paste0(
                    "the purpose-built Term / GeneID / Pvalue, NO Padj frames defined by ",
                    "T1-17 V3b (they do not exist in Wave 0, which is why T3-03 cannot ",
                    "snapshot this tag)"),
                  df_names        = "unset (NULL)",
                  min_terms       = "5 (documented default)",
                  min_value       = "0.1",
                  distance_metric = "kappa (documented default)",
                  distance_cutoff = "0.5 (documented default)",
                  linkage_method  = "average (documented default)",
                  linkage_cutoff  = "0.5 (documented default)",
                  filter_on       = "unset")),
    stop(sprintf("unknown artifact shape '%s'", shape), call. = FALSE))
}

rc_args_to_string <- function(args) {
  paste(sprintf("%s = %s", names(args), unlist(args, use.names = FALSE)), collapse = "; ")
}

rc_check_tag_string <- function(tag) {
  if (!is.character(tag) || length(tag) != 1L || is.na(tag) || !nzchar(tag)) {
    stop("tag must be a single non-empty character string", call. = FALSE)
  }
  tag
}

#' Resolve a concrete tag against the registry.
#'
#' Returns a list with the tag's entry point, the component of that entry
#' point's return value it names, the build state it was taken on, and its FULL
#' argument set.  Errors on any tag the registry does not define -- an
#' unregistered tag is a typo, and a typo that silently becomes a new baseline
#' is the defect this registry exists to prevent.
rc_artifact_registry <- function(tag) {
  tag <- rc_check_tag_string(tag)

  if (identical(tag, "v102-cl-nopadj")) {
    sh <- rc_shape_spec("cl-nopadj")
    return(rc_registry_record(tag, "v102", sh,
      build = paste0("pre-fix numeric behaviour: src/ unmodified ",
                     "(see the record's package_version field for ",
                     "the DESCRIPTION version at emission time)",
                     "; EMITTED BY T1-17, not by T3-03"),
      note  = paste0("cluster()-level snapshot on frames T1-17 V3b builds. Those frames ",
                     "do not exist in Wave 0, so T3-03 cannot snapshot them.")))
  }

  m <- regmatches(tag, regexec("^(v102|w1|w2|v110)-(.+)$", tag))[[1]]
  if (length(m) == 3L) {
    family <- m[2]; shape_str <- m[3]
    build <- switch(family,
      v102 = paste0("pre-fix numeric behaviour: src/ unmodified ",
                    "(see the record's package_version field for ",
                    "the DESCRIPTION version at emission time)"),
      w1   = "end of Wave 1 (after T2-05)",
      w2   = "end of Wave 3 (emitted in one pass immediately after T1-15)",
      v110 = "corrected DAVID path, canonical merge order (SPEC-RC-009 DS-03/DS-04 plus C12)")
    allowed <- if (family == "v110") c("david", "dc")
               else if (family == "v102") c("dm", "clusters", "david", "dc", "merge")
               else c("dm", "clusters", "cl-dm", "david", "dc", "merge")
    sh <- rc_parse_shape(shape_str, tag, allowed, family)
    return(rc_registry_record(tag, family, sh, build = build))
  }

  m <- regmatches(tag, regexec("^t0-([0-9]{2})-(.+)$", tag))[[1]]
  if (length(m) == 3L) {
    sh <- rc_parse_shape(m[3], tag, c("dm", "clusters", "cl-dm", "dc", "merge"), "t0")
    return(rc_registry_record(tag, "t0", sh,
      build = sprintf("immediately after item T0-%s lands (Wave 1 chain)", m[2]),
      note  = paste0("Wave 1 is a chain: w1-* is not emitted until Wave 2, so a Wave-1 item ",
                     "asserts against the tag its predecessor emitted at the same entry point.")))
  }

  m <- regmatches(tag, regexec("^pre-(T[0-4]-[0-9]{2}[a-z]?)-(.+)$", tag))[[1]]
  if (length(m) == 3L) {
    sh <- rc_parse_shape(m[3], tag, c("dm", "clusters", "cl-dm", "david", "dc", "merge"), "pre")
    return(rc_registry_record(tag, "pre", sh,
      build = sprintf("immediately before item %s lands", m[2]),
      note  = paste0("A no-op proof is a before/after identity, and mid-wave the 'before' ",
                     "state is frequently not one of the frozen wave tags.")))
  }

  stop(sprintf(paste0("tag '%s' is not in the T3-01(d) artifact registry.\n",
                      "Registered patterns: v102-<shape>, v102-cl-nopadj, w1-<shape>, ",
                      "w2-<shape>, t0-<nn>-<shape>, pre-<item>-<shape>, v110-<shape>\n",
                      "with <shape> in dm-<metric>-<linkage>, clusters-<metric>-<linkage>, ",
                      "cl-dm-<metric>-<linkage>, david, dc, merge."), tag), call. = FALSE)
}

rc_parse_shape <- function(shape_str, tag, allowed, family = NA_character_) {
  if (shape_str %in% c("david", "dc", "merge")) {
    if (!(shape_str %in% allowed)) {
      stop(sprintf("shape '%s' is not defined for tag '%s'; allowed shapes: %s",
                   shape_str, tag, paste(allowed, collapse = ", ")), call. = FALSE)
    }
    return(rc_shape_spec(shape_str))
  }
  m <- regmatches(shape_str, regexec("^(cl-dm|dm|clusters)-([a-z]+)-([a-z]+)$", shape_str))[[1]]
  if (length(m) != 4L) {
    stop(sprintf("cannot parse shape '%s' of tag '%s'", shape_str, tag), call. = FALSE)
  }
  shape <- m[2]; metric <- m[3]; linkage <- m[4]
  if (!(shape %in% allowed)) {
    stop(sprintf("shape '%s' is not defined for tag '%s'; allowed shapes: %s",
                 shape, tag, paste(allowed, collapse = ", ")), call. = FALSE)
  }
  if (!(metric %in% RC_METRICS)) {
    stop(sprintf("unknown distance metric '%s' in tag '%s'; known: %s",
                 metric, tag, paste(RC_METRICS, collapse = ", ")), call. = FALSE)
  }
  # SPEC-RC-005, finding A7: dice is ADDED by SPEC-RC-002.  It does not exist in
  # the stock build the v102 family characterises, so there is no stock
  # behaviour for a v102-* tag to capture.  Valid in post-fix families only.
  #
  # The wording says "the stock build", not "the stock <version> build": finding
  # B3 (REQ-RC-084) forbids asserting a version number in prose, and AC-RC-084
  # greps this whole file for the old phrase.  REQ-RC-079's own text likewise
  # reads "the stock build".
  if (identical(family, "v102") && identical(metric, "dice")) {
    stop(sprintf(paste0("tag '%s' names metric 'dice' in the v102 family, but ",
                        "'dice' does not exist in the stock build -- it is ",
                        "ADDED by SPEC-RC-002.  There is no stock behaviour for ",
                        "'dice' to characterise.  Use dice only in the post-fix ",
                        "families (w1, w2, t0, pre)."), tag), call. = FALSE)
  }
  if (!(linkage %in% RC_LINKAGES)) {
    stop(sprintf("unknown linkage method '%s' in tag '%s'; known: %s",
                 linkage, tag, paste(RC_LINKAGES, collapse = ", ")), call. = FALSE)
  }
  rc_shape_spec(shape, metric = metric, linkage = linkage)
}

rc_registry_record <- function(tag, family, sh, build, note = NA_character_) {
  list(tag         = tag,
       family      = family,
       shape       = sub("-<metric>-<linkage>$", "", sh$shape_pattern),
       entry_point = sh$entry_point,
       component   = sh$component,
       build       = build,
       args        = sh$args,
       note        = note)
}


# --- Emission and comparison ----------------------------------------------

#' Emit an artifact under a registry tag.
#'
#' Writes tests/testthat/fixtures/<tag>.rds.  The stored record carries the
#' object TOGETHER WITH its entry point, its component, the build state and its
#' full argument set, so the artifact is self-describing on disk.
#'
#' entry_point is the entry point the CALLER used to produce `object`; it must
#' match the registry's entry point for `tag`, or this errors.
#'
#' Refuses to overwrite an existing artifact whose object differs unless
#' overwrite = TRUE.  Re-emitting identical content is a no-op.
#'
#' @section FROZEN CONTRACT (SPEC-RC-005, reviewer finding A1):
#' This FOUR-argument signature is the specification.  An earlier plan document
#' specified rc_emit(tag, object); that document no longer exists and the
#' shipped implementation supersedes it.  The ruling is recorded here, and
#' asserted by tests/testthat/test-instrument-contracts.R, so it cannot drift
#' again.
#'
#' entry_point carries NO DEFAULT, deliberately.  T3-01 V6 requires
#' rc_expect_identical() to ERROR on an entry-point mismatch, and a
#' caller-declared entry point is the only thing the helper can compare "the
#' caller's" against.  A defaulted or omitted entry_point would make the guard
#' bypassable by omission -- the opposite of "impossible to write by accident".
#'
#' Do not add, remove, reorder, or default these arguments.
rc_emit <- function(tag, object, entry_point, overwrite = FALSE) {
  reg <- rc_artifact_registry(tag)
  ep  <- rc_normalise_entry_point(entry_point)
  if (!identical(ep, reg$entry_point)) {
    stop(sprintf(paste0("ENTRY-POINT MISMATCH while emitting '%s'.\n",
                        "  registry entry point : %s\n",
                        "  caller entry point   : %s\n",
                        "Two artifacts produced at different entry points are never ",
                        "comparable; emit under the tag whose entry point you used."),
                 tag, reg$entry_point, ep), call. = FALSE)
  }
  dir <- rc_fixture_dir()
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  path <- rc_artifact_path(tag)
  if (file.exists(path) && !isTRUE(overwrite)) {
    prev <- readRDS(path)
    if (!identical(prev$object, object)) {
      stop(sprintf(paste0("Refusing to overwrite artifact '%s': the object differs from the ",
                          "one already frozen at '%s'.  Silently re-baselining a frozen tag ",
                          "invalidates every identity asserted against it.  Pass ",
                          "overwrite = TRUE if that is genuinely intended."), tag, path),
           call. = FALSE)
    }
    return(invisible(path))
  }
  record <- list(
    rc_artifact_version = 1L,
    tag             = tag,
    entry_point     = reg$entry_point,
    component       = reg$component,
    build           = reg$build,
    args            = reg$args,
    note            = reg$note,
    object          = object,
    created_at      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    r_version       = R.version.string,
    package_version = tryCatch(as.character(utils::packageVersion("richCluster")),
                               error = function(e) NA_character_)
  )
  saveRDS(record, path)
  invisible(path)
}

#' Read the stored record for a tag (object plus its self-description).
rc_artifact_record <- function(tag) {
  path <- rc_artifact_path(tag)
  if (file.exists(path)) return(readRDS(path))          # stored file wins
  # SPEC-RC-004 (author ruling 2026-08-26): generate the v110-* family at check
  # time rather than skipping.  Stored-file precedence above keeps the source
  # tree bit-for-bit on the frozen record; only the tarball reaches this branch.
  if (rc_artifact_generatable(tag)) return(rc_artifact_generate(tag))
  stop(sprintf(paste0("artifact '%s' has not been emitted (expected '%s').\n",
                      "Its emission point is: %s"),
               tag, path, rc_artifact_registry(tag)$build), call. = FALSE)
}

#' Read just the stored object for a tag.
#'
#' This is what a cross-team identical() comparison uses: the record also
#' carries created_at / r_version, which legitimately differ between teams.
rc_artifact_object <- function(tag) rc_artifact_record(tag)$object

#' Assert that `object` is identical() to the artifact frozen under `tag`.
#'
#' ERRORS (rather than reporting a mismatch) when the stored entry point differs
#' from the caller's, or when the stored record has drifted from the registry.
#' That is what makes an identical() assertion across two different entry points
#' impossible to write by accident (T3-01 V6).
rc_expect_identical <- function(tag, object, entry_point, args = NULL) {
  reg    <- rc_artifact_registry(tag)
  ep     <- rc_normalise_entry_point(entry_point)
  record <- rc_artifact_record(tag)

  if (!identical(ep, record$entry_point)) {
    stop(sprintf(paste0("ENTRY-POINT MISMATCH for artifact '%s'.\n",
                        "  stored entry point : %s\n",
                        "  caller entry point : %s\n",
                        "  stored arguments   : %s\n",
                        "These are different objects and are never comparable.  Cite the tag ",
                        "whose entry point you actually called."),
                 tag, record$entry_point, ep, rc_args_to_string(record$args)), call. = FALSE)
  }
  if (!identical(record$entry_point, reg$entry_point) ||
      !identical(record$args, reg$args)) {
    stop(sprintf(paste0("REGISTRY DRIFT for artifact '%s': the stored record no longer ",
                        "matches the T3-01(d) registry.\n",
                        "  stored : %s | %s\n  registry: %s | %s\n",
                        "Re-emit the artifact at its declared emission point (%s)."),
                 tag, record$entry_point, rc_args_to_string(record$args),
                 reg$entry_point, rc_args_to_string(reg$args), reg$build), call. = FALSE)
  }
  if (!is.null(args) && !identical(args, record$args)) {
    stop(sprintf(paste0("ARGUMENT-SET MISMATCH for artifact '%s'.\n",
                        "  stored : %s\n  caller : %s"),
                 tag, rc_args_to_string(record$args), rc_args_to_string(args)), call. = FALSE)
  }

  testthat::expect_true(
    identical(object, record$object),
    info = sprintf(paste0("artifact identity FAILED for tag '%s'\n",
                          "  entry point : %s\n  component   : %s\n",
                          "  build       : %s\n  arguments   : %s"),
                   tag, record$entry_point, record$component,
                   record$build, rc_args_to_string(record$args)))
}
