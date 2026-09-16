load_cluster_result <- function() {
  path <- system.file("extdata", "cluster_result.rds", package = "richCluster")
  if (!nzchar(path)) {
    # SPEC-RC-004 REQ-RC004-007 / REQ-RC004-015: inst/extdata/cluster_result.rds
    # SHIPS in the tarball, so its absence can never be a legitimate absence --
    # it is a real failure.  An instrument that disappears quietly is the defect
    # class SPEC-RC-004 closes.
    stop("SPEC-RC-004: inst/extdata/cluster_result.rds is missing.  It ships ",
         "in the tarball; its absence is a failure, not a skip.", call. = FALSE)
  }
  readRDS(path)
}
