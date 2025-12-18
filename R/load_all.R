#' Load all R scripts in subdirectories
#'
#' @name load_all
#' @keywords internal
NULL

source_files <- function(path) {
  files <- list.files(path, pattern = "\\.[Rr]$", full.names = TRUE, recursive = TRUE)
  for (file in files) {
    source(file)
  }
}

# Source all R scripts in /modules (shiny), /scripts (R)
source_files(system.file("modules", package = "richCluster"))
# source_files(system.file("scripts", package = "richCluster"))
