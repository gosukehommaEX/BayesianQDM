#' Decision from the Go and NoGo criteria
#'
#' @param met_go,met_nogo Logical vectors: whether the Go and NoGo criteria
#'   are met.
#'
#' @return A character vector of \code{"Go"}, \code{"NoGo"}, \code{"Gray"}
#'   (neither criterion met) or \code{"Miss"} (both criteria met).
#' @noRd
decision_label <- function(met_go, met_nogo) {
  ifelse(met_go, ifelse(met_nogo, "Miss", "Go"),
         ifelse(met_nogo, "NoGo", "Gray"))
}
