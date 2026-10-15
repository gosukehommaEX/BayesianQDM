#' Check the favorable direction of the treatment effect
#'
#' @param direction Character vector whose elements are \code{"greater"}
#'   (larger effects are favorable) or \code{"less"} (smaller effects are
#'   favorable).
#' @param n_endpoints Number of endpoints (1 or 2).
#'
#' @return A character vector of length \code{n_endpoints}.
#' @noRd
check_direction <- function(direction, n_endpoints) {
  ok <- is.character(direction) && !anyNA(direction) &&
    length(direction) %in% unique(c(1L, n_endpoints)) &&
    all(direction %in% c("greater", "less"))
  if (!ok) {
    stop(sprintf("'direction' must be \"greater\" or \"less\"%s.",
                 if (n_endpoints == 2L) {
                   ", or one of them for each of the two endpoints"
                 } else {
                   ""
                 }), call. = FALSE)
  }
  rep_len(direction, n_endpoints)
}
