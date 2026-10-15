#' Validate Go and NoGo region sets
#'
#' Checks that the Go and NoGo region sets for two endpoints are non-empty,
#' contain valid region numbers without duplicates, and are disjoint.
#'
#' @param go_regions,nogo_regions Region numbers.
#' @param n_regions Number of regions (9 for posterior probability, 4 for
#'   predictive probability).
#'
#' @return A list with sorted integer vectors \code{go} and \code{nogo}.
#' @noRd
check_regions <- function(go_regions, nogo_regions, n_regions) {
  check_numeric(go_regions, "go_regions", len = NULL, lower = 1,
                upper = n_regions, integer = TRUE)
  check_numeric(nogo_regions, "nogo_regions", len = NULL, lower = 1,
                upper = n_regions, integer = TRUE)
  go <- as.integer(round(go_regions))
  nogo <- as.integer(round(nogo_regions))
  if (anyDuplicated(go) > 0L) {
    stop("'go_regions' must not contain duplicated values.", call. = FALSE)
  }
  if (anyDuplicated(nogo) > 0L) {
    stop("'nogo_regions' must not contain duplicated values.", call. = FALSE)
  }
  if (length(intersect(go, nogo)) > 0L) {
    stop("'go_regions' and 'nogo_regions' must not share any region.",
         call. = FALSE)
  }
  list(go = sort(go), nogo = sort(nogo))
}
