#' Regions used by the Go and NoGo criteria of two endpoints
#'
#' @param criterion An object of class \code{"qdm_criterion"} for two
#'   endpoints.
#'
#' @return A list with elements \code{go} and \code{nogo}, each a list of
#'   integer vectors of region numbers: one vector (the union of regions) for
#'   \code{rule = "joint"}, and one per endpoint for the marginal rules
#'   (Yamaguchi et al., submitted, equations 4, 5, 9 and 10).
#' @noRd
region_sets <- function(criterion) {
  if (criterion$rule == "joint") {
    return(list(go = list(criterion$go_regions),
                nogo = list(criterion$nogo_regions)))
  }
  if (criterion$type == "posterior") {
    list(go = list(1:3, c(1L, 4L, 7L)), nogo = list(7:9, c(3L, 6L, 9L)))
  } else {
    list(go = list(1:2, c(1L, 3L)), nogo = list(3:4, c(2L, 4L)))
  }
}
