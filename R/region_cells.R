#' Cells of the nine regions of two continuous endpoints
#'
#' The probabilities of two continuous endpoints are computed on the nine
#' regions of a posterior criterion (Yamaguchi et al., submitted, Table 1).
#' The four regions of a predictive criterion (Table 2) are the regions 1,
#' 3, 7 and 9 of this grid with the target and minimum acceptable values
#' both equal to the null value.
#'
#' @param regions Region numbers of the criterion.
#' @param type \code{"posterior"} or \code{"predictive"}.
#'
#' @return The numbers of the corresponding regions of the nine-region grid.
#' @noRd
region_cells <- function(regions, type) {
  if (type == "posterior") {
    return(as.integer(regions))
  }
  c(1L, 3L, 7L, 9L)[regions]
}
