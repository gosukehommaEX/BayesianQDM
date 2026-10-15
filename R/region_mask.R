#' Bit mask of a set of regions
#'
#' @param regions Integer vector of region numbers (from 1).
#'
#' @return An integer with bit \code{r - 1} set for each region \code{r}.
#' @noRd
region_mask <- function(regions) {
  as.integer(sum(2 ^ (unique(regions) - 1)))
}
