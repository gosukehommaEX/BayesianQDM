#' Probabilities of unions of regions for two continuous endpoints
#'
#' @param criterion A criterion for two endpoints.
#' @param dist_t,dist_c Distributions of the treatment and control means (see
#'   \code{mean_dist_continuous2()}), with the same number of rows or one
#'   row.
#' @param sets List of integer vectors of region numbers of the criterion.
#' @param method \code{"MM"} or \code{"MC"}.
#' @param nMC Number of Monte Carlo draws (MC only).
#'
#' @return A numeric matrix with one row per pair of distributions and one
#'   column per element of \code{sets}.
#' @noRd
region_prob_continuous2 <- function(criterion, dist_t, dist_c, sets, method,
                                    nMC = 10000L) {
  if (criterion$type == "posterior") {
    tv <- as.numeric(criterion$theta_tv)
    mav <- as.numeric(criterion$theta_mav)
  } else {
    tv <- as.numeric(criterion$theta_null)
    mav <- tv
  }
  masks <- vapply(sets, function(s) {
    region_mask(region_cells(s, criterion$type))
  }, integer(1))
  n <- max(nrow(dist_t$loc), nrow(dist_c$loc))
  rows <- function(x) {
    x <- x[rep_len(seq_len(nrow(x)), n), , drop = FALSE]
    storage.mode(x) <- "double"
    x
  }
  args <- list(rows(dist_t$loc), rows(dist_t$scale),
               rep_len(as.numeric(dist_t$df), n), rows(dist_c$loc),
               rows(dist_c$scale), rep_len(as.numeric(dist_c$df), n), tv,
               mav, masks)
  if (method == "MC") {
    return(do.call(cpp_regions_continuous2_mc,
                   c(args, list(as.integer(nMC)))))
  }
  do.call(cpp_regions_continuous2_mm, args)
}
