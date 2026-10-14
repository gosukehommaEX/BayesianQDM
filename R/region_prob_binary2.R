#' Probabilities of regions for two binary endpoints
#'
#' Computes, for pairs of treatment and control Dirichlet distributions, the
#' probabilities of unions of the regions of the criterion: posterior
#' probabilities by Monte Carlo simulation with \code{nMC} draws shared by
#' all pairs (Yamaguchi et al., submitted, equation A.2), or posterior
#' predictive probabilities computed exactly (equation A.4).
#'
#' @param criterion An object of class \code{"qdm_criterion"} for two
#'   endpoints.
#' @param alpha_t,alpha_c Matrices of Dirichlet parameters (one row per
#'   distribution).
#' @param pair_t,pair_c Integer vectors of row indices (from 1) of the pairs.
#' @param sets List of integer vectors of region numbers.
#' @param nMC Number of Monte Carlo draws (posterior criteria).
#'
#' @return A numeric matrix with one row per pair and one column per element
#'   of \code{sets}.
#' @noRd
region_prob_binary2 <- function(criterion, alpha_t, alpha_c, pair_t, pair_c,
                                sets, nMC) {
  masks <- vapply(sets, region_mask, integer(1))
  storage.mode(alpha_t) <- "double"
  storage.mode(alpha_c) <- "double"
  pair_t <- as.integer(pair_t) - 1L
  pair_c <- as.integer(pair_c) - 1L
  if (criterion$type == "posterior") {
    cpp_post_regions_binary2(alpha_t, alpha_c, pair_t, pair_c,
                             as.numeric(criterion$theta_tv),
                             as.numeric(criterion$theta_mav), masks,
                             as.integer(nMC))
  } else {
    cpp_pred_regions_binary2(alpha_t, alpha_c, pair_t, pair_c,
                             criterion$m_t, criterion$m_c,
                             as.numeric(criterion$theta_null), masks)
  }
}
