#' Posterior Dirichlet parameters for two binary endpoints
#'
#' Combines a Dirichlet prior, external data borrowed through a power prior
#' and hypothetical control data with the counts of the current data. With
#' prior Dir(alpha), external counts x_e with weight a_0, hypothetical counts
#' x_h and current counts x, the posterior distribution is
#' Dir(alpha + a_0 x_e + x_h + x) (Yamaguchi et al., submitted, Appendix
#' A.1.5).
#'
#' @param prior An object of class \code{"qdm_prior_dirichlet"}.
#' @param counts Integer matrix of current counts with one row per outcome
#'   (or a vector of length 4).
#' @param external An object of class \code{"qdm_borrow"} or \code{NULL}.
#' @param hypothetical An object of class \code{"qdm_hypothetical"} with data
#'   for two binary endpoints, or \code{NULL}.
#'
#' @return A numeric matrix of Dirichlet parameters with one row per row of
#'   \code{counts} and the columns \code{"00"}, \code{"01"}, \code{"10"} and
#'   \code{"11"}.
#' @noRd
posterior_binary2 <- function(prior, counts, external = NULL,
                              hypothetical = NULL) {
  counts <- matrix(as.numeric(counts), ncol = 4L)
  add <- prior$alpha
  if (!is.null(external)) {
    add <- add + external$a0 * external$data$counts
  }
  if (!is.null(hypothetical)) {
    add <- add + hypothetical$data$counts
  }
  out <- counts + rep(as.numeric(add), each = nrow(counts))
  dimnames(out) <- list(NULL, c("00", "01", "10", "11"))
  out
}
