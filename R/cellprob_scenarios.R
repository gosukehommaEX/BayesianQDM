#' Cell probabilities of scenarios for two binary endpoints
#'
#' Applies \code{\link{cellprob_binary2}} to each scenario. A scenario whose
#' correlation is outside the feasible range for its response rates is not
#' evaluable and gives missing probabilities.
#'
#' @param pi1,pi2,rho Numeric vectors of response rates and correlations
#'   (checked), one element per scenario.
#'
#' @return A numeric matrix with one row per scenario and the columns
#'   \code{"00"}, \code{"01"}, \code{"10"} and \code{"11"}; rows of
#'   infeasible scenarios are \code{NA}.
#' @noRd
cellprob_scenarios <- function(pi1, pi2, rho) {
  out <- t(vapply(seq_along(pi1), function(s) {
    tryCatch(cellprob_binary2(pi1[s], pi2[s], rho[s]),
             error = function(e) rep(NA_real_, 4L))
  }, numeric(4)))
  dimnames(out) <- list(NULL, c("00", "01", "10", "11"))
  out
}
