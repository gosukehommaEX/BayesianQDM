#' Decision with the highest probability in each scenario
#'
#' For the overlay heat map of two endpoints: the decision with the highest
#' probability in each scenario (ties, up to 1e-10, go to the first decision
#' of \code{decisions}), its probability and the step of that probability:
#' 1 below 0.5, 2 from 0.5 to below 0.6, 3 from 0.6, 4 from 0.7 and 5 from
#' 0.8 (a probability below a step by less than 1e-10 counts as reaching
#' it).
#'
#' @param x A data frame with a column of probabilities for each element of
#'   \code{decisions}.
#' @param decisions Character vector of decisions, in the order used for
#'   ties.
#'
#' @return A data frame with the columns \code{decision} (\code{NA} for a
#'   scenario that is not evaluable), \code{prob} and \code{step}.
#' @noRd
dominant_decision <- function(x, decisions) {
  p <- as.matrix(x[, decisions, drop = FALSE])
  k <- apply(p, 1L, function(r) {
    if (anyNA(r)) NA_integer_ else which(r >= max(r) - 1e-10)[1L]
  })
  prob <- p[cbind(seq_len(nrow(p)), k)]
  data.frame(decision = decisions[k], prob = prob,
             step = findInterval(prob + 1e-10, c(0.5, 0.6, 0.7, 0.8)) + 1L,
             stringsAsFactors = FALSE)
}
