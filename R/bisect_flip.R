#' Locate the points at which a logical function changes its value
#'
#' Vectorized bisection: for each element, the logical function \code{pred}
#' has different values at \code{lo} and \code{hi}, and the interval is
#' halved \code{n_iter} times, keeping the half whose ends have different
#' values.
#'
#' @param pred Function of a vector of condition indices and a vector of
#'   values that returns a logical vector.
#' @param idx Condition indices.
#' @param lo,hi Ends of the intervals (\code{lo < hi}).
#' @param n_iter Number of bisection steps.
#'
#' @return The midpoints of the final intervals, of width
#'   \code{(hi - lo) / 2^n_iter}.
#' @noRd
bisect_flip <- function(pred, idx, lo, hi, n_iter = 40L) {
  if (length(idx) == 0L) {
    return(numeric(0))
  }
  at_lo <- pred(idx, lo)
  for (k in seq_len(n_iter)) {
    mid <- (lo + hi) / 2
    same <- pred(idx, mid) == at_lo
    lo <- ifelse(same, mid, lo)
    hi <- ifelse(same, hi, mid)
  }
  (lo + hi) / 2
}
