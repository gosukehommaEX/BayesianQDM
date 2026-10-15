#' Probabilities of the Go and NoGo criteria from region probabilities
#'
#' @param p Numeric matrix with one column per set of regions: the Go sets
#'   followed by the NoGo sets (see \code{region_sets()}).
#' @param n_go Number of Go sets (1 for \code{rule = "joint"}, 2 for the
#'   marginal rules).
#' @param dim Optional dimensions given to each column (for example the
#'   numbers of treatment and control outcomes).
#'
#' @return A list with elements \code{go} and \code{nogo}: a vector (or array
#'   with dimensions \code{dim}) for one set, or a list of two for the
#'   marginal rules.
#' @noRd
split_sets <- function(p, n_go, dim = NULL) {
  col <- function(k) {
    v <- p[, k]
    if (!is.null(dim)) {
      v <- array(v, dim = dim)
    }
    v
  }
  pick <- function(ks) if (length(ks) == 1L) col(ks) else lapply(ks, col)
  list(go = pick(seq_len(n_go)), nogo = pick(n_go + seq_len(n_go)))
}
