#' Print a decision boundary
#'
#' Prints the number of grid points with each decision and the first rows of
#' the grid.
#'
#' @param x An object of class \code{"qdm_boundary"}.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' des <- design_binary(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_beta(0.5, 0.5), prior_c = prior_beta(0.5, 0.5),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                   gamma_go = 0.5, gamma_nogo = 0.8)
#' )
#' qdm_boundary(des)
#' @export
print.qdm_boundary <- function(x, ...) {
  counts <- table(x$decision)
  cat(sprintf("Decision boundary on a grid of %d point%s (%s)\n", nrow(x),
              if (nrow(x) == 1L) "" else "s",
              paste(names(counts), counts, sep = ": ", collapse = ", ")))
  df <- as.data.frame(x)
  attr(df, "design") <- NULL
  attr(df, "axes") <- NULL
  print(df[seq_len(min(6L, nrow(df))), , drop = FALSE])
  k <- nrow(df) - 6L
  if (k > 0L) {
    cat(sprintf("... %d more row%s\n", k, if (k == 1L) "" else "s"))
  }
  invisible(x)
}
