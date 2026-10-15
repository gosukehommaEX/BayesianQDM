#' Print a decision boundary
#'
#' Prints the decision intervals of \code{\link{qdm_boundary}} as up to
#' three tables, each with one row per condition (for example, the number of
#' control responders or the treatment group size) and one column per
#' decision. The first table gives the range of the observed quantity
#' (numbers of responders or events, or observed means), the second the
#' corresponding observed response rates or milestone survival
#' probabilities, and the last the observed treatment effects. For example,
#' the entry \code{"10-20"} in the column \code{Go} means that 10 to 20
#' responders lead to a Go decision. The object itself
#' is a data frame with one row per interval (see \code{\link{qdm_boundary}}).
#'
#' @param x An object of class \code{"qdm_boundary"}.
#' @param digits Number of significant digits for observed means, milestone
#'   survival probabilities, effects and conditions.
#' @param n Maximum number of conditions shown in each table.
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
#' b <- qdm_boundary(des)
#' print(b, n = 5)
#' # The intervals as a data frame
#' head(as.data.frame(b))
#' @export
print.qdm_boundary <- function(x, digits = 4L, n = 20L, ...) {
  if (is.null(attr(x, "axes"))) {
    # Attributes are lost when the rows are subset
    print(as.data.frame(x))
    return(invisible(x))
  }
  check_numeric(digits, "digits", lower = 1, integer = TRUE)
  check_numeric(n, "n", lower = 1, integer = TRUE)
  cat(format(x, digits = digits, n = n), sep = "\n")
  invisible(x)
}
