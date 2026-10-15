#' Distribution of the mean of a group with two continuous endpoints
#'
#' The posterior distribution of the mean is
#' \eqn{t_{\nu - 1}(\mu, \Lambda / (\kappa (\nu - 1)))}, and the posterior
#' predictive distribution of the mean of \eqn{m} future subjects is
#' \eqn{t_{\nu - 1}(\mu, (1 / \kappa + 1 / m) \Lambda / (\nu - 1))}.
#'
#' @param post Parameters of the group (see \code{posterior_continuous2()}).
#' @param type \code{"posterior"} or \code{"predictive"}.
#' @param m Number of future subjects (predictive only).
#'
#' @return A list with elements \code{loc} (matrix with two columns),
#'   \code{scale} (matrix with the columns \code{V11}, \code{V12} and
#'   \code{V22}) and \code{df} (a number).
#' @noRd
mean_dist_continuous2 <- function(post, type, m = NULL) {
  df <- post$nu - 1
  f <- if (type == "posterior") {
    1 / (post$kappa * df)
  } else {
    (1 / post$kappa + 1 / m) / df
  }
  scale <- f * post$Lambda
  dimnames(scale) <- list(NULL, c("V11", "V12", "V22"))
  list(loc = post$mu, scale = scale, df = df)
}
