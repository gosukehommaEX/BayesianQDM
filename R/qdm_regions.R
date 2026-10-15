#' Region numbers for two endpoints
#'
#' Shows how the regions used by the decision criteria for two endpoints are
#' numbered (Yamaguchi et al., submitted, Tables 1 and 2).
#'
#' For posterior probability, the target value (TV) and the minimum acceptable
#' value (MAV) of each endpoint divide the plane of \eqn{(\theta_1, \theta_2)}
#' into nine regions. Columns correspond to endpoint 1 and rows to endpoint 2.
#' For posterior predictive probability, the null values divide the plane into
#' four regions.
#'
#' For an endpoint on which smaller effects are favorable
#' (\code{direction = "less"} in the criterion), the inequalities of that
#' endpoint are reversed, so that region 1 is always the most favorable one:
#' for example, \code{"theta2 < TV2"}, \code{"TV2 <= theta2 < MAV2"} and
#' \code{"theta2 >= MAV2"}.
#'
#' @param prob Either \code{"posterior"} (nine regions) or \code{"predictive"}
#'   (four regions).
#' @param direction \code{"greater"} (the default) or \code{"less"}, as in the
#'   criterion: a single value for both endpoints or one value per endpoint.
#'
#' @return A character matrix of region labels (printed without quotes), with
#'   the conditions on endpoint 1 as column names and those on endpoint 2 as
#'   row names.
#'
#' @references
#' Yamaguchi, Y., Hu, H., Huang, Y., Han, C., Homma, G. and Wang, A.
#' (submitted). Bivariate Bayesian quantitative decision-making for two
#' endpoints in rare disease drug development.
#'
#' @examples
#' qdm_regions("posterior")
#' qdm_regions("predictive")
#' # Endpoint 2 is a safety endpoint (smaller effects are favorable)
#' qdm_regions("posterior", direction = c("greater", "less"))
#' @export
qdm_regions <- function(prob = c("posterior", "predictive"),
                        direction = "greater") {
  prob <- match.arg(prob)
  direction <- check_direction(direction, 2L)
  labels <- function(k) {
    th <- paste0("theta", k)
    if (prob == "posterior") {
      tv <- paste0("TV", k)
      mav <- paste0("MAV", k)
      if (direction[k] == "greater") {
        c(paste(th, ">", tv), paste(mav, "<", th, "<=", tv),
          paste(th, "<=", mav))
      } else {
        c(paste(th, "<", tv), paste(tv, "<=", th, "<", mav),
          paste(th, ">=", mav))
      }
    } else {
      null <- paste0("NULL", k)
      if (direction[k] == "greater") {
        c(paste(th, ">", null), paste(th, "<=", null))
      } else {
        c(paste(th, "<", null), paste(th, ">=", null))
      }
    }
  }
  n <- if (prob == "posterior") 3L else 2L
  out <- matrix(paste0("R", seq_len(n * n)), nrow = n,
                dimnames = list(labels(2L), labels(1L)))
  noquote(out)
}
