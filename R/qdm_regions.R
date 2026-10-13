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
#' @param prob Either \code{"posterior"} (nine regions) or \code{"predictive"}
#'   (four regions).
#'
#' @return A character matrix of region labels (printed without quotes), with
#'   the conditions on endpoint 1 as column names and those on endpoint 2 as
#'   row names.
#'
#' @examples
#' qdm_regions("posterior")
#' qdm_regions("predictive")
#' @export
qdm_regions <- function(prob = c("posterior", "predictive")) {
  prob <- match.arg(prob)
  if (prob == "posterior") {
    out <- matrix(paste0("R", 1:9), nrow = 3L,
                  dimnames = list(c("theta2 > TV2", "MAV2 < theta2 <= TV2",
                                    "theta2 <= MAV2"),
                                  c("theta1 > TV1", "MAV1 < theta1 <= TV1",
                                    "theta1 <= MAV1")))
  } else {
    out <- matrix(paste0("R", 1:4), nrow = 2L,
                  dimnames = list(c("theta2 > NULL2", "theta2 <= NULL2"),
                                  c("theta1 > NULL1", "theta1 <= NULL1")))
  }
  noquote(out)
}
