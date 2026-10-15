#' Labels of the Go and NoGo events of a criterion
#'
#' @param criterion An object of class \code{"qdm_criterion"}.
#' @param digits Number of significant digits for numeric values.
#'
#' @return A list with character vectors \code{go} and \code{nogo}, such as
#'   \code{"P(theta > 0.2 | data)"}; for two endpoints with a marginal rule,
#'   each has one label per endpoint.
#' @noRd
event_labels <- function(criterion, digits = 4L) {
  num <- function(x) vapply(x, format_num, character(1), digits = digits)
  if (criterion$n_endpoints == 2L && criterion$rule == "joint") {
    go <- sprintf("P(regions {%s} | data)",
                  paste(criterion$go_regions, collapse = ", "))
    nogo <- sprintf("P(regions {%s} | data)",
                    paste(criterion$nogo_regions, collapse = ", "))
  } else {
    # One endpoint: "theta"; two endpoints (marginal rules): "theta1" and
    # "theta2"
    sym <- if (criterion$n_endpoints == 2L) c("theta1", "theta2") else "theta"
    if (criterion$type == "posterior") {
      go <- sprintf("P(%s > %s | data)", sym, num(criterion$theta_tv))
      nogo <- sprintf("P(%s <= %s | data)", sym, num(criterion$theta_mav))
    } else {
      null <- num(criterion$theta_null)
      go <- sprintf("P(future %s > %s | data)", sym, null)
      nogo <- sprintf("P(future %s <= %s | data)", sym, null)
    }
  }
  list(go = unname(go), nogo = unname(nogo))
}
