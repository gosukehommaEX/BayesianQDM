#' Labels of the Go and NoGo events of a criterion
#'
#' @param criterion An object of class \code{"qdm_criterion"}.
#' @param digits Number of significant digits for numeric values.
#'
#' @return A character vector with elements \code{go} and \code{nogo}, such
#'   as \code{"P(theta > 0.2 | data)"}.
#' @noRd
event_labels <- function(criterion, digits = 4L) {
  if (criterion$n_endpoints == 2L) {
    go <- sprintf("P(regions {%s} | data)",
                  paste(criterion$go_regions, collapse = ", "))
    nogo <- sprintf("P(regions {%s} | data)",
                    paste(criterion$nogo_regions, collapse = ", "))
  } else if (criterion$type == "posterior") {
    go <- sprintf("P(theta > %s | data)",
                  format_num(criterion$theta_tv, digits))
    nogo <- sprintf("P(theta <= %s | data)",
                    format_num(criterion$theta_mav, digits))
  } else {
    null <- format_num(criterion$theta_null, digits)
    go <- sprintf("P(future theta > %s | data)", null)
    nogo <- sprintf("P(future theta <= %s | data)", null)
  }
  c(go = go, nogo = nogo)
}
