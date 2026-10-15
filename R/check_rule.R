#' Check the thresholds, regions and rule of a decision criterion
#'
#' For one endpoint, the thresholds are single numbers and the regions and
#' the rule are not used. For two endpoints with \code{rule = "joint"}, the
#' thresholds are single numbers and the Go and NoGo regions are unions of
#' the regions of \code{\link{qdm_regions}}. With \code{rule = "and"} or
#' \code{"or"}, the criteria are based on the marginal probabilities of each
#' endpoint; the thresholds may differ between the endpoints, and the regions
#' are not used.
#'
#' @param rule One of \code{"joint"}, \code{"and"} or \code{"or"}.
#' @param n_endpoints Number of endpoints (1 or 2).
#' @param gamma_go,gamma_nogo Thresholds.
#' @param go_regions,nogo_regions Region numbers.
#' @param regions_given Logical; whether the regions were specified by the
#'   user.
#' @param n_regions Number of regions (9 or 4).
#'
#' @return A list with elements \code{rule}, \code{gamma_go},
#'   \code{gamma_nogo} (of length 2 for the marginal rules), \code{go} and
#'   \code{nogo} (sorted region numbers, or \code{NULL}).
#' @noRd
check_rule <- function(rule, n_endpoints, gamma_go, gamma_nogo, go_regions,
                       nogo_regions, regions_given, n_regions) {
  marginal <- n_endpoints == 2L && rule != "joint"
  len <- if (marginal) NULL else 1L
  check_numeric(gamma_go, "gamma_go", len = len, lower = 0, upper = 1,
                lower_open = TRUE, upper_open = TRUE)
  check_numeric(gamma_nogo, "gamma_nogo", len = len, lower = 0, upper = 1,
                lower_open = TRUE, upper_open = TRUE)
  if (n_endpoints == 1L) {
    return(list(rule = "joint", gamma_go = gamma_go, gamma_nogo = gamma_nogo,
                go = NULL, nogo = NULL))
  }
  if (!marginal) {
    regions <- check_regions(go_regions, nogo_regions, n_regions = n_regions)
    return(list(rule = "joint", gamma_go = gamma_go, gamma_nogo = gamma_nogo,
                go = regions$go, nogo = regions$nogo))
  }
  if (regions_given) {
    stop("'go_regions' and 'nogo_regions' are used only with rule = \"joint\".",
         call. = FALSE)
  }
  if (length(gamma_go) > 2L || length(gamma_nogo) > 2L) {
    stop(sprintf(paste("With rule = \"%s\", 'gamma_go' and 'gamma_nogo' must",
                       "have length 1 or 2 (one threshold per endpoint)."),
                 rule), call. = FALSE)
  }
  list(rule = rule, gamma_go = rep_len(gamma_go, 2L),
       gamma_nogo = rep_len(gamma_nogo, 2L), go = NULL, nogo = NULL)
}
