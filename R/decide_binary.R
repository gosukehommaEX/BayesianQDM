#' Posterior distributions and region probabilities for observed binary data
#'
#' @param design An object of class \code{"qdm_design_binary"}.
#' @param data_t,data_c Observed data of the treatment and control groups;
#'   \code{data_c} is \code{NULL} for a design without a concurrent control
#'   group.
#'
#' @return A list with elements \code{prob_go}, \code{prob_nogo},
#'   \code{posterior_t} and \code{posterior_c} (named vectors of the Beta
#'   parameters \code{a} and \code{b}).
#' @noRd
decide_binary <- function(design, data_t, data_c) {
  if (!inherits(data_t, "qdm_data_binary")) {
    stop("'data_t' must be created by data_binary().", call. = FALSE)
  }
  if (design$control == "concurrent") {
    if (!inherits(data_c, "qdm_data_binary")) {
      stop(paste("'data_c' must be created by data_binary() for a design",
                 "with a concurrent control group."), call. = FALSE)
    }
    post_c <- posterior_binary(design$prior_c, data_c$n, data_c$y,
                               external = design$external_c)
  } else {
    if (!is.null(data_c)) {
      stop(paste("'data_c' must be NULL for a design without a concurrent",
                 "control group."), call. = FALSE)
    }
    post_c <- posterior_binary(design$prior_c, 0, 0,
                               external = design$external_c,
                               hypothetical = design$hypothetical_c)
  }
  post_t <- posterior_binary(design$prior_t, data_t$n, data_t$y,
                             external = design$external_t)
  p <- region_prob_binary(design$criterion, post_t$a, post_t$b, post_c$a,
                          post_c$b)
  list(prob_go = p$go, prob_nogo = p$nogo,
       posterior_t = c(a = post_t$a, b = post_t$b),
       posterior_c = c(a = post_c$a, b = post_c$b))
}
