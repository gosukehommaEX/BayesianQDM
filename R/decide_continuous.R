#' Posterior distributions and region probabilities for observed continuous
#' data
#'
#' @param design An object of class \code{"qdm_design_continuous"}.
#' @param data_t,data_c Observed data of the treatment and control groups;
#'   \code{data_c} is \code{NULL} for a design without a concurrent control
#'   group.
#'
#' @return A list with elements \code{prob_go}, \code{prob_nogo},
#'   \code{posterior_t} and \code{posterior_c} (named vectors of the
#'   parameters \code{mu}, \code{kappa}, \code{nu} and \code{sigma}; for a
#'   fixed control, \code{posterior_c} is the fixed value named
#'   \code{value}).
#' @noRd
decide_continuous <- function(design, data_t, data_c) {
  if (!inherits(data_t, "qdm_data_continuous")) {
    stop("'data_t' must be created by data_continuous().", call. = FALSE)
  }
  post_t <- posterior_continuous(design$prior_t, data_t$n, data_t$mean,
                                 data_t$sd, external = design$external_t)
  if (design$control == "concurrent") {
    if (!inherits(data_c, "qdm_data_continuous")) {
      stop(paste("'data_c' must be created by data_continuous() for a",
                 "design with a concurrent control group."), call. = FALSE)
    }
    post_c <- control_posterior_continuous(design, post_t, data_c$n,
                                           data_c$mean, data_c$sd)
  } else {
    if (!is.null(data_c)) {
      stop(paste("'data_c' must be NULL for a design without a concurrent",
                 "control group."), call. = FALSE)
    }
    post_c <- control_posterior_continuous(design, post_t)
  }
  if (!is.null(post_c)) {
    check_df_mm(post_t$nu, post_c$nu, design$method)
  }
  p <- region_prob_continuous(design$criterion, post_t, post_c,
                              design$method, design$nMC,
                              fixed = design$fixed_c)
  params <- function(post) {
    if (is.null(post)) {
      return(c(value = design$fixed_c))
    }
    c(mu = post$mu, kappa = post$kappa, nu = post$nu,
      sigma = sqrt(post$sigma2))
  }
  list(prob_go = p$go, prob_nogo = p$nogo, posterior_t = params(post_t),
       posterior_c = params(post_c))
}
