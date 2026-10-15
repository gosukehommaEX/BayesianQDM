#' Region probabilities for observed data of two continuous endpoints
#'
#' @param design An object of class \code{"qdm_design_continuous2"}.
#' @param data_t,data_c Observed data of the treatment and control groups;
#'   \code{data_c} is \code{NULL} for a design without a concurrent control
#'   group.
#'
#' @return A list with elements \code{prob_go}, \code{prob_nogo} (one number,
#'   or one per endpoint for the marginal rules), \code{prob_regions} (the
#'   probabilities of all regions, \code{R1} to \code{R9} or \code{R1} to
#'   \code{R4}), \code{posterior_t} and \code{posterior_c} (named vectors of
#'   the parameters \code{mu1}, \code{mu2}, \code{kappa}, \code{nu},
#'   \code{L11}, \code{L12} and \code{L22}; for a fixed control,
#'   \code{posterior_c} is the fixed values named \code{value1} and
#'   \code{value2}). With \code{method = "MC"}, all probabilities are
#'   estimated from the same draws.
#' @noRd
decide_continuous2 <- function(design, data_t, data_c) {
  if (!inherits(data_t, "qdm_data_continuous2")) {
    stop("'data_t' must be created by data_continuous2().", call. = FALSE)
  }
  post_t <- posterior_continuous2(design$prior_t, data_t$n, data_t$mean,
                                  data_t$cov, external = design$external_t)
  if (design$control == "concurrent") {
    if (!inherits(data_c, "qdm_data_continuous2")) {
      stop(paste("'data_c' must be created by data_continuous2() for a",
                 "design with a concurrent control group."), call. = FALSE)
    }
    post_c <- control_posterior_continuous2(design, post_t, data_c$n,
                                            data_c$mean, data_c$cov)
  } else {
    if (!is.null(data_c)) {
      stop(paste("'data_c' must be NULL for a design without a concurrent",
                 "control group."), call. = FALSE)
    }
    post_c <- control_posterior_continuous2(design, post_t)
  }
  crit <- design$criterion
  dist_c <- NULL
  if (!is.null(post_c)) {
    check_df_mm(post_t$nu - 1, post_c$nu - 1, design$method,
                "method = \"MC\"")
    dist_c <- mean_dist_continuous2(post_c, crit$type, crit$m_c)
  }
  n_regions <- if (crit$type == "posterior") 9L else 4L
  sets <- region_sets(crit)
  dist_t <- mean_dist_continuous2(post_t, crit$type, crit$m_t)
  p <- region_prob_continuous2(crit, dist_t, dist_c,
                               c(as.list(seq_len(n_regions)), sets$go,
                                 sets$nogo), design$method, design$nMC,
                               fixed = design$fixed_c)
  n_go <- length(sets$go)
  params <- function(post) {
    if (is.null(post)) {
      return(c(value1 = design$fixed_c[1L], value2 = design$fixed_c[2L]))
    }
    c(mu1 = post$mu[1L, 1L], mu2 = post$mu[1L, 2L], kappa = post$kappa,
      nu = post$nu, post$Lambda[1L, ])
  }
  list(prob_go = p[1L, n_regions + seq_len(n_go)],
       prob_nogo = p[1L, n_regions + n_go + seq_len(n_go)],
       prob_regions = stats::setNames(p[1L, seq_len(n_regions)],
                                      paste0("R", seq_len(n_regions))),
       posterior_t = params(post_t), posterior_c = params(post_c))
}
