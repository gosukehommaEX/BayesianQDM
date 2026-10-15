#' Region probabilities for observed data of two binary endpoints
#'
#' @param design An object of class \code{"qdm_design_binary2"}.
#' @param data_t,data_c Observed data of the treatment and control groups;
#'   \code{data_c} is \code{NULL} for a design without a concurrent control
#'   group.
#'
#' @return A list with elements \code{prob_go}, \code{prob_nogo} (one number,
#'   or one per endpoint for the marginal rules), \code{prob_regions} (the
#'   probabilities of all regions, \code{R1} to \code{R9} or \code{R1} to
#'   \code{R4}), \code{posterior_t} and \code{posterior_c} (named vectors of
#'   Dirichlet parameters; for a fixed control, \code{posterior_c} is the
#'   fixed values named \code{value1} and \code{value2}). All probabilities
#'   are computed from the same Monte Carlo draws.
#' @noRd
decide_binary2 <- function(design, data_t, data_c) {
  if (!inherits(data_t, "qdm_data_binary2")) {
    stop("'data_t' must be created by data_binary2().", call. = FALSE)
  }
  alpha_t <- posterior_binary2(design$prior_t, data_t$counts,
                               external = design$external_t)
  if (design$control == "concurrent") {
    if (!inherits(data_c, "qdm_data_binary2")) {
      stop(paste("'data_c' must be created by data_binary2() for a design",
                 "with a concurrent control group."), call. = FALSE)
    }
    alpha_c <- posterior_binary2(design$prior_c, data_c$counts,
                                 external = design$external_c)
  } else {
    if (!is.null(data_c)) {
      stop(paste("'data_c' must be NULL for a design without a concurrent",
                 "control group."), call. = FALSE)
    }
    alpha_c <- if (design$control == "fixed") {
      NULL
    } else {
      posterior_binary2(design$prior_c, c(0, 0, 0, 0),
                        external = design$external_c,
                        hypothetical = design$hypothetical_c)
    }
  }
  crit <- design$criterion
  n_regions <- if (crit$type == "posterior") 9L else 4L
  sets <- region_sets(crit)
  p <- region_prob_binary2(crit, alpha_t, alpha_c, 1L, 1L,
                           c(as.list(seq_len(n_regions)), sets$go,
                             sets$nogo), design$nMC, fixed = design$fixed_c)
  n_go <- length(sets$go)
  list(prob_go = p[1L, n_regions + seq_len(n_go)],
       prob_nogo = p[1L, n_regions + n_go + seq_len(n_go)],
       prob_regions = stats::setNames(p[1L, seq_len(n_regions)],
                                      paste0("R", seq_len(n_regions))),
       posterior_t = alpha_t[1L, ],
       posterior_c = if (is.null(alpha_c)) {
         c(value1 = design$fixed_c[1L], value2 = design$fixed_c[2L])
       } else {
         alpha_c[1L, ]
       })
}
