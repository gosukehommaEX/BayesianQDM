#' Decision boundary of a design
#'
#' Computes the decision of a design (Go, NoGo, Gray or Miss) for every point
#' of a grid of possible observed data of the treatment group, so that the
#' regions of the data that lead to each decision can be plotted
#' (\code{\link{plot.qdm_boundary}}). Unlike the operating characteristics,
#' the result does not depend on the true parameters: it shows which data
#' lead to which decision.
#'
#' The grid depends on the endpoint of the design:
#' \describe{
#'   \item{One binary endpoint}{With a concurrent control group, all numbers
#'     of responders \eqn{y_t = 0, \ldots, n_t} and \eqn{y_c = 0, \ldots, n_c}
#'     (\code{n_t} must be a single number; \eqn{n_c} is that of the design).
#'     Without one, all \eqn{(n, y_t)} with \eqn{n} in \code{n_t} and
#'     \eqn{y_t = 0, \ldots, n}, which shows how the boundary changes with the
#'     sample size.}
#'   \item{One continuous endpoint}{All combinations of \code{mean_t},
#'     \code{sd_t} and \code{n_t}, of which at most one of \code{sd_t} and
#'     \code{n_t} may have more than one value (the horizontal axis). With a
#'     concurrent control group, \code{data_c} gives the observed control
#'     data and the vertical axis is the observed difference of the means.}
#'   \item{Two continuous endpoints}{All combinations of the observed
#'     treatment means of the two endpoints, \code{mean_t[[1]]} and
#'     \code{mean_t[[2]]}, with the sample covariance matrix \code{cov_t} and
#'     \code{n_t} subjects (a single number). With a concurrent control
#'     group, \code{data_c} gives the observed control data and the axes are
#'     the observed differences of the means.}
#'   \item{One time-to-event endpoint}{All combinations of the numbers of
#'     events \code{events_t} and the total observed times
#'     \code{exposure_t} of the treatment group. With a concurrent control
#'     group, \code{data_c} gives the observed control data.}
#' }
#' Two binary endpoints are not supported, because their decision depends on
#' the four cell counts and cannot be shown on two axes.
#'
#' The decision at each point is made as in \code{\link{qdm_decide}}.
#' Probabilities computed by Monte Carlo simulation (\code{method = "MC"},
#' or posterior predictive probabilities of a time-to-event endpoint) make
#' the boundary irregular; deterministic methods give smooth boundaries.
#'
#' @param design A design created by \code{\link{design_binary}},
#'   \code{\link{design_continuous}}, \code{\link{design_continuous2}} or
#'   \code{\link{design_tte}}.
#' @param n_t Numbers of subjects in the treatment group (whole numbers); by
#'   default that of the design. Several values are allowed for one binary
#'   endpoint without a concurrent control group and for one continuous
#'   endpoint. For a time-to-event endpoint, \code{n_t} only sets the default
#'   numbers of events.
#' @param mean_t Observed treatment means: a numeric vector for one
#'   continuous endpoint, or a list of two numeric vectors (endpoints 1 and
#'   2) for two continuous endpoints.
#' @param sd_t Observed treatment standard deviations (one continuous
#'   endpoint; positive).
#' @param cov_t Observed treatment covariance matrix (two continuous
#'   endpoints; a 2 x 2 positive definite matrix).
#' @param events_t Numbers of events of the treatment group (time-to-event
#'   endpoint; whole numbers); by default \code{0:n_t}.
#' @param exposure_t Total observed times of the treatment group
#'   (time-to-event endpoint; positive).
#' @param data_c Observed data of the concurrent control group, created by
#'   \code{\link{data_continuous}}, \code{\link{data_continuous2}} or
#'   \code{\link{data_tte}}; required for those endpoints with a concurrent
#'   control group and \code{NULL} otherwise.
#' @param seed Random seed for probabilities computed by Monte Carlo
#'   simulation, or \code{NULL} to use the current random number stream. If
#'   given, the previous state of the random number generator is restored on
#'   exit.
#'
#' @return An object of class \code{c("qdm_boundary", "data.frame")}: the
#'   grid (columns named as the arguments, and \code{y_t} and \code{y_c} for
#'   binary data or \code{diff}, \code{diff1} and \code{diff2} for
#'   differences of means), the columns \code{prob_go} and \code{prob_nogo}
#'   (for two endpoints with a marginal rule, \code{prob_go1},
#'   \code{prob_go2}, \code{prob_nogo1} and \code{prob_nogo2}) and the
#'   factor \code{decision}. The attributes \code{design} and \code{axes}
#'   record the design and the columns and labels of the axes of the plot.
#'
#' @family operations
#' @examples
#' # Single-arm study with a fixed control response rate of 0.2
#' des <- design_binary(
#'   n_t = 20, n_c = 0, prior_t = prior_beta(1, 1),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.05,
#'                                   gamma_go = 0.8, gamma_nogo = 0.8),
#'   fixed_c = 0.2
#' )
#' b <- qdm_boundary(des, n_t = 10:30)
#' plot(b)
#'
#' # Continuous endpoint with a concurrent control group
#' des2 <- design_continuous(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_vague(), prior_c = prior_vague(),
#'   criterion = criterion_posterior(theta_tv = 30, theta_mav = 15,
#'                                   gamma_go = 0.6, gamma_nogo = 0.7),
#'   method = "MM"
#' )
#' b2 <- qdm_boundary(des2, mean_t = seq(-20, 60, by = 2),
#'                    sd_t = seq(10, 40, by = 2),
#'                    data_c = data_continuous(12, 0, 25))
#' plot(b2)
#' @export
qdm_boundary <- function(design, n_t = NULL, mean_t = NULL, sd_t = NULL,
                         cov_t = NULL, events_t = NULL, exposure_t = NULL,
                         data_c = NULL, seed = NULL) {
  check_design(design)
  if (design$endpoint == "binary2") {
    stop(paste("qdm_boundary() does not support two binary endpoints,",
               "whose decision depends on the four cell counts."),
         call. = FALSE)
  }
  if (is.null(n_t)) {
    n_t <- design$n_t
  }
  used <- switch(design$endpoint,
                 binary = "n_t",
                 continuous = c("n_t", "mean_t", "sd_t", "data_c"),
                 continuous2 = c("n_t", "mean_t", "cov_t", "data_c"),
                 tte = c("n_t", "events_t", "exposure_t", "data_c"))
  given <- c(mean_t = !is.null(mean_t), sd_t = !is.null(sd_t),
             cov_t = !is.null(cov_t), events_t = !is.null(events_t),
             exposure_t = !is.null(exposure_t), data_c = !is.null(data_c))
  unused <- setdiff(names(given)[given], used)
  if (length(unused) > 0L) {
    stop(sprintf("%s %s not used for this endpoint.",
                 paste0("'", unused, "'", collapse = ", "),
                 if (length(unused) > 1L) "are" else "is"), call. = FALSE)
  }
  concurrent <- design$control == "concurrent"
  if ("data_c" %in% used && concurrent && is.null(data_c)) {
    stop("'data_c' must be given for a design with a concurrent control group.",
         call. = FALSE)
  }
  if (!concurrent && !is.null(data_c)) {
    stop(paste("'data_c' must be NULL for a design without a concurrent",
               "control group."), call. = FALSE)
  }
  res <- with_seed(seed, switch(
    design$endpoint,
    binary = boundary_binary(design, n_t),
    continuous = boundary_continuous(design, n_t, mean_t, sd_t, data_c),
    continuous2 = boundary_continuous2(design, n_t, mean_t, cov_t, data_c),
    tte = boundary_tte(design, n_t, events_t, exposure_t, data_c)
  ))
  crit <- design$criterion
  met_go <- meets_criterion(res$p_go, crit$gamma_go, crit$rule)
  met_nogo <- meets_criterion(res$p_nogo, crit$gamma_nogo, crit$rule)
  out <- res$grid
  if (is.list(res$p_go)) {
    out$prob_go1 <- res$p_go[[1L]]
    out$prob_go2 <- res$p_go[[2L]]
    out$prob_nogo1 <- res$p_nogo[[1L]]
    out$prob_nogo2 <- res$p_nogo[[2L]]
  } else {
    out$prob_go <- res$p_go
    out$prob_nogo <- res$p_nogo
  }
  out$decision <- factor(decision_label(met_go, met_nogo),
                         levels = c("Go", "NoGo", "Gray", "Miss"))
  rownames(out) <- NULL
  structure(out, class = c("qdm_boundary", "data.frame"), design = design,
            axes = res$axes)
}
