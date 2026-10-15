#' Decision boundary of a design
#'
#' Finds the observed data that lead to each decision (Go, NoGo, Gray, or
#' Miss) of a design: for example, how many responders out of \eqn{n} are
#' needed for a Go decision, or how large the observed difference in means
#' must be. Unlike the operating characteristics, the result does not depend
#' on the true parameters. It is a table of decision intervals of an
#' observed quantity (the numbers of responders or events, or an observed
#' mean), one set of intervals for each value of a second quantity (the
#' condition), together with the corresponding observed response rates,
#' milestone survival probabilities and treatment effects. The intervals are
#' printed as tables (\code{\link{print.qdm_boundary}}) and drawn as the
#' regions of the decisions with the boundary values written on them
#' (\code{\link{plot.qdm_boundary}}).
#'
#' The observed quantity and the conditions depend on the endpoint:
#' \describe{
#'   \item{One binary endpoint}{The number of responders of the treatment
#'     group \eqn{y_t = 0, \ldots, n_t}. With a concurrent control group, one
#'     set of intervals for each number of control responders
#'     \eqn{y_c = 0, \ldots, n_c} (\code{n_t} must be a single number;
#'     \eqn{n_c} is that of the design). Without one, one set for each
#'     sample size in \code{n_t}, which shows how the boundary changes with
#'     the sample size.}
#'   \item{One continuous endpoint}{The observed treatment mean, for each
#'     combination of \code{sd_t} and \code{n_t} (at most one of them may
#'     have more than one value). With a concurrent control group,
#'     \code{data_c} gives the observed control data. The boundary values
#'     are located by bisection between the values of a grid, \code{mean_t}
#'     or by default a grid of 201 values that covers the thresholds of the
#'     criterion.}
#'   \item{Two continuous endpoints}{The observed treatment mean of endpoint
#'     2, located by bisection on the grid \code{mean_t[[2]]}, for each
#'     observed treatment mean of endpoint 1 in \code{mean_t[[1]]}, with the
#'     sample covariance matrix \code{cov_t} and \code{n_t} subjects (a
#'     single number). With a concurrent control group, \code{data_c} gives
#'     the observed control data.}
#'   \item{One time-to-event endpoint}{The number of events of the treatment
#'     group (\code{events_t}), for each total observed time in
#'     \code{exposure_t}. With a concurrent control group, \code{data_c}
#'     gives the observed control data.}
#' }
#' Two binary endpoints are not supported, because their decision depends on
#' the four cell counts and cannot be summarized by one observed quantity.
#'
#' The decision for each value is made as in \code{\link{qdm_decide}}. For an
#' observed mean, the first and last intervals extend to minus and plus
#' infinity; changes of the decision beyond the grid (whose ends are
#' extended by two values far beyond them) or two changes of the same
#' criterion between two neighboring grid values are not detected. With a
#' conjugate (Normal-Inverse-chi-squared) prior, the posterior scale
#' increases with the distance between the observed mean and the prior
#' mean, so the probability of the Go region need not increase with the
#' observed mean; the intervals then show where it decreases again.
#' Probabilities computed by Monte Carlo simulation (\code{method = "MC"},
#' or posterior predictive probabilities of a time-to-event endpoint) make
#' the boundary values approximate.
#'
#' @param design A design created by \code{\link{design_binary}},
#'   \code{\link{design_continuous}}, \code{\link{design_continuous2}} or
#'   \code{\link{design_tte}}.
#' @param n_t Numbers of subjects in the treatment group (whole numbers); by
#'   default that of the design. Several values are allowed for one binary
#'   endpoint without a concurrent control group and for one continuous
#'   endpoint. For a time-to-event endpoint, \code{n_t} only sets the default
#'   numbers of events.
#' @param mean_t For one continuous endpoint, an optional numeric vector:
#'   the grid of observed treatment means on which the boundary values are
#'   searched. For two continuous endpoints, a list of two numeric vectors:
#'   the observed treatment means of endpoint 1 (the conditions) and the
#'   grid of observed treatment means of endpoint 2.
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
#' @return An object of class \code{c("qdm_boundary", "data.frame")} with one
#'   row per decision interval and the columns
#'   \describe{
#'     \item{Conditions}{\code{n_t}, \code{n_c} and \code{y_c} (binary;
#'       only \code{n_t} without a concurrent control group),
#'       \code{n_t} and \code{sd_t} (continuous), \code{n_t},
#'       \code{mean_t1} and \code{effect1} (two continuous endpoints; the
#'       observed treatment mean of endpoint 1 and its difference from the
#'       control) or \code{exposure_t} (time-to-event).}
#'     \item{decision}{The decision (a factor with levels \code{"Go"},
#'       \code{"NoGo"}, \code{"Gray"} and \code{"Miss"}).}
#'     \item{from, to}{The ends of the interval of the observed quantity: the
#'       number of treatment responders or events (both ends included), or
#'       the observed treatment mean (of endpoint 2).}
#'     \item{rate_from, rate_to}{Binary endpoint: the observed response
#'       rates \eqn{y_t / n_t} at the ends.}
#'     \item{surv_from, surv_to}{Time-to-event endpoint evaluated by the
#'       milestone survival: the observed milestone survival probabilities
#'       \eqn{\exp(-t_0 d_t / E_t)} at the ends.}
#'     \item{effect_from, effect_to}{The observed treatment effect at the
#'       ends: the observed rate, mean or milestone survival minus that of
#'       the concurrent control group, or minus the fixed control value or
#'       the mean (location) of the distribution of the control group; for
#'       the hazard ratio, the observed hazard \eqn{d_t / E_t} divided by the
#'       observed control hazard (\code{NA} if the control group has no
#'       events), the fixed control hazard or the mean of the distribution
#'       of the control hazard. The values correspond to
#'       \code{from} and \code{to}, so \code{effect_from} is larger than
#'       \code{effect_to} when larger values of the observed quantity give
#'       smaller effects (for example, more events and the milestone
#'       survival).}
#'   }
#'   The attributes \code{design}, \code{axes}, \code{titles},
#'   \code{increasing}, \code{discrete} and \code{range} record the design
#'   and the information used by the print and plot methods.
#'
#' @family operations
#' @examples
#' # Single-arm study with a fixed control response rate of 0.2: responders
#' # needed out of 10 to 30 subjects
#' des <- design_binary(
#'   n_t = 20, n_c = 0, prior_t = prior_beta(1, 1),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.05,
#'                                   gamma_go = 0.8, gamma_nogo = 0.8),
#'   fixed_c = 0.2
#' )
#' b <- qdm_boundary(des, n_t = c(10, 20, 30))
#' b
#' plot(b)
#'
#' # Continuous endpoint with a concurrent control group: observed
#' # differences in means needed for each observed standard deviation
#' des2 <- design_continuous(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_vague(), prior_c = prior_vague(),
#'   criterion = criterion_posterior(theta_tv = 30, theta_mav = 15,
#'                                   gamma_go = 0.6, gamma_nogo = 0.7),
#'   method = "MM"
#' )
#' b2 <- qdm_boundary(des2, sd_t = c(20, 25, 30),
#'                    data_c = data_continuous(12, 0, 25))
#' b2
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
  spec <- switch(
    design$endpoint,
    binary = boundary_binary(design, n_t),
    continuous = boundary_continuous(design, n_t, mean_t, sd_t, data_c),
    continuous2 = boundary_continuous2(design, n_t, mean_t, cov_t, data_c),
    tte = boundary_tte(design, n_t, events_t, exposure_t, data_c)
  )
  iv <- with_seed(seed, boundary_intervals(spec, design$criterion))
  out <- spec$cond[iv$row, , drop = FALSE]
  out$decision <- factor(iv$decision, levels = c("Go", "NoGo", "Gray", "Miss"))
  out$from <- iv$from
  out$to <- iv$to
  at_from <- spec$derived(iv$row, iv$from)
  at_to <- spec$derived(iv$row, iv$to)
  for (nm in names(at_from)) {
    out[[paste0(nm, "_from")]] <- at_from[[nm]]
    out[[paste0(nm, "_to")]] <- at_to[[nm]]
  }
  rownames(out) <- NULL
  structure(out, class = c("qdm_boundary", "data.frame"), design = design,
            axes = spec$axes, titles = spec$titles,
            increasing = spec$increasing, discrete = spec$discrete,
            range = spec$range)
}
