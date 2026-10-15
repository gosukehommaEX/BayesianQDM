#' Design with two continuous endpoints
#'
#' Specifies a proof-of-concept study with two continuous endpoints: the
#' sample sizes, the prior distributions, the external data, the form of the
#' control group, the decision criterion for two endpoints and the method
#' used to compute the probabilities. The design is used by
#' \code{\link{qdm_decide}}, \code{\link{qdm_oc}} and
#' \code{\link{qdm_calibrate}}.
#'
#' The responses of a subject on the two endpoints have a bivariate normal
#' distribution with mean vector \eqn{\mu_k} and covariance matrix
#' \eqn{\Sigma_k} in group \eqn{k} (Yamaguchi et al., submitted, Appendix
#' A.2). The prior distribution of \eqn{(\mu_k, \Sigma_k)} is
#' Normal-Inverse-Wishart (\code{\link{prior_niw}}) or vague
#' (\code{\link{prior_vague}}, the density \eqn{|\Sigma_k|^{-3/2}}). External
#' data enter through a power prior with weight \eqn{a_0}, which has the
#' conjugate form, and the posterior parameters are updated with the external
#' data and then with the current data. The vague prior corresponds to
#' \eqn{\kappa_0 = 0}, \eqn{\nu_0 = -1} and \eqn{\Lambda_0 = 0}, so the
#' posterior distribution of the mean is
#' \eqn{t_{n - 2}(\bar y, S / (n (n - 2)))} without external data, where
#' \eqn{S} is the matrix of sums of squares and cross-products, \eqn{n - 1}
#' times the sample covariance matrix (Yamaguchi et al., submitted,
#' equation A.9).
#'
#' The treatment effects are the differences of the means,
#' \eqn{\theta_k = \mu_{tk} - \mu_{ck}}. The mean of each group has a
#' bivariate t distribution with \eqn{\nu_n - 1} degrees of freedom,
#' location \eqn{\mu_n} and scale matrix
#' \eqn{\Lambda_n / (\kappa_n (\nu_n - 1))}. For a predictive criterion,
#' the mean of \eqn{m} future subjects has the same distribution with the
#' scale matrix \eqn{(1 / \kappa_n + 1 / m) \Lambda_n / (\nu_n - 1)},
#' because the future subjects share the same mean and covariance matrix
#' (version 0.1.0 and Yamaguchi et al. used the scale matrix
#' \eqn{(1 + \kappa_n) \Lambda_n / (\kappa_n (\nu_n - 1) m)}, which treats
#' them as independent draws from the predictive distribution of one
#' subject). A hypothetical control (\code{hypothetical(mean = c(m1, m2),
#' r = r)}) gives the control mean the distribution of the treatment mean
#' with location \code{mean} and the scale matrix multiplied by \code{r}
#' (Yamaguchi et al., Section A.2.6, with the hyperparameters of the
#' treatment group). With external control data only (\code{n_c = 0} and
#' \code{external_c}) and the vague prior, \eqn{a_0 n_e} must be greater
#' than 2. A fixed control (\code{n_c = 0} and \code{fixed_c}, the known
#' control means of the two endpoints) shifts the distribution of the
#' treatment mean, so the effect has a bivariate t distribution and the
#' probabilities of the regions are computed exactly from its distribution
#' function; \code{method}, \code{nMC} and \code{m_c} are then not used,
#' and \code{prior_c} must be \code{NULL}.
#'
#' The probabilities of the regions are computed by one of two methods:
#' \describe{
#'   \item{\code{"MM"}}{moment-matching approximation: the effect is
#'     approximated by a bivariate t distribution with the same mean,
#'     covariance matrix and radial fourth moment (Yamaguchi, Homma, Takeda
#'     and Maruo, submitted, Theorem 3), whose rectangle
#'     probabilities are computed by numerical integration for real degrees
#'     of freedom. It requires more than 4 degrees of freedom for both
#'     groups and is the fastest method.}
#'   \item{\code{"MC"}}{Monte Carlo simulation with \code{nMC} draws from
#'     each distribution; the probabilities of all regions are estimated
#'     from the same draws. The computing time is proportional to
#'     \code{nMC}, so \code{"MM"} is much faster for operating
#'     characteristics, which need the probabilities for many simulated
#'     studies.}
#' }
#'
#' The operating characteristics are computed by simulating the sample means
#' and the sums of squares and cross-products of \code{nsim} studies (see
#' \code{\link{qdm_oc}}).
#'
#' @param n_t Number of subjects in the treatment group (a whole number of at
#'   least 3).
#' @param n_c Number of subjects in the concurrent control group: 0 or a
#'   whole number of at least 3.
#' @param prior_t,prior_c Prior distributions of the treatment and control
#'   groups, created by \code{\link{prior_niw}} or \code{\link{prior_vague}}.
#'   For a hypothetical or fixed control, \code{prior_c} must be \code{NULL}
#'   (the default).
#' @param criterion Decision criterion for two endpoints, created by
#'   \code{\link{criterion_posterior}} or \code{\link{criterion_predictive}}.
#' @param external_t,external_c External data of the treatment and control
#'   groups, created by \code{\link{borrow}} with
#'   \code{\link{data_continuous2}}, or \code{NULL} (no external data).
#' @param hypothetical_c Hypothetical control, created by
#'   \code{\link{hypothetical}} with \code{mean} of length 2 and \code{r},
#'   or \code{NULL}. It can be used only when \code{n_c = 0}.
#' @param fixed_c Fixed control means of endpoints 1 and 2, a numeric vector
#'   of length 2, or \code{NULL}. It can be used only when \code{n_c = 0}.
#' @param method Method used to compute the probabilities: \code{"MM"}
#'   (default) or \code{"MC"}.
#' @param nMC Number of Monte Carlo draws for \code{method = "MC"}.
#'
#' @return An object of class \code{c("qdm_design_continuous2",
#'   "qdm_design")}, a list with the arguments and the element
#'   \code{control} (\code{"concurrent"}, \code{"external"},
#'   \code{"hypothetical"} or \code{"fixed"}).
#'
#' @references
#' Yamaguchi, Y., Hu, H., Huang, Y., Han, C., Homma, G. and Wang, A.
#' (submitted). Bivariate Bayesian quantitative decision-making for two
#' endpoints in rare disease drug development.
#'
#' Yamaguchi, Y., Homma, G., Takeda, K. and Maruo, K. (submitted).
#' Moment-matching approximation for t-differences in clinical trials with
#' continuous endpoints.
#'
#' @family designs
#' @examples
#' # Case study 2 of Yamaguchi et al. (submitted)
#' L0 <- matrix(c(80^2, 15, 15, 6^2), 2, 2)
#' des <- design_continuous2(
#'   n_t = 30, n_c = 30,
#'   prior_t = prior_niw(mu0 = c(0, 1), kappa0 = 0.1, nu0 = 4, Lambda0 = L0),
#'   prior_c = prior_niw(mu0 = c(-10, -1), kappa0 = 0.1, nu0 = 4,
#'                       Lambda0 = L0),
#'   criterion = criterion_posterior(theta_tv = c(10, 2),
#'                                   theta_mav = c(5, 1),
#'                                   gamma_go = 0.296, gamma_nogo = 0.244)
#' )
#' des
#' qdm_decide(des,
#'            data_t = data_continuous2(30, c(5, 2),
#'                                      matrix(c(6400, 190, 190, 36), 2, 2)),
#'            data_c = data_continuous2(30, c(-10, -1),
#'                                      matrix(c(5600, 150, 150, 42), 2, 2)))
#' @export
design_continuous2 <- function(n_t, n_c, prior_t, prior_c = NULL, criterion,
                               external_t = NULL, external_c = NULL,
                               hypothetical_c = NULL, fixed_c = NULL,
                               method = c("MM", "MC"), nMC = 10000L) {
  check_numeric(n_t, "n_t", lower = 3, integer = TRUE)
  check_numeric(n_c, "n_c", lower = 0, integer = TRUE)
  if (round(n_c) %in% c(1, 2)) {
    stop("'n_c' must be 0 or a whole number of at least 3.", call. = FALSE)
  }
  check_prior_continuous2(prior_t, "prior_t")
  if (!inherits(criterion, "qdm_criterion") || criterion$n_endpoints != 2L) {
    stop(paste("'criterion' must be created by criterion_posterior() or",
               "criterion_predictive() for two endpoints."), call. = FALSE)
  }
  check_external_continuous2(external_t, "external_t")
  check_external_continuous2(external_c, "external_c")
  if (!is.null(hypothetical_c) &&
      (!inherits(hypothetical_c, "qdm_hypothetical") ||
       hypothetical_c$type != "location" ||
       length(hypothetical_c$mean) != 2L)) {
    stop(paste("'hypothetical_c' must be NULL or created by hypothetical()",
               "with 'mean' of length 2 and 'r'."), call. = FALSE)
  }
  fixed_c <- check_fixed_c(fixed_c, 2L, probability = FALSE)
  control <- control_form(n_c, external_c, hypothetical_c, fixed_c)
  check_prior_control(prior_c, control, check_prior_continuous2,
                      unused = c("hypothetical", "fixed"))
  method <- match.arg(method)
  check_numeric(nMC, "nMC", lower = 1, upper = .Machine$integer.max,
                integer = TRUE)
  des <- structure(list(endpoint = "continuous2",
                        n_t = as.integer(round(n_t)),
                        n_c = as.integer(round(n_c)), prior_t = prior_t,
                        prior_c = prior_c, criterion = criterion,
                        external_t = external_t, external_c = external_c,
                        hypothetical_c = hypothetical_c, fixed_c = fixed_c,
                        control = control, method = method,
                        nMC = as.integer(round(nMC))),
                   class = c("qdm_design_continuous2", "qdm_design"))
  if (control == "fixed") {
    # The probabilities are computed exactly, whatever the method
    return(des)
  }
  # Degrees of freedom of the distributions of the group means for the
  # sample sizes of the design (they do not depend on the data values)
  post_t <- posterior_continuous2(prior_t, des$n_t, external = external_t)
  post_c <- control_posterior_continuous2(des, post_t, des$n_c)
  if (post_c$nu <= 1) {
    stop(paste("With prior_vague() for the control group and external",
               "control data only, 'a0' times the number of external",
               "subjects must be greater than 2, so that the distribution",
               "of the control mean is proper."), call. = FALSE)
  }
  check_df_mm(post_t$nu - 1, post_c$nu - 1, method, "method = \"MC\"")
  des
}
