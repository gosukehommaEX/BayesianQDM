#' Bayesian Posterior or Posterior Predictive Probability for a Single
#' Continuous Endpoint
#'
#' Computes the Bayesian posterior probability or posterior predictive
#' probability for continuous-outcome clinical trials under a
#' Normal-Inverse-Chi-squared (or vague Jeffreys) conjugate model. The
#' function supports controlled, uncontrolled, and external designs, with
#' optional incorporation of external data through power priors. Vector
#' inputs for \code{bar_y_t}, \code{s_t}, \code{bar_y_c}, and \code{s_c}
#' are supported for efficient batch processing (e.g., across simulation
#' replicates in \code{\link{pbayesdecisionprob1cont}}).
#'
#' @param prob A character string specifying the probability type.
#'        Must be \code{'posterior'} or \code{'predictive'}.
#' @param design A character string specifying the trial design.
#'        Must be \code{'controlled'}, \code{'uncontrolled'}, or
#'        \code{'external'}.
#' @param prior A character string specifying the prior type.
#'        Must be \code{'vague'} (Jeffreys) or \code{'N-Inv-Chisq'}
#'        (Normal-Inverse-Chi-squared conjugate).
#' @param CalcMethod A character string specifying the calculation method.
#'        Must be \code{'NI'} (numerical integration), \code{'MC'}
#'        (Monte Carlo), or \code{'MM'} (Moment-Matching approximation).
#'        For large \code{nsim} or multi-scenario use, \code{'MM'} is
#'        strongly recommended; see Details.
#' @param theta0 A numeric scalar giving the threshold for the treatment
#'        effect (difference in means).
#' @param nMC A positive integer giving the number of Monte Carlo draws.
#'        Required when \code{CalcMethod = 'MC'}; otherwise \code{NULL}.
#' @param n_t A positive integer giving the number of patients in the
#'        treatment group in the proof-of-concept (PoC) trial.
#' @param n_c A positive integer giving the number of patients in the
#'        control group in the PoC trial. Required for
#'        \code{design = 'controlled'} or \code{'external'}; set to
#'        \code{NULL} for \code{design = 'uncontrolled'}.
#' @param m_t A positive integer giving the number of patients in the
#'        treatment group for the future trial. Required when
#'        \code{prob = 'predictive'}; otherwise set to \code{NULL}.
#' @param m_c A positive integer giving the number of patients in the
#'        control group for the future trial. Required when
#'        \code{prob = 'predictive'}; otherwise set to \code{NULL}.
#' @param kappa0_t A positive numeric scalar giving the prior precision
#'        parameter for the treatment group. Required when
#'        \code{prior = 'N-Inv-Chisq'}; otherwise \code{NULL}.
#' @param kappa0_c A positive numeric scalar giving the prior precision
#'        parameter for the control group. Required when
#'        \code{prior = 'N-Inv-Chisq'} and \code{design = 'controlled'};
#'        otherwise \code{NULL}.
#' @param nu0_t A positive numeric scalar giving the prior degrees of freedom
#'        for the treatment group. Required when \code{prior = 'N-Inv-Chisq'};
#'        otherwise \code{NULL}.
#' @param nu0_c A positive numeric scalar giving the prior degrees of freedom
#'        for the control group. Required when \code{prior = 'N-Inv-Chisq'}
#'        and \code{design = 'controlled'}; otherwise \code{NULL}.
#' @param mu0_t A numeric scalar giving the prior mean for the treatment
#'        group. Required when \code{prior = 'N-Inv-Chisq'}; otherwise
#'        \code{NULL}.
#' @param mu0_c A numeric scalar giving the prior mean for the control group
#'        (\code{design = 'controlled'} with \code{prior = 'N-Inv-Chisq'})
#'        or the hypothetical control mean (\code{design = 'uncontrolled'}).
#'        Required for uncontrolled design; otherwise \code{NULL}.
#' @param sigma0_t A positive numeric scalar giving the prior scale for the
#'        treatment group. Required when \code{prior = 'N-Inv-Chisq'};
#'        otherwise \code{NULL}.
#' @param sigma0_c A positive numeric scalar giving the prior scale for the
#'        control group. Required when \code{prior = 'N-Inv-Chisq'} and
#'        \code{design = 'controlled'}; otherwise \code{NULL}.
#' @param bar_y_t A numeric scalar or vector giving the sample mean for the
#'        treatment group. When a vector of length \code{nsim} is supplied,
#'        all posterior parameters are computed simultaneously for all
#'        replicates.
#' @param bar_y_c A numeric scalar or vector giving the sample mean for the
#'        control group. Required for \code{design = 'controlled'} or
#'        \code{'external'}; set to \code{NULL} for uncontrolled design.
#' @param s_t A positive numeric scalar or vector giving the sample standard
#'        deviation for the treatment group.
#' @param s_c A positive numeric scalar or vector giving the sample standard
#'        deviation for the control group. Required for
#'        \code{design = 'controlled'} or \code{'external'}; otherwise
#'        \code{NULL}.
#' @param r A positive numeric scalar giving the variance scaling factor for
#'        the hypothetical control. Required for \code{design = 'uncontrolled'};
#'        otherwise \code{NULL}. The hypothetical control scale is
#'        \eqn{\mathrm{sd.control} = \sqrt{r} \cdot \mathrm{sd.treatment}}.
#' @param ne_t A positive integer giving the number of patients in the
#'        treatment group of the external data set. Required when
#'        \code{design = 'external'} and external treatment data are
#'        available; otherwise set to \code{NULL}.
#' @param ne_c A positive integer giving the number of patients in the
#'        control group of the external data set. Required when
#'        \code{design = 'external'} and external control data are available;
#'        otherwise set to \code{NULL}.
#' @param alpha0e_t A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for the external treatment data. Required when external
#'        treatment data are used; otherwise set to \code{NULL}.
#' @param alpha0e_c A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for the external control data. Required when external
#'        control data are used; otherwise set to \code{NULL}.
#' @param bar_ye_t A numeric scalar giving the external treatment group
#'        sample mean. Required when \code{ne_t} is provided; otherwise
#'        \code{NULL}.
#' @param bar_ye_c A numeric scalar giving the external control group sample
#'        mean. Required when \code{ne_c} is provided; otherwise \code{NULL}.
#' @param se_t A positive numeric scalar giving the external treatment group
#'        sample SD. Required when \code{ne_t} is provided; otherwise
#'        \code{NULL}.
#' @param se_c A positive numeric scalar giving the external control group
#'        sample SD. Required when \code{ne_c} is provided; otherwise
#'        \code{NULL}.
#' @param lower.tail A logical scalar; if \code{TRUE} (default), the function
#'        returns \eqn{P(\mathrm{effect} \le \theta_0)}, otherwise
#'        \eqn{P(\mathrm{effect} > \theta_0)}.
#'
#' @return A numeric scalar or vector in \code{[0, 1]}. When the input data
#'         parameters (\code{bar_y_t}, etc.) are vectors of length \eqn{n},
#'         a vector of length \eqn{n} is returned.
#'
#' @details
#' Under the Normal-Inverse-Chi-squared (N-Inv-ChiSq) conjugate model, the
#' marginal posterior/predictive distribution of the group mean follows a
#' non-standardised t-distribution, parameterised by a location
#' \eqn{\mu_j}, a scale \eqn{\sigma_j}, and degrees of freedom \eqn{\nu_j}
#' that depend on the design and probability type. The final probability is
#' computed by one of three methods:
#' \itemize{
#'   \item \code{NI}: exact numerical integration via \code{\link{ptdiff_NI}}.
#'   \item \code{MC}: Monte Carlo simulation via \code{\link{ptdiff_MC}}.
#'   \item \code{MM}: Moment-Matching approximation via
#'         \code{\link{ptdiff_MM}}.
#' }
#'
#' \strong{Performance note}: \code{CalcMethod = 'NI'} calls
#' \code{integrate()} once per element of the input vector, which is
#' prohibitively slow for large vectors (e.g., \code{nsim x n_scenarios}).
#' \code{CalcMethod = 'MC'} creates an \code{nMC x n} matrix internally,
#' consuming \eqn{O(nMC \cdot n)} memory. For large-scale simulation use
#' \code{CalcMethod = 'MM'}, which is fully vectorised and exact in the
#' normal limit.
#'
#' @examples
#' # Example 1: Controlled design - posterior probability with N-Inv-Chisq
#' pbayespostpred1cont(
#'   prob = 'posterior', design = 'controlled', prior = 'N-Inv-Chisq',
#'   CalcMethod = 'NI', theta0 = 2, n_t = 12, n_c = 12,
#'   kappa0_t = 5, kappa0_c = 5, nu0_t = 5, nu0_c = 5,
#'   mu0_t = 5, mu0_c = 5, sigma0_t = sqrt(5), sigma0_c = sqrt(5),
#'   bar_y_t = 2, bar_y_c = 0, s_t = 1, s_c = 1,
#'   m_t = NULL, m_c = NULL, r = NULL,
#'   ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
#'   bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
#'   lower.tail = FALSE
#' )
#'
#' # Example 2: Uncontrolled design - posterior probability with vague prior
#' pbayespostpred1cont(
#'   prob = 'posterior', design = 'uncontrolled', prior = 'vague',
#'   CalcMethod = 'MM', theta0 = 1.5, n_t = 15, n_c = NULL,
#'   bar_y_t = 3.5, bar_y_c = NULL, s_t = 1.2, s_c = NULL,
#'   mu0_c = 1.5, r = 1.2,
#'   kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
#'   mu0_t = NULL, sigma0_t = NULL, sigma0_c = NULL,
#'   m_t = NULL, m_c = NULL,
#'   ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
#'   bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
#'   lower.tail = FALSE
#' )
#'
#' # Example 3: External design - posterior probability
#' pbayespostpred1cont(
#'   prob = 'posterior', design = 'external', prior = 'N-Inv-Chisq',
#'   CalcMethod = 'MM', theta0 = 2, n_t = 12, n_c = 12,
#'   kappa0_t = 5, kappa0_c = 5, nu0_t = 5, nu0_c = 5,
#'   mu0_t = 5, mu0_c = 5, sigma0_t = sqrt(5), sigma0_c = sqrt(5),
#'   bar_y_t = 2.5, bar_y_c = 1, s_t = 1.1, s_c = 0.9,
#'   m_t = NULL, m_c = NULL, r = NULL,
#'   ne_t = 10, ne_c = 10, alpha0e_t = 0.5, alpha0e_c = 0.5,
#'   bar_ye_t = 2, bar_ye_c = 0.5, se_t = 1, se_c = 0.8,
#'   lower.tail = FALSE
#' )
#'
#' # Example 4: Controlled design - posterior predictive probability
#' pbayespostpred1cont(
#'   prob = 'predictive', design = 'controlled', prior = 'N-Inv-Chisq',
#'   CalcMethod = 'MM', theta0 = 1, n_t = 12, n_c = 12,
#'   kappa0_t = 5, kappa0_c = 5, nu0_t = 5, nu0_c = 5,
#'   mu0_t = 5, mu0_c = 5, sigma0_t = sqrt(5), sigma0_c = sqrt(5),
#'   bar_y_t = 2, bar_y_c = 0, s_t = 1, s_c = 1,
#'   m_t = 20, m_c = 20, r = NULL,
#'   ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
#'   bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
#'   lower.tail = FALSE
#' )
#'
#' # Example 5: Uncontrolled design - posterior predictive probability
#' pbayespostpred1cont(
#'   prob = 'predictive', design = 'uncontrolled', prior = 'vague',
#'   CalcMethod = 'MM', theta0 = 1.5, n_t = 15, n_c = NULL,
#'   bar_y_t = 3.5, bar_y_c = NULL, s_t = 1.2, s_c = NULL,
#'   mu0_c = 1.5, r = 1.2,
#'   kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
#'   mu0_t = NULL, sigma0_t = NULL, sigma0_c = NULL,
#'   m_t = 20, m_c = 20,
#'   ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
#'   bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
#'   lower.tail = FALSE
#' )
#'
#' # Example 6: External design - posterior predictive probability
#' pbayespostpred1cont(
#'   prob = 'predictive', design = 'external', prior = 'N-Inv-Chisq',
#'   CalcMethod = 'MM', theta0 = 1, n_t = 12, n_c = 12,
#'   kappa0_t = 5, kappa0_c = 5, nu0_t = 5, nu0_c = 5,
#'   mu0_t = 5, mu0_c = 5, sigma0_t = sqrt(5), sigma0_c = sqrt(5),
#'   bar_y_t = 2.5, bar_y_c = 1, s_t = 1.1, s_c = 0.9,
#'   m_t = 20, m_c = 20, r = NULL,
#'   ne_t = 10, ne_c = 10, alpha0e_t = 0.5, alpha0e_c = 0.5,
#'   bar_ye_t = 2, bar_ye_c = 0.5, se_t = 1, se_c = 0.8,
#'   lower.tail = FALSE
#' )
#'
#' @export
pbayespostpred1cont <- function(prob = "posterior", design = "controlled",
                                prior = "vague", CalcMethod = "NI",
                                theta0, nMC = NULL, n_t, n_c = NULL,
                                m_t = NULL, m_c = NULL,
                                kappa0_t = NULL, kappa0_c = NULL,
                                nu0_t = NULL, nu0_c = NULL,
                                mu0_t = NULL, mu0_c = NULL,
                                sigma0_t = NULL, sigma0_c = NULL,
                                bar_y_t, bar_y_c = NULL, s_t, s_c = NULL,
                                r = NULL,
                                ne_t = NULL, ne_c = NULL,
                                alpha0e_t = NULL, alpha0e_c = NULL,
                                bar_ye_t = NULL, bar_ye_c = NULL,
                                se_t = NULL, se_c = NULL,
                                lower.tail = TRUE) {

  # --- Input validation ---
  if (!is.character(prob) || length(prob) != 1L ||
      !prob %in% c("posterior", "predictive")) {
    stop("'prob' must be either 'posterior' or 'predictive'")
  }

  if (!is.character(design) || length(design) != 1L ||
      !design %in% c("controlled", "uncontrolled", "external")) {
    stop("'design' must be 'controlled', 'uncontrolled', or 'external'")
  }

  if (!is.character(prior) || length(prior) != 1L ||
      !prior %in% c("vague", "N-Inv-Chisq")) {
    stop("'prior' must be either 'vague' or 'N-Inv-Chisq'")
  }

  if (!is.character(CalcMethod) || length(CalcMethod) != 1L ||
      !CalcMethod %in% c("NI", "MC", "MM")) {
    stop("'CalcMethod' must be 'NI', 'MC', or 'MM'")
  }

  if (!is.numeric(theta0) || length(theta0) != 1L || is.na(theta0)) {
    stop("'theta0' must be a single numeric value")
  }

  if (!is.numeric(n_t) || length(n_t) != 1L || is.na(n_t) ||
      n_t != floor(n_t) || n_t < 1L) {
    stop("'n_t' must be a single positive integer")
  }

  if (design != "uncontrolled") {
    if (is.null(n_c) || !is.numeric(n_c) || length(n_c) != 1L || is.na(n_c) ||
        n_c != floor(n_c) || n_c < 1L) {
      stop("'n_c' must be a single positive integer for controlled or external design")
    }
  }

  if (!is.logical(lower.tail) || length(lower.tail) != 1L || is.na(lower.tail)) {
    stop("'lower.tail' must be a single logical value (TRUE or FALSE)")
  }

  # Validate vector data arguments
  if (!is.numeric(bar_y_t) || length(bar_y_t) < 1L || any(is.na(bar_y_t))) {
    stop("'bar_y_t' must be a numeric scalar or vector with no missing values")
  }

  if (!is.numeric(s_t) || length(s_t) < 1L || any(is.na(s_t)) || any(s_t <= 0)) {
    stop("'s_t' must be a positive numeric scalar or vector with no missing values")
  }

  # Validate CalcMethod-specific parameters
  if (CalcMethod == "MC") {
    if (is.null(nMC)) {
      stop("'nMC' must be specified when CalcMethod = 'MC'")
    }
    if (!is.numeric(nMC) || length(nMC) != 1L || is.na(nMC) ||
        nMC != floor(nMC) || nMC < 1L) {
      stop("'nMC' must be a single positive integer")
    }
  }

  # Validate prob-specific parameters
  if (prob == "predictive") {
    if (is.null(m_t) || is.null(m_c)) {
      stop("'m_t' and 'm_c' must be non-NULL when prob = 'predictive'")
    }
    if (!is.numeric(m_t) || length(m_t) != 1L || is.na(m_t) ||
        m_t != floor(m_t) || m_t < 1L) {
      stop("'m_t' must be a single positive integer")
    }
    if (!is.numeric(m_c) || length(m_c) != 1L || is.na(m_c) ||
        m_c != floor(m_c) || m_c < 1L) {
      stop("'m_c' must be a single positive integer")
    }
  }

  # Validate prior-specific parameters for N-Inv-Chisq
  if (prior == "N-Inv-Chisq") {
    for (nm in c("kappa0_t", "nu0_t", "mu0_t", "sigma0_t")) {
      val <- get(nm)
      if (is.null(val)) {
        stop(paste0("'", nm, "' must be non-NULL when prior = 'N-Inv-Chisq'"))
      }
    }
    if (!is.numeric(kappa0_t) || length(kappa0_t) != 1L || is.na(kappa0_t) || kappa0_t <= 0) {
      stop("'kappa0_t' must be a single positive numeric value")
    }
    if (!is.numeric(nu0_t) || length(nu0_t) != 1L || is.na(nu0_t) || nu0_t <= 0) {
      stop("'nu0_t' must be a single positive numeric value")
    }
    if (!is.numeric(mu0_t) || length(mu0_t) != 1L || is.na(mu0_t)) {
      stop("'mu0_t' must be a single numeric value")
    }
    if (!is.numeric(sigma0_t) || length(sigma0_t) != 1L || is.na(sigma0_t) || sigma0_t <= 0) {
      stop("'sigma0_t' must be a single positive numeric value")
    }
    if (design %in% c("controlled", "external")) {
      for (nm in c("kappa0_c", "nu0_c", "mu0_c", "sigma0_c")) {
        val <- get(nm)
        if (is.null(val)) {
          stop(paste0("'", nm, "' must be non-NULL when prior = 'N-Inv-Chisq' and design = '",
                      design, "'"))
        }
      }
      if (!is.numeric(kappa0_c) || length(kappa0_c) != 1L || is.na(kappa0_c) || kappa0_c <= 0) {
        stop("'kappa0_c' must be a single positive numeric value")
      }
      if (!is.numeric(nu0_c) || length(nu0_c) != 1L || is.na(nu0_c) || nu0_c <= 0) {
        stop("'nu0_c' must be a single positive numeric value")
      }
      if (!is.numeric(sigma0_c) || length(sigma0_c) != 1L || is.na(sigma0_c) || sigma0_c <= 0) {
        stop("'sigma0_c' must be a single positive numeric value")
      }
    }
  }

  # Validate design-specific parameters
  if (design == "uncontrolled") {
    if (is.null(mu0_c)) {
      stop("'mu0_c' (hypothetical control mean) must be non-NULL when design = 'uncontrolled'")
    }
    if (!is.numeric(mu0_c) || length(mu0_c) != 1L || is.na(mu0_c)) {
      stop("'mu0_c' must be a single numeric value")
    }
    if (is.null(r)) {
      stop("'r' (variance scaling factor) must be non-NULL when design = 'uncontrolled'")
    }
    if (!is.numeric(r) || length(r) != 1L || is.na(r) || r <= 0) {
      stop("'r' must be a single positive numeric value")
    }
  }

  if (design == "external") {
    if (is.null(ne_t) && is.null(ne_c)) {
      stop("For design = 'external', at least one of 'ne_t' or 'ne_c' must be non-NULL")
    }
    if (!is.null(ne_t)) {
      if (is.null(alpha0e_t) || is.null(bar_ye_t) || is.null(se_t)) {
        stop("'alpha0e_t', 'bar_ye_t', and 'se_t' must all be non-NULL when 'ne_t' is provided")
      }
      if (!is.numeric(ne_t) || length(ne_t) != 1L || is.na(ne_t) ||
          ne_t != floor(ne_t) || ne_t < 1L) {
        stop("'ne_t' must be a single positive integer")
      }
      if (!is.numeric(alpha0e_t) || length(alpha0e_t) != 1L || is.na(alpha0e_t) ||
          alpha0e_t <= 0 || alpha0e_t > 1) {
        stop("'alpha0e_t' must be a single numeric value in (0, 1]")
      }
      if (!is.numeric(bar_ye_t) || length(bar_ye_t) != 1L || is.na(bar_ye_t)) {
        stop("'bar_ye_t' must be a single numeric value")
      }
      if (!is.numeric(se_t) || length(se_t) != 1L || is.na(se_t) || se_t <= 0) {
        stop("'se_t' must be a single positive numeric value")
      }
    }
    if (!is.null(ne_c)) {
      if (is.null(alpha0e_c) || is.null(bar_ye_c) || is.null(se_c)) {
        stop("'alpha0e_c', 'bar_ye_c', and 'se_c' must all be non-NULL when 'ne_c' is provided")
      }
      if (!is.numeric(ne_c) || length(ne_c) != 1L || is.na(ne_c) ||
          ne_c != floor(ne_c) || ne_c < 1L) {
        stop("'ne_c' must be a single positive integer")
      }
      if (!is.numeric(alpha0e_c) || length(alpha0e_c) != 1L || is.na(alpha0e_c) ||
          alpha0e_c <= 0 || alpha0e_c > 1) {
        stop("'alpha0e_c' must be a single numeric value in (0, 1]")
      }
      if (!is.numeric(bar_ye_c) || length(bar_ye_c) != 1L || is.na(bar_ye_c)) {
        stop("'bar_ye_c' must be a single numeric value")
      }
      if (!is.numeric(se_c) || length(se_c) != 1L || is.na(se_c) || se_c <= 0) {
        stop("'se_c' must be a single positive numeric value")
      }
    }
    if (!is.null(bar_y_c) && (!is.numeric(bar_y_c) || any(is.na(bar_y_c)))) {
      stop("'bar_y_c' must be a numeric scalar or vector with no missing values")
    }
    if (!is.null(s_c) && (!is.numeric(s_c) || any(is.na(s_c)) || any(s_c <= 0))) {
      stop("'s_c' must be a positive numeric scalar or vector with no missing values")
    }
  }

  if (design == "controlled") {
    if (is.null(bar_y_c) || is.null(s_c)) {
      stop("'bar_y_c' and 's_c' must be non-NULL when design = 'controlled'")
    }
    if (!is.numeric(bar_y_c) || any(is.na(bar_y_c))) {
      stop("'bar_y_c' must be a numeric scalar or vector with no missing values")
    }
    if (!is.numeric(s_c) || any(is.na(s_c)) || any(s_c <= 0)) {
      stop("'s_c' must be a positive numeric scalar or vector with no missing values")
    }
  }

  # ---------------------------------------------------------------------------
  # Compute posterior t-distribution parameters (vectorised over bar_y_t, s_t,
  # bar_y_c, s_c; all arithmetic operates element-wise on vectors).
  # ---------------------------------------------------------------------------

  if (design == 'external') {

    if (prior == 'N-Inv-Chisq') {

      # --- External design with N-Inv-Chisq prior ---
      # Power prior for normal data is conjugate with the N-Inv-Chisq prior.
      # The posterior hyperparameters incorporate both the N-Inv-Chisq prior
      # and the external data via the power prior weight.

      # Group 1 (Treatment): incorporate external data if available
      if (!is.null(ne_t) && !is.null(alpha0e_t)) {
        # Intermediate: N-Inv-Chisq prior updated with power-weighted external data
        kappa_e_t  <- kappa0_t + alpha0e_t * ne_t
        nu_e_t     <- nu0_t + alpha0e_t * ne_t
        mu_e_t     <- (kappa0_t * mu0_t + alpha0e_t * ne_t * bar_ye_t) / kappa_e_t
        var_e_t    <- (nu0_t * sigma0_t ^ 2 + alpha0e_t * (ne_t - 1) * se_t ^ 2 +
                         alpha0e_t * ne_t * kappa0_t * (bar_ye_t - mu0_t) ^ 2 / kappa_e_t) / nu_e_t
        # Final update with PoC data
        kappa_n_t  <- kappa_e_t + n_t
        nu_t     <- nu_e_t + n_t
        mu_t     <- (kappa_e_t * mu_e_t + n_t * bar_y_t) / kappa_n_t
        var_n_t    <- (nu_e_t * var_e_t + (n_t - 1) * s_t ^ 2 +
                         n_t * kappa_e_t * (mu_e_t - bar_y_t) ^ 2 / kappa_n_t) / nu_t
      } else {
        # No external treatment data: N-Inv-Chisq prior updated with PoC data only
        kappa_n_t  <- kappa0_t + n_t
        nu_t     <- nu0_t + n_t
        mu_t     <- (kappa0_t * mu0_t + n_t * bar_y_t) / kappa_n_t
        var_n_t    <- (nu0_t * sigma0_t ^ 2 + (n_t - 1) * s_t ^ 2 +
                         n_t * kappa0_t * (mu0_t - bar_y_t) ^ 2 / kappa_n_t) / nu_t
      }

      # Group 2 (Control): incorporate external data if available
      if (!is.null(ne_c) && !is.null(alpha0e_c)) {
        # Intermediate: N-Inv-Chisq prior updated with power-weighted external data
        kappa_e_c  <- kappa0_c + alpha0e_c * ne_c
        nu_e_c     <- nu0_c + alpha0e_c * ne_c
        mu_e_c     <- (kappa0_c * mu0_c + alpha0e_c * ne_c * bar_ye_c) / kappa_e_c
        var_e_c    <- (nu0_c * sigma0_c ^ 2 + alpha0e_c * (ne_c - 1) * se_c ^ 2 +
                         alpha0e_c * ne_c * kappa0_c * (bar_ye_c - mu0_c) ^ 2 / kappa_e_c) / nu_e_c
        # Final update with PoC data
        kappa_n_c  <- kappa_e_c + n_c
        nu_c     <- nu_e_c + n_c
        mu_c     <- (kappa_e_c * mu_e_c + n_c * bar_y_c) / kappa_n_c
        var_n_c    <- (nu_e_c * var_e_c + (n_c - 1) * s_c ^ 2 +
                         n_c * kappa_e_c * (mu_e_c - bar_y_c) ^ 2 / kappa_n_c) / nu_c
      } else {
        # No external control data: N-Inv-Chisq prior updated with PoC data only
        kappa_n_c  <- kappa0_c + n_c
        nu_c     <- nu0_c + n_c
        mu_c     <- (kappa0_c * mu0_c + n_c * bar_y_c) / kappa_n_c
        var_n_c    <- (nu0_c * sigma0_c ^ 2 + (n_c - 1) * s_c ^ 2 +
                         n_c * kappa0_c * (mu0_c - bar_y_c) ^ 2 / kappa_n_c) / nu_c
      }

      # Scale parameters depend on probability type
      if (prob == 'posterior') {
        sd_t <- sqrt(var_n_t / kappa_n_t)
        sd_c <- sqrt(var_n_c / kappa_n_c)
      } else {
        sd_t <- sqrt((1 + 1 / kappa_n_t) * var_n_t / m_t)
        sd_c <- sqrt((1 + 1 / kappa_n_c) * var_n_c / m_c)
      }

    } else {

      # --- External design with vague prior ---
      # Group 1 (Treatment): apply power prior if external data available
      if (!is.null(ne_t) && !is.null(alpha0e_t)) {
        mu_t          <- (alpha0e_t * ne_t * bar_ye_t + n_t * bar_y_t) / (alpha0e_t * ne_t + n_t)
        kappa_star_n_t  <- alpha0e_t * ne_t + n_t
        nu_t          <- alpha0e_t * ne_t + n_t - 1
        sigma2_star_n_t <- (alpha0e_t * (ne_t - 1) * se_t ^ 2 + (n_t - 1) * s_t ^ 2 +
                              (alpha0e_t * ne_t * n_t * (bar_ye_t - bar_y_t) ^ 2) /
                              (alpha0e_t * ne_t + n_t)) / (alpha0e_t * ne_t + n_t)
      } else {
        # No external treatment data: use vague prior
        mu_t          <- bar_y_t
        kappa_star_n_t  <- n_t
        nu_t          <- n_t - 1
        sigma2_star_n_t <- s_t ^ 2
      }

      # Group 2 (Control): apply power prior if external data available
      if (!is.null(ne_c) && !is.null(alpha0e_c)) {
        mu_c          <- (alpha0e_c * ne_c * bar_ye_c + n_c * bar_y_c) / (alpha0e_c * ne_c + n_c)
        kappa_star_n_c  <- alpha0e_c * ne_c + n_c
        nu_c          <- alpha0e_c * ne_c + n_c - 1
        sigma2_star_n_c <- (alpha0e_c * (ne_c - 1) * se_c ^ 2 + (n_c - 1) * s_c ^ 2 +
                              (alpha0e_c * ne_c * n_c * (bar_ye_c - bar_y_c) ^ 2) /
                              (alpha0e_c * ne_c + n_c)) / (alpha0e_c * ne_c + n_c)
      } else {
        # No external control data: use vague prior
        mu_c          <- bar_y_c
        kappa_star_n_c  <- n_c
        nu_c          <- n_c - 1
        sigma2_star_n_c <- s_c ^ 2
      }

      # Scale parameters depend on probability type
      if (prob == 'posterior') {
        sd_t <- sqrt(sigma2_star_n_t / kappa_star_n_t)
        sd_c <- sqrt(sigma2_star_n_c / kappa_star_n_c)
      } else {
        sd_t <- sqrt((1 + 1 / kappa_star_n_t) * sigma2_star_n_t / m_t)
        sd_c <- sqrt((1 + 1 / kappa_star_n_c) * sigma2_star_n_c / m_c)
      }
    }

  } else {

    # --- Controlled or uncontrolled designs (vague or N-Inv-Chisq prior) ---
    if (prior == 'N-Inv-Chisq') {

      # Updated precision parameters
      kappa_n_t <- kappa0_t + n_t
      kappa_n_c <- kappa0_c + n_c

      # Updated degrees of freedom
      nu_t <- nu0_t + n_t
      if (design == 'controlled') {
        nu_c <- nu0_c + n_c
      } else {
        nu_c <- nu_t
      }

      # Posterior means
      mu_t <- (kappa0_t * mu0_t + n_t * bar_y_t) / kappa_n_t
      if (design == 'controlled') {
        mu_c <- (kappa0_c * mu0_c + n_c * bar_y_c) / kappa_n_c
      } else {
        mu_c <- mu0_c
      }

      # Posterior variance for the treatment group (vectorised over bar_y_t, s_t)
      var_n_t <- (nu0_t * sigma0_t ^ 2 + (n_t - 1) * s_t ^ 2 +
                    n_t * kappa0_t * (mu0_t - bar_y_t) ^ 2 / kappa_n_t) / nu_t

      # Posterior variance for the control group
      if (design == 'controlled') {
        var_n_c <- (nu0_c * sigma0_c ^ 2 + (n_c - 1) * s_c ^ 2 +
                      n_c * kappa0_c * (mu0_c - bar_y_c) ^ 2 / kappa_n_c) / nu_c
      } else {
        var_n_c <- NULL
      }

      # Scale parameters depend on probability type
      if (prob == 'posterior') {
        sd_t <- sqrt(var_n_t / kappa_n_t)
        if (design == 'controlled') {
          sd_c <- sqrt(var_n_c / kappa_n_c)
        } else {
          sd_c <- sqrt(r) * sd_t
        }
      } else {
        sd_t <- sqrt((1 + 1 / kappa_n_t) * var_n_t / m_t)
        if (design == 'controlled') {
          sd_c <- sqrt((1 + 1 / kappa_n_c) * var_n_c / m_c)
        } else {
          sd_c <- sqrt(r) * sd_t
        }
      }

    } else {
      # Vague (Jeffreys) prior

      # Degrees of freedom
      nu_t <- n_t - 1
      if (design == 'controlled') {
        nu_c <- n_c - 1
      } else {
        nu_c <- nu_t
      }

      # Posterior means (vectorised over bar_y_t, bar_y_c)
      mu_t <- bar_y_t
      if (design == 'controlled') {
        mu_c <- bar_y_c
      } else {
        mu_c <- mu0_c
      }

      # Scale parameters (vectorised over s_t, s_c)
      if (prob == 'posterior') {
        sd_t <- sqrt(s_t ^ 2 / n_t)
        if (design == 'controlled') {
          sd_c <- sqrt(s_c ^ 2 / n_c)
        } else {
          sd_c <- sqrt(r) * sd_t
        }
      } else {
        sd_t <- sqrt((1 + 1 / n_t) * s_t ^ 2 / m_t)
        if (design == 'controlled') {
          sd_c <- sqrt((1 + 1 / n_c) * s_c ^ 2 / m_c)
        } else {
          sd_c <- sqrt(r) * sd_t
        }
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Compute probability using the specified method.
  # mu_t, mu_c, sd_t, sd_c may be vectors of length nsim when called from
  # pbayesdecisionprob1cont; each method receives the full vector and returns a
  # vector of the same length (no outer loop required).
  # ---------------------------------------------------------------------------
  if (CalcMethod == 'NI') {
    results <- ptdiff_NI(theta0, mu_t, mu_c, sd_t, sd_c, nu_t, nu_c, lower.tail)
  } else if (CalcMethod == 'MC') {
    results <- ptdiff_MC(nMC, theta0, mu_t, mu_c, sd_t, sd_c, nu_t, nu_c, lower.tail)
  } else {
    results <- ptdiff_MM(theta0, mu_t, mu_c, sd_t, sd_c, nu_t, nu_c, lower.tail)
  }

  return(results)
}
