#' Bayesian Posterior or Posterior Predictive Probability for a Clinical Trial
#' with a Single Continuous Endpoint
#'
#' Computes the Bayesian posterior probability or posterior predictive probability
#' for continuous-outcome clinical trials under a Normal-Inverse-Chi-squared
#' (or vague Jeffreys) conjugate model. The function supports controlled,
#' uncontrolled, and external-control designs, with optional incorporation of
#' external data through power priors. Vector inputs for \code{bar.y1}, \code{s1},
#' \code{bar.y2}, and \code{s2} are supported for efficient batch processing
#' (e.g., across simulation replicates in \code{\link{pbayesdecisionprob1cont}}).
#'
#' @param prob A character string specifying the probability type.
#'        Must be \code{'posterior'} or \code{'predictive'}.
#' @param design A character string specifying the trial design.
#'        Must be \code{'controlled'}, \code{'uncontrolled'}, or \code{'external'}.
#' @param prior A character string specifying the prior type.
#'        Must be \code{'vague'} (Jeffreys) or \code{'N-Inv-Chisq'}
#'        (Normal-Inverse-Chi-squared conjugate).
#' @param CalcMethod A character string specifying the calculation method.
#'        Must be \code{'NI'} (numerical integration), \code{'MC'} (Monte Carlo),
#'        or \code{'MM'} (Moment-Matching approximation).
#'        For large \code{nsim} or multi-scenario use, \code{'MM'} is strongly
#'        recommended; see Details.
#' @param theta0 A numeric scalar giving the threshold for the treatment effect
#'        (difference in means).
#' @param nMC A positive integer giving the number of Monte Carlo draws.
#'        Required when \code{CalcMethod = 'MC'}; otherwise \code{NULL}.
#' @param n1 A positive integer giving the number of patients in group 1
#'        (treatment) in the proof-of-concept (PoC) trial.
#' @param n2 A positive integer giving the number of patients in group 2
#'        in the PoC trial. For \code{design = 'uncontrolled'}, set equal to
#'        \code{n1} (not used in calculations but required for consistency).
#' @param m1 A positive integer giving the future sample size for group 1.
#'        Required when \code{prob = 'predictive'}; otherwise \code{NULL}.
#' @param m2 A positive integer giving the future sample size for group 2.
#'        Required when \code{prob = 'predictive'}; otherwise \code{NULL}.
#' @param kappa01 A positive numeric scalar giving the prior precision parameter
#'        for group 1. Required when \code{prior = 'N-Inv-Chisq'}; otherwise
#'        \code{NULL}.
#' @param kappa02 A positive numeric scalar giving the prior precision parameter
#'        for group 2. Required when \code{prior = 'N-Inv-Chisq'} and
#'        \code{design = 'controlled'}; otherwise \code{NULL}.
#' @param nu01 A positive numeric scalar giving the prior degrees of freedom
#'        for group 1. Required when \code{prior = 'N-Inv-Chisq'}; otherwise
#'        \code{NULL}.
#' @param nu02 A positive numeric scalar giving the prior degrees of freedom
#'        for group 2. Required when \code{prior = 'N-Inv-Chisq'} and
#'        \code{design = 'controlled'}; otherwise \code{NULL}.
#' @param mu01 A numeric scalar giving the prior mean for group 1.
#'        Required when \code{prior = 'N-Inv-Chisq'}; otherwise \code{NULL}.
#' @param mu02 A numeric scalar giving the prior mean for group 2
#'        (\code{design = 'controlled'} with \code{prior = 'N-Inv-Chisq'}) or the
#'        hypothetical control mean (\code{design = 'uncontrolled'}).
#'        Required for uncontrolled design; otherwise \code{NULL}.
#' @param sigma01 A positive numeric scalar giving the prior scale for group 1.
#'        Required when \code{prior = 'N-Inv-Chisq'}; otherwise \code{NULL}.
#' @param sigma02 A positive numeric scalar giving the prior scale for group 2.
#'        Required when \code{prior = 'N-Inv-Chisq'} and
#'        \code{design = 'controlled'}; otherwise \code{NULL}.
#' @param bar.y1 A numeric scalar or vector giving the sample mean for group 1.
#'        When a vector of length \code{nsim} is supplied, all posterior parameters
#'        are computed simultaneously for all replicates.
#' @param bar.y2 A numeric scalar or vector giving the sample mean for group 2.
#'        Required for \code{design = 'controlled'} or \code{'external'}; set to
#'        \code{NULL} for uncontrolled design.
#' @param s1 A positive numeric scalar or vector giving the sample standard
#'        deviation for group 1.
#' @param s2 A positive numeric scalar or vector giving the sample standard
#'        deviation for group 2. Required for \code{design = 'controlled'} or
#'        \code{'external'}; otherwise \code{NULL}.
#' @param r A positive numeric scalar giving the variance scaling factor for
#'        the hypothetical control. Required for \code{design = 'uncontrolled'};
#'        otherwise \code{NULL}. The hypothetical control scale is
#'        \eqn{\mathrm{sd.control} = \sqrt{r} \cdot \mathrm{sd.treatment}}.
#' @param ne1 A positive integer giving the external group 1 sample size.
#'        Required for \code{design = 'external'} with external treatment data;
#'        otherwise \code{NULL}.
#' @param ne2 A positive integer giving the external group 2 sample size.
#'        Required for \code{design = 'external'} with external control data;
#'        otherwise \code{NULL}.
#' @param alpha01 A numeric scalar in \code{(0, 1]} giving the power prior weight
#'        for group 1. Required when \code{ne1} is provided; otherwise \code{NULL}.
#' @param alpha02 A numeric scalar in \code{(0, 1]} giving the power prior weight
#'        for group 2. Required when \code{ne2} is provided; otherwise \code{NULL}.
#' @param bar.ye1 A numeric scalar giving the external group 1 sample mean.
#'        Required when \code{ne1} is provided; otherwise \code{NULL}.
#' @param bar.ye2 A numeric scalar giving the external group 2 sample mean.
#'        Required when \code{ne2} is provided; otherwise \code{NULL}.
#' @param se1 A positive numeric scalar giving the external group 1 sample SD.
#'        Required when \code{ne1} is provided; otherwise \code{NULL}.
#' @param se2 A positive numeric scalar giving the external group 2 sample SD.
#'        Required when \code{ne2} is provided; otherwise \code{NULL}.
#' @param lower.tail A logical scalar; if \code{TRUE} (default), returns
#'        \eqn{P(\mathrm{effect} \le \theta_0)}, otherwise
#'        \eqn{P(\mathrm{effect} > \theta_0)}.
#'
#' @return A numeric scalar or vector in \code{[0, 1]}.  When the input data
#'         parameters (\code{bar.y1}, etc.) are vectors of length \eqn{n},
#'         a vector of length \eqn{n} is returned.
#'
#' @details
#' Under the Normal-Inverse-Chi-squared (N-Inv-ChiSq) conjugate model, the
#' marginal posterior/predictive distribution of the group mean follows a
#' non-standardised t-distribution, parameterised by a location \eqn{\mu_k},
#' a scale \eqn{\sigma_k}, and degrees of freedom \eqn{\nu_k} that depend on
#' the design and probability type. The final probability is computed by one of
#' three methods:
#' \itemize{
#'   \item \code{NI}: exact numerical integration via \code{\link{ptdiff_NI}}.
#'   \item \code{MC}: Monte Carlo simulation via \code{\link{ptdiff_MC}}.
#'   \item \code{MM}: Moment-Matching approximation via \code{\link{ptdiff_MM}}.
#' }
#'
#' \strong{Performance note}: \code{CalcMethod = 'NI'} calls \code{integrate()}
#' once per element of the input vector, which is prohibitively slow for large
#' vectors (e.g., \code{nsim} x \code{n_scenarios}). \code{CalcMethod = 'MC'}
#' creates an \code{nMC x n} matrix internally, consuming \eqn{O(nMC \cdot n)}
#' memory. For large-scale simulation use \code{CalcMethod = 'MM'}, which is
#' fully vectorised and exact in the normal limit.
#'
#' @examples
#' # Example 1: Controlled design - posterior probability with N-Inv-Chisq prior
#' pbayespostpred1cont(
#'   prob = 'posterior', design = 'controlled', prior = 'N-Inv-Chisq', CalcMethod = 'NI',
#'   theta0 = 2, n1 = 12, n2 = 12,
#'   kappa01 = 5, kappa02 = 5, nu01 = 5, nu02 = 5,
#'   mu01 = 5, mu02 = 5, sigma01 = sqrt(5), sigma02 = sqrt(5),
#'   bar.y1 = 2, bar.y2 = 0, s1 = 1, s2 = 1,
#'   m1 = NULL, m2 = NULL, r = NULL,
#'   ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
#'   bar.ye1 = NULL, bar.ye2 = NULL, se1 = NULL, se2 = NULL,
#'   lower.tail = FALSE
#' )
#'
#' # Example 2: Uncontrolled design with hypothetical control and vague prior
#' pbayespostpred1cont(
#'   prob = 'posterior', design = 'uncontrolled', prior = 'vague', CalcMethod = 'MM',
#'   theta0 = 1.5, n1 = 15, n2 = 15,
#'   bar.y1 = 3.5, bar.y2 = NULL, s1 = 1.2, s2 = NULL,
#'   mu02 = 1.5, r = 1.2,
#'   kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
#'   mu01 = NULL, sigma01 = NULL, sigma02 = NULL,
#'   m1 = NULL, m2 = NULL,
#'   ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
#'   bar.ye1 = NULL, bar.ye2 = NULL, se1 = NULL, se2 = NULL,
#'   lower.tail = FALSE
#' )
#'
#' # Example 3: Posterior predictive probability - controlled design
#' pbayespostpred1cont(
#'   prob = 'predictive', design = 'controlled', prior = 'N-Inv-Chisq', CalcMethod = 'MM',
#'   theta0 = 0.5, n1 = 15, n2 = 15, m1 = 100, m2 = 100,
#'   kappa01 = 3, kappa02 = 3, nu01 = 4, nu02 = 4,
#'   mu01 = 2, mu02 = 2, sigma01 = 2, sigma02 = 2,
#'   bar.y1 = 2.5, bar.y2 = 1.8, s1 = 1.8, s2 = 1.6,
#'   r = NULL, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
#'   bar.ye1 = NULL, bar.ye2 = NULL, se1 = NULL, se2 = NULL,
#'   lower.tail = FALSE
#' )
#'
#' # Example 4: External control design with power prior
#' pbayespostpred1cont(
#'   prob = 'posterior', design = 'external', prior = 'vague', CalcMethod = 'NI',
#'   theta0 = 1.5, n1 = 12, n2 = 12,
#'   bar.y1 = 4, bar.y2 = 2, s1 = 1.2, s2 = 1.1,
#'   ne2 = 20, alpha02 = 0.5, bar.ye2 = 1.8, se2 = 1.0,
#'   kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
#'   mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
#'   m1 = NULL, m2 = NULL, r = NULL,
#'   ne1 = NULL, alpha01 = NULL, bar.ye1 = NULL, se1 = NULL,
#'   lower.tail = FALSE
#' )
#'
#' # Example 5: Vectorised usage across simulation replicates
#' set.seed(1)
#' bar.y1_vec <- rnorm(100, mean = 3, sd = 1 / sqrt(12))
#' s1_vec     <- sqrt(rchisq(100, df = 11) * 1 / 11)
#' bar.y2_vec <- rnorm(100, mean = 1, sd = 1 / sqrt(12))
#' s2_vec     <- sqrt(rchisq(100, df = 11) * 1 / 11)
#' pbayespostpred1cont(
#'   prob = 'posterior', design = 'controlled', prior = 'vague', CalcMethod = 'MM',
#'   theta0 = 1.5, n1 = 12, n2 = 12,
#'   bar.y1 = bar.y1_vec, bar.y2 = bar.y2_vec, s1 = s1_vec, s2 = s2_vec,
#'   kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
#'   mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
#'   m1 = NULL, m2 = NULL, r = NULL,
#'   ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
#'   bar.ye1 = NULL, bar.ye2 = NULL, se1 = NULL, se2 = NULL,
#'   lower.tail = FALSE
#' )
#'
#' @export
pbayespostpred1cont <- function(prob = "posterior", design = "controlled",
                                prior = "vague", CalcMethod = "NI",
                                theta0, nMC = NULL, n1, n2,
                                m1 = NULL, m2 = NULL,
                                kappa01 = NULL, kappa02 = NULL,
                                nu01 = NULL, nu02 = NULL,
                                mu01 = NULL, mu02 = NULL,
                                sigma01 = NULL, sigma02 = NULL,
                                bar.y1, bar.y2 = NULL, s1, s2 = NULL,
                                r = NULL,
                                ne1 = NULL, ne2 = NULL,
                                alpha01 = NULL, alpha02 = NULL,
                                bar.ye1 = NULL, bar.ye2 = NULL,
                                se1 = NULL, se2 = NULL,
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

  if (!is.numeric(n1) || length(n1) != 1L || is.na(n1) ||
      n1 != floor(n1) || n1 < 1L) {
    stop("'n1' must be a single positive integer")
  }

  if (!is.numeric(n2) || length(n2) != 1L || is.na(n2) ||
      n2 != floor(n2) || n2 < 1L) {
    stop("'n2' must be a single positive integer")
  }

  if (!is.logical(lower.tail) || length(lower.tail) != 1L || is.na(lower.tail)) {
    stop("'lower.tail' must be a single logical value (TRUE or FALSE)")
  }

  # Validate vector data arguments
  if (!is.numeric(bar.y1) || length(bar.y1) < 1L || any(is.na(bar.y1))) {
    stop("'bar.y1' must be a numeric scalar or vector with no missing values")
  }

  if (!is.numeric(s1) || length(s1) < 1L || any(is.na(s1)) || any(s1 <= 0)) {
    stop("'s1' must be a positive numeric scalar or vector with no missing values")
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
    if (is.null(m1) || is.null(m2)) {
      stop("'m1' and 'm2' must be non-NULL when prob = 'predictive'")
    }
    if (!is.numeric(m1) || length(m1) != 1L || is.na(m1) ||
        m1 != floor(m1) || m1 < 1L) {
      stop("'m1' must be a single positive integer")
    }
    if (!is.numeric(m2) || length(m2) != 1L || is.na(m2) ||
        m2 != floor(m2) || m2 < 1L) {
      stop("'m2' must be a single positive integer")
    }
  }

  # Validate prior-specific parameters for N-Inv-Chisq
  if (prior == "N-Inv-Chisq") {
    for (nm in c("kappa01", "nu01", "mu01", "sigma01")) {
      val <- get(nm)
      if (is.null(val)) {
        stop(paste0("'", nm, "' must be non-NULL when prior = 'N-Inv-Chisq'"))
      }
    }
    if (!is.numeric(kappa01) || length(kappa01) != 1L || is.na(kappa01) || kappa01 <= 0) {
      stop("'kappa01' must be a single positive numeric value")
    }
    if (!is.numeric(nu01) || length(nu01) != 1L || is.na(nu01) || nu01 <= 0) {
      stop("'nu01' must be a single positive numeric value")
    }
    if (!is.numeric(mu01) || length(mu01) != 1L || is.na(mu01)) {
      stop("'mu01' must be a single numeric value")
    }
    if (!is.numeric(sigma01) || length(sigma01) != 1L || is.na(sigma01) || sigma01 <= 0) {
      stop("'sigma01' must be a single positive numeric value")
    }
    if (design == "controlled") {
      for (nm in c("kappa02", "nu02", "mu02", "sigma02")) {
        val <- get(nm)
        if (is.null(val)) {
          stop(paste0("'", nm, "' must be non-NULL when prior = 'N-Inv-Chisq' and design = 'controlled'"))
        }
      }
      if (!is.numeric(kappa02) || length(kappa02) != 1L || is.na(kappa02) || kappa02 <= 0) {
        stop("'kappa02' must be a single positive numeric value")
      }
      if (!is.numeric(nu02) || length(nu02) != 1L || is.na(nu02) || nu02 <= 0) {
        stop("'nu02' must be a single positive numeric value")
      }
      if (!is.numeric(sigma02) || length(sigma02) != 1L || is.na(sigma02) || sigma02 <= 0) {
        stop("'sigma02' must be a single positive numeric value")
      }
    }
  }

  # Validate design-specific parameters
  if (design == "uncontrolled") {
    if (is.null(mu02)) {
      stop("'mu02' (hypothetical control mean) must be non-NULL when design = 'uncontrolled'")
    }
    if (!is.numeric(mu02) || length(mu02) != 1L || is.na(mu02)) {
      stop("'mu02' must be a single numeric value")
    }
    if (is.null(r)) {
      stop("'r' (variance scaling factor) must be non-NULL when design = 'uncontrolled'")
    }
    if (!is.numeric(r) || length(r) != 1L || is.na(r) || r <= 0) {
      stop("'r' must be a single positive numeric value")
    }
  }

  if (design == "external") {
    if (is.null(ne1) && is.null(ne2)) {
      stop("For design = 'external', at least one of 'ne1' or 'ne2' must be non-NULL")
    }
    if (!is.null(ne1)) {
      if (is.null(alpha01) || is.null(bar.ye1) || is.null(se1)) {
        stop("'alpha01', 'bar.ye1', and 'se1' must all be non-NULL when 'ne1' is provided")
      }
      if (!is.numeric(ne1) || length(ne1) != 1L || is.na(ne1) ||
          ne1 != floor(ne1) || ne1 < 1L) {
        stop("'ne1' must be a single positive integer")
      }
      if (!is.numeric(alpha01) || length(alpha01) != 1L || is.na(alpha01) ||
          alpha01 <= 0 || alpha01 > 1) {
        stop("'alpha01' must be a single numeric value in (0, 1]")
      }
      if (!is.numeric(bar.ye1) || length(bar.ye1) != 1L || is.na(bar.ye1)) {
        stop("'bar.ye1' must be a single numeric value")
      }
      if (!is.numeric(se1) || length(se1) != 1L || is.na(se1) || se1 <= 0) {
        stop("'se1' must be a single positive numeric value")
      }
    }
    if (!is.null(ne2)) {
      if (is.null(alpha02) || is.null(bar.ye2) || is.null(se2)) {
        stop("'alpha02', 'bar.ye2', and 'se2' must all be non-NULL when 'ne2' is provided")
      }
      if (!is.numeric(ne2) || length(ne2) != 1L || is.na(ne2) ||
          ne2 != floor(ne2) || ne2 < 1L) {
        stop("'ne2' must be a single positive integer")
      }
      if (!is.numeric(alpha02) || length(alpha02) != 1L || is.na(alpha02) ||
          alpha02 <= 0 || alpha02 > 1) {
        stop("'alpha02' must be a single numeric value in (0, 1]")
      }
      if (!is.numeric(bar.ye2) || length(bar.ye2) != 1L || is.na(bar.ye2)) {
        stop("'bar.ye2' must be a single numeric value")
      }
      if (!is.numeric(se2) || length(se2) != 1L || is.na(se2) || se2 <= 0) {
        stop("'se2' must be a single positive numeric value")
      }
    }
    if (!is.null(bar.y2) && (!is.numeric(bar.y2) || any(is.na(bar.y2)))) {
      stop("'bar.y2' must be a numeric scalar or vector with no missing values")
    }
    if (!is.null(s2) && (!is.numeric(s2) || any(is.na(s2)) || any(s2 <= 0))) {
      stop("'s2' must be a positive numeric scalar or vector with no missing values")
    }
  }

  if (design == "controlled") {
    if (is.null(bar.y2) || is.null(s2)) {
      stop("'bar.y2' and 's2' must be non-NULL when design = 'controlled'")
    }
    if (!is.numeric(bar.y2) || any(is.na(bar.y2))) {
      stop("'bar.y2' must be a numeric scalar or vector with no missing values")
    }
    if (!is.numeric(s2) || any(is.na(s2)) || any(s2 <= 0)) {
      stop("'s2' must be a positive numeric scalar or vector with no missing values")
    }
  }

  # ---------------------------------------------------------------------------
  # Compute posterior t-distribution parameters (vectorised over bar.y1, s1,
  # bar.y2, s2; all arithmetic operates element-wise on vectors).
  # ---------------------------------------------------------------------------

  if (design == 'external' && prior == 'vague') {

    # --- External design with vague prior ---
    # Group 1 (Treatment): apply power prior if external data available
    if (!is.null(ne1) && !is.null(alpha01)) {
      mu.t1          <- (alpha01 * ne1 * bar.ye1 + n1 * bar.y1) / (alpha01 * ne1 + n1)
      kappa.star.n1  <- alpha01 * ne1 + n1
      nu.t1          <- alpha01 * ne1 + n1 - 1
      sigma2.star.n1 <- (alpha01 * (ne1 - 1) * se1 ^ 2 + (n1 - 1) * s1 ^ 2 +
                           (alpha01 * ne1 * n1 * (bar.ye1 - bar.y1) ^ 2) /
                           (alpha01 * ne1 + n1)) / (alpha01 * ne1 + n1)
    } else {
      # No external treatment data: use vague prior
      mu.t1          <- bar.y1
      kappa.star.n1  <- n1
      nu.t1          <- n1 - 1
      sigma2.star.n1 <- s1 ^ 2
    }

    # Group 2 (Control): apply power prior if external data available
    if (!is.null(ne2) && !is.null(alpha02)) {
      mu.t2          <- (alpha02 * ne2 * bar.ye2 + n2 * bar.y2) / (alpha02 * ne2 + n2)
      kappa.star.n2  <- alpha02 * ne2 + n2
      nu.t2          <- alpha02 * ne2 + n2 - 1
      sigma2.star.n2 <- (alpha02 * (ne2 - 1) * se2 ^ 2 + (n2 - 1) * s2 ^ 2 +
                           (alpha02 * ne2 * n2 * (bar.ye2 - bar.y2) ^ 2) /
                           (alpha02 * ne2 + n2)) / (alpha02 * ne2 + n2)
    } else {
      # No external control data: use vague prior
      mu.t2          <- bar.y2
      kappa.star.n2  <- n2
      nu.t2          <- n2 - 1
      sigma2.star.n2 <- s2 ^ 2
    }

    # Scale parameters depend on probability type
    if (prob == 'posterior') {
      sd.t1 <- sqrt(sigma2.star.n1 / kappa.star.n1)
      sd.t2 <- sqrt(sigma2.star.n2 / kappa.star.n2)
    } else {
      sd.t1 <- sqrt((1 + 1 / kappa.star.n1) * sigma2.star.n1 / m1)
      sd.t2 <- sqrt((1 + 1 / kappa.star.n2) * sigma2.star.n2 / m2)
    }

  } else {

    # --- Controlled or uncontrolled designs (vague or N-Inv-Chisq prior) ---
    if (prior == 'N-Inv-Chisq') {

      # Updated precision parameters
      kappa.n1 <- kappa01 + n1
      kappa.n2 <- kappa02 + n2

      # Updated degrees of freedom
      nu.t1 <- nu01 + n1
      if (design == 'controlled') {
        nu.t2 <- nu02 + n2
      } else {
        nu.t2 <- nu.t1
      }

      # Posterior means
      mu.t1 <- (kappa01 * mu01 + n1 * bar.y1) / kappa.n1
      if (design == 'controlled') {
        mu.t2 <- (kappa02 * mu02 + n2 * bar.y2) / kappa.n2
      } else {
        mu.t2 <- mu02
      }

      # Posterior variance for group 1 (vectorised over bar.y1, s1)
      var.n1 <- (nu01 * sigma01 ^ 2 + (n1 - 1) * s1 ^ 2 +
                   n1 * kappa01 * (mu01 - bar.y1) ^ 2 / kappa.n1) / nu.t1

      # Posterior variance for group 2
      if (design == 'controlled') {
        var.n2 <- (nu02 * sigma02 ^ 2 + (n2 - 1) * s2 ^ 2 +
                     n2 * kappa02 * (mu02 - bar.y2) ^ 2 / kappa.n2) / nu.t2
      } else {
        var.n2 <- NULL
      }

      # Scale parameters depend on probability type
      if (prob == 'posterior') {
        sd.t1 <- sqrt(var.n1 / kappa.n1)
        if (design == 'controlled') {
          sd.t2 <- sqrt(var.n2 / kappa.n2)
        } else {
          sd.t2 <- sqrt(r) * sd.t1
        }
      } else {
        sd.t1 <- sqrt((1 + 1 / kappa.n1) * var.n1 / m1)
        if (design == 'controlled') {
          sd.t2 <- sqrt((1 + 1 / kappa.n2) * var.n2 / m2)
        } else {
          sd.t2 <- sqrt(r) * sd.t1
        }
      }

    } else {
      # Vague (Jeffreys) prior

      # Degrees of freedom
      nu.t1 <- n1 - 1
      if (design == 'controlled') {
        nu.t2 <- n2 - 1
      } else {
        nu.t2 <- nu.t1
      }

      # Posterior means (vectorised over bar.y1, bar.y2)
      mu.t1 <- bar.y1
      if (design == 'controlled') {
        mu.t2 <- bar.y2
      } else {
        mu.t2 <- mu02
      }

      # Scale parameters (vectorised over s1, s2)
      if (prob == 'posterior') {
        sd.t1 <- sqrt(s1 ^ 2 / n1)
        if (design == 'controlled') {
          sd.t2 <- sqrt(s2 ^ 2 / n2)
        } else {
          sd.t2 <- sqrt(r) * sd.t1
        }
      } else {
        sd.t1 <- sqrt((1 + 1 / n1) * s1 ^ 2 / m1)
        if (design == 'controlled') {
          sd.t2 <- sqrt((1 + 1 / n2) * s2 ^ 2 / m2)
        } else {
          sd.t2 <- sqrt(r) * sd.t1
        }
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Compute probability using the specified method.
  # mu.t1, mu.t2, sd.t1, sd.t2 may be vectors of length nsim when called from
  # pbayesdecisionprob1cont; each method receives the full vector and returns a
  # vector of the same length (no outer loop required).
  # ---------------------------------------------------------------------------
  if (CalcMethod == 'NI') {
    results <- ptdiff_NI(theta0, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2, lower.tail)
  } else if (CalcMethod == 'MC') {
    results <- ptdiff_MC(nMC, theta0, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2, lower.tail)
  } else {
    results <- ptdiff_MM(theta0, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2, lower.tail)
  }

  return(results)
}
