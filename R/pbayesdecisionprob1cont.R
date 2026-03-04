#' Bayesian Go/NoGo/Gray Decision Probabilities for Single Continuous Endpoint
#'
#' Evaluates Go/NoGo/Gray decision probabilities for a single continuous endpoint
#' via Monte Carlo simulation. Supports controlled (parallel control), uncontrolled
#' (single-arm with informative priors), and external control (power prior borrowing)
#' designs with both posterior and predictive probability approaches.
#'
#' @param nsim A positive integer specifying the number of Monte Carlo simulation replicates.
#' @param prob A character string specifying the probability type: \code{'posterior'} or
#'        \code{'predictive'}.
#' @param design A character string specifying the trial design: \code{'controlled'},
#'        \code{'uncontrolled'}, or \code{'external'}.
#' @param prior A character string specifying the prior distribution: \code{'vague'} or
#'        \code{'N-Inv-Chisq'}.
#' @param CalcMethod A character string specifying the computation method: \code{'NI'}
#'        (Numerical Integration), \code{'MC'} (Monte Carlo), or \code{'MM'}
#'        (Moment Matching).
#' @param theta.TV A numeric value representing the target value (TV) threshold for the
#'        Go decision. Required if \code{prob = 'posterior'}.
#' @param theta.MAV A numeric value representing the minimum acceptable value (MAV) threshold
#'        for the NoGo decision. Required if \code{prob = 'posterior'}.
#' @param theta.NULL A numeric value representing the null hypothesis threshold.
#'        Required if \code{prob = 'predictive'}.
#' @param nMC A positive integer specifying the number of Monte Carlo draws for computing
#'        posterior probabilities. Required if \code{CalcMethod = 'MC'}.
#' @param gamma1 A numeric value in (0, 1) representing the Go decision threshold.
#' @param gamma2 A numeric value in (0, 1) representing the NoGo decision threshold.
#'        Must be strictly less than \code{gamma1}.
#' @param n1 A positive integer representing the sample size for group 1 (treatment).
#' @param n2 A positive integer representing the sample size for group 2 (control).
#'        Required if \code{design = 'controlled'} or \code{design = 'external'}.
#'        Set to \code{NULL} if \code{design = 'uncontrolled'}.
#' @param m1 A positive integer representing the future sample size for group 1.
#'        Required if \code{prob = 'predictive'}.
#' @param m2 A positive integer representing the future sample size for group 2.
#'        Required if \code{prob = 'predictive'}.
#' @param kappa01 A positive numeric value representing the prior hyperparameter kappa for
#'        group 1. Required if \code{prior = 'N-Inv-Chisq'}.
#' @param kappa02 A positive numeric value representing the prior hyperparameter kappa for
#'        group 2. Required if \code{prior = 'N-Inv-Chisq'} and
#'        \code{design \%in\% c('controlled', 'external')}.
#' @param nu01 A positive numeric value representing the prior hyperparameter nu for
#'        group 1. Required if \code{prior = 'N-Inv-Chisq'}.
#' @param nu02 A positive numeric value representing the prior hyperparameter nu for
#'        group 2. Required if \code{prior = 'N-Inv-Chisq'} and
#'        \code{design \%in\% c('controlled', 'external')}.
#' @param mu01 A numeric value representing the prior mean for group 1. Required if
#'        \code{prior = 'N-Inv-Chisq'}.
#' @param mu02 A numeric value representing the prior mean for group 2. For
#'        \code{design = 'uncontrolled'}, this represents the hypothetical control mean.
#'        Required if \code{prior = 'N-Inv-Chisq'} and
#'        \code{design \%in\% c('controlled', 'external')}, or if
#'        \code{design = 'uncontrolled'}.
#' @param sigma01 A positive numeric value representing the prior standard deviation for
#'        group 1. Required if \code{prior = 'N-Inv-Chisq'}.
#' @param sigma02 A positive numeric value representing the prior standard deviation for
#'        group 2. Required if \code{prior = 'N-Inv-Chisq'} and
#'        \code{design \%in\% c('controlled', 'external')}.
#' @param mu1 A numeric value representing the true mean for group 1 in the simulation.
#' @param mu2 A numeric value representing the true mean for group 2 in the simulation.
#'        For uncontrolled design, this represents the historical control mean.
#' @param sigma1 A positive numeric value representing the true standard deviation
#'        for group 1 in the simulation.
#' @param sigma2 A positive numeric value representing the true standard deviation
#'        for group 2 in the simulation. For uncontrolled design, this represents
#'        the historical control standard deviation.
#' @param r A positive numeric value representing the variance scaling factor that allows
#'        the scale of hypothetical control to be different from treatment. Specifically,
#'        \code{sd.control = sqrt(r) * sd.treatment}. Required if \code{design = 'uncontrolled'}.
#'        When \code{r = 1}, the control and treatment have the same variance scale.
#' @param ne1 A positive integer representing the number of patients in group 1 for
#'        the external data. Required if \code{design = 'external'} and external
#'        treatment data are available.
#' @param ne2 A positive integer representing the number of patients in group 2 for
#'        the external data. Required if \code{design = 'external'} and external
#'        control data are available.
#' @param alpha01 A numeric value in (0, 1] representing the power prior scale parameter
#'        for group 1. Controls the degree of borrowing from external treatment data:
#'        0 = no borrowing, 1 = full borrowing. Required if \code{ne1} is specified.
#' @param alpha02 A numeric value in (0, 1] representing the power prior scale parameter
#'        for group 2. Controls the degree of borrowing from external control data:
#'        0 = no borrowing, 1 = full borrowing. Required if \code{ne2} is specified.
#' @param bar.ye1 A numeric value representing the sample mean of the external data
#'        for group 1. Required if \code{ne1} is specified.
#' @param bar.ye2 A numeric value representing the sample mean of the external data
#'        for group 2. Required if \code{ne2} is specified.
#' @param se1 A positive numeric value representing the sample standard deviation
#'        of the external data for group 1. Required if \code{ne1} is specified.
#' @param se2 A positive numeric value representing the sample standard deviation
#'        of the external data for group 2. Required if \code{ne2} is specified.
#' @param error_if_Miss A logical value; if \code{TRUE} (default), the function stops
#'        with an error when positive Miss probability is obtained, indicating poorly
#'        chosen thresholds. If \code{FALSE}, the function proceeds and reports Miss
#'        probability based on \code{Gray_inc_Miss} setting.
#' @param Gray_inc_Miss A logical value; if \code{TRUE}, Miss probability is included
#'        in Gray probability (Miss is not reported separately). If \code{FALSE}
#'        (default), Miss probability is reported as a separate category. This parameter
#'        is only active when \code{error_if_Miss = FALSE}.
#' @param seed A numeric value representing the seed number for reproducible random number generation.
#'
#' @return A data frame containing the true means for both groups and the Go, NoGo, and
#'         Gray probabilities. When \code{error_if_Miss = FALSE} and \code{Gray_inc_Miss = FALSE},
#'         Miss probability is also included as a separate column. For uncontrolled design,
#'         only mu1 is included (not mu2).
#'
#' @details
#' The function performs Monte Carlo simulation to evaluate operating characteristics by:
#' \itemize{
#'   \item Generating random trial data based on specified true parameters
#'   \item Computing posterior or predictive probabilities for each simulated trial
#'   \item Classifying each trial as Go, NoGo, or Gray based on decision thresholds
#' }
#'
#' Posterior parameter calculations (mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2) are
#' fully vectorized over the nsim simulated datasets. pbayespostpred1cont is called
#' twice (once per threshold), with each call receiving vectors of length nsim and
#' returning a vector of nsim probabilities - no inner loop over simulation replicates
#' is required.
#'
#' For external control designs, power priors are incorporated using exact conjugate representation:
#' \itemize{
#'   \item Power priors for normal data are mathematically equivalent to Normal-Inverse-Chi-squared distributions
#'   \item This enables closed-form computation without MCMC sampling
#'   \item Alpha parameters control the degree of borrowing (0 = no borrowing, 1 = full borrowing)
#' }
#'
#' **Decision rules**:
#' \itemize{
#'   \item **Go**: P(treatment effect > threshold) >= gamma1
#'   \item **NoGo**: P(treatment effect > threshold) <= gamma2
#'   \item **Gray**: gamma2 < P(treatment effect > threshold) < gamma1
#'   \item **Miss**: Both Go and NoGo criteria are met simultaneously (indicates
#'                   poorly chosen thresholds)
#' }
#'
#' **Handling Miss probability**:
#' \itemize{
#'   \item When \code{error_if_Miss = TRUE} (default): Function stops with error if
#'         Miss probability > 0, prompting reconsideration of thresholds
#'   \item When \code{error_if_Miss = FALSE} and \code{Gray_inc_Miss = TRUE}: Miss
#'         probability is added to Gray probability
#'   \item When \code{error_if_Miss = FALSE} and \code{Gray_inc_Miss = FALSE}: Miss
#'         probability is reported as a separate category
#' }
#'
#' The function can be used for:
#' \itemize{
#'   \item **Controlled design**: Two-arm randomized trial
#'   \item **Uncontrolled design**: Single-arm trial with informative priors (historical control)
#'   \item **External design**: Incorporating historical data through power priors
#' }
#'
#' @examples
#' # Example 1: Controlled design with vague prior and NI method
#' # (default: error_if_Miss = TRUE, Gray_inc_Miss = FALSE)
#' pbayesdecisionprob1cont(
#'   nsim = 100, prob = 'posterior', design = 'controlled', prior = 'vague', CalcMethod = 'NI',
#'   theta.TV = 1.5, theta.MAV = -0.5, theta.NULL = NULL,
#'   nMC = NULL, gamma1 = 0.7, gamma2 = 0.2,
#'   n1 = 15, n2 = 15, m1 = NULL, m2 = NULL,
#'   kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
#'   mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
#'   mu1 = 3, mu2 = 1, sigma1 = 1.2, sigma2 = 1.1,
#'   r = NULL, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
#'   bar.ye1 = NULL, bar.ye2 = NULL, se1 = NULL, se2 = NULL,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed = 2
#' )
#'
#' # Example 2: Uncontrolled design with informative prior
#' pbayesdecisionprob1cont(
#'   nsim = 100, prob = 'posterior', design = 'uncontrolled', prior = 'N-Inv-Chisq', CalcMethod = 'NI',
#'   theta.TV = 1.0, theta.MAV = 0.0, theta.NULL = NULL,
#'   nMC = NULL, gamma1 = 0.8, gamma2 = 0.2,
#'   n1 = 20, n2 = NULL, m1 = NULL, m2 = NULL,
#'   kappa01 = 2, kappa02 = NULL, nu01 = 5, nu02 = NULL,
#'   mu01 = 3.0, mu02 = 1.5, sigma01 = 1.5, sigma02 = NULL,
#'   mu1 = 3.5, mu2 = NULL, sigma1 = 1.3, sigma2 = NULL,
#'   r = 1, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
#'   bar.ye1 = NULL, bar.ye2 = NULL, se1 = NULL, se2 = NULL,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed = 3
#' )
#'
#' # Example 3: External design with control data using MM approximation
#' \dontrun{
#' pbayesdecisionprob1cont(
#'   nsim = 100, prob = 'posterior', design = 'external', prior = 'vague', CalcMethod = 'MM',
#'   theta.TV = 1.0, theta.MAV = 0.0, theta.NULL = NULL,
#'   nMC = NULL, gamma1 = 0.8, gamma2 = 0.2,
#'   n1 = 12, n2 = 12, m1 = NULL, m2 = NULL,
#'   kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
#'   mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
#'   mu1 = 2, mu2 = 0, sigma1 = 1, sigma2 = 1,
#'   r = NULL, ne1 = NULL, ne2 = 20, alpha01 = NULL, alpha02 = 0.5,
#'   bar.ye1 = NULL, bar.ye2 = 0, se1 = NULL, se2 = 1,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed = 4
#' )
#' }
#'
#' # Example 4: Controlled design with predictive probability
#' pbayesdecisionprob1cont(
#'   nsim = 100, prob = 'predictive', design = 'controlled', prior = 'N-Inv-Chisq', CalcMethod = 'NI',
#'   theta.TV = NULL, theta.MAV = NULL, theta.NULL = 2.0,
#'   nMC = NULL, gamma1 = 0.75, gamma2 = 0.35,
#'   n1 = 15, n2 = 15, m1 = 50, m2 = 50,
#'   kappa01 = 3, kappa02 = 3, nu01 = 4, nu02 = 4,
#'   mu01 = 3.5, mu02 = 1.5, sigma01 = 1.5, sigma02 = 1.5,
#'   mu1 = 3.2, mu2 = 1.3, sigma1 = 1.4, sigma2 = 1.2,
#'   r = NULL, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
#'   bar.ye1 = NULL, bar.ye2 = NULL, se1 = NULL, se2 = NULL,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed = 5
#' )
#'
#' # Example 5: External design with predictive probability using MC method
#' \dontrun{
#' pbayesdecisionprob1cont(
#'   nsim = 100, prob = 'predictive', design = 'external', prior = 'vague', CalcMethod = 'MC',
#'   theta.TV = NULL, theta.MAV = NULL, theta.NULL = 1.5,
#'   nMC = 5000, gamma1 = 0.7, gamma2 = 0.4,
#'   n1 = 12, n2 = 12, m1 = 30, m2 = 30,
#'   kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
#'   mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
#'   mu1 = 2.5, mu2 = 1.0, sigma1 = 1.3, sigma2 = 1.1,
#'   r = NULL, ne1 = 15, ne2 = 18, alpha01 = 0.6, alpha02 = 0.7,
#'   bar.ye1 = 2.3, bar.ye2 = 0.9, se1 = 1.2, se2 = 1.0,
#'   error_if_Miss = FALSE, Gray_inc_Miss = FALSE, seed = 6
#' )
#' }
#'
#' # Example 6: Report Miss probability separately when thresholds may be suboptimal
#' pbayesdecisionprob1cont(
#'   nsim = 100, prob = 'posterior', design = 'controlled', prior = 'vague', CalcMethod = 'NI',
#'   theta.TV = 1.0, theta.MAV = 0.8, theta.NULL = NULL,
#'   nMC = NULL, gamma1 = 0.65, gamma2 = 0.55,
#'   n1 = 10, n2 = 10, m1 = NULL, m2 = NULL,
#'   kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
#'   mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
#'   mu1 = 2.5, mu2 = 1.5, sigma1 = 1.0, sigma2 = 1.0,
#'   r = NULL, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
#'   bar.ye1 = NULL, bar.ye2 = NULL, se1 = NULL, se2 = NULL,
#'   error_if_Miss = FALSE, Gray_inc_Miss = FALSE, seed = 7
#' )
#'
#' # Example 7: Include Miss probability in Gray when error_if_Miss = FALSE
#' pbayesdecisionprob1cont(
#'   nsim = 100, prob = 'posterior', design = 'controlled', prior = 'vague', CalcMethod = 'NI',
#'   theta.TV = 1.0, theta.MAV = 0.8, theta.NULL = NULL,
#'   nMC = NULL, gamma1 = 0.65, gamma2 = 0.55,
#'   n1 = 10, n2 = 10, m1 = NULL, m2 = NULL,
#'   kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
#'   mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
#'   mu1 = 2.5, mu2 = 1.5, sigma1 = 1.0, sigma2 = 1.0,
#'   r = NULL, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
#'   bar.ye1 = NULL, bar.ye2 = NULL, se1 = NULL, se2 = NULL,
#'   error_if_Miss = FALSE, Gray_inc_Miss = TRUE, seed = 8
#' )
#'
#' @importFrom stats rnorm
#' @export
pbayesdecisionprob1cont <- function(nsim, prob, design, prior, CalcMethod, theta.TV, theta.MAV, theta.NULL,
                                    nMC = NULL, gamma1, gamma2, n1, n2 = NULL, m1, m2, kappa01, kappa02, nu01, nu02,
                                    mu01, mu02, sigma01, sigma02, mu1, mu2 = NULL, sigma1, sigma2 = NULL,
                                    r = NULL, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
                                    bar.ye1 = NULL, bar.ye2 = NULL, se1 = NULL, se2 = NULL,
                                    error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed) {

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

  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) ||
      nsim != floor(nsim) || nsim < 1L) {
    stop("'nsim' must be a single positive integer")
  }

  # Validate n1 (always required)
  if (!is.numeric(n1) || length(n1) != 1L || is.na(n1) ||
      n1 != floor(n1) || n1 < 1L) {
    stop("'n1' must be a single positive integer")
  }

  # Validate n2 (required only for controlled and external designs)
  if (design %in% c("controlled", "external")) {
    if (!is.numeric(n2) || length(n2) != 1L || is.na(n2) ||
        n2 != floor(n2) || n2 < 1L) {
      stop("'n2' must be a single positive integer")
    }
  }

  if (!is.numeric(gamma1) || length(gamma1) != 1L || is.na(gamma1) ||
      gamma1 <= 0 || gamma1 >= 1) {
    stop("'gamma1' must be a single numeric value in (0, 1)")
  }

  if (!is.numeric(gamma2) || length(gamma2) != 1L || is.na(gamma2) ||
      gamma2 <= 0 || gamma2 >= 1) {
    stop("'gamma2' must be a single numeric value in (0, 1)")
  }

  if (gamma2 >= gamma1) {
    stop("'gamma2' must be strictly less than 'gamma1'")
  }

  if (!is.numeric(sigma1) || length(sigma1) != 1L || is.na(sigma1) || sigma1 <= 0) {
    stop("'sigma1' must be a single positive numeric value")
  }

  if (!is.logical(error_if_Miss) || length(error_if_Miss) != 1L || is.na(error_if_Miss)) {
    stop("'error_if_Miss' must be a single logical value (TRUE or FALSE)")
  }

  if (!is.logical(Gray_inc_Miss) || length(Gray_inc_Miss) != 1L || is.na(Gray_inc_Miss)) {
    stop("'Gray_inc_Miss' must be a single logical value (TRUE or FALSE)")
  }

  # Validate CalcMethod-specific parameter
  if (CalcMethod == "MC") {
    if (is.null(nMC)) {
      stop("'nMC' must be non-NULL when CalcMethod = 'MC'")
    }
    if (!is.numeric(nMC) || length(nMC) != 1L || is.na(nMC) ||
        nMC != floor(nMC) || nMC < 1L) {
      stop("'nMC' must be a single positive integer")
    }
  }

  # Validate prob-specific parameters
  if (prob == "posterior") {
    if (is.null(theta.TV) || is.null(theta.MAV)) {
      stop("'theta.TV' and 'theta.MAV' must be non-NULL when prob = 'posterior'")
    }
    if (!is.numeric(theta.TV) || length(theta.TV) != 1L || is.na(theta.TV)) {
      stop("'theta.TV' must be a single numeric value")
    }
    if (!is.numeric(theta.MAV) || length(theta.MAV) != 1L || is.na(theta.MAV)) {
      stop("'theta.MAV' must be a single numeric value")
    }
    if (theta.TV <= theta.MAV) {
      stop("'theta.TV' must be strictly greater than 'theta.MAV'")
    }
  } else {
    if (is.null(theta.NULL)) {
      stop("'theta.NULL' must be non-NULL when prob = 'predictive'")
    }
    if (!is.numeric(theta.NULL) || length(theta.NULL) != 1L || is.na(theta.NULL)) {
      stop("'theta.NULL' must be a single numeric value")
    }
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

  # Validate prior-specific parameters
  if (prior == "N-Inv-Chisq") {
    if (is.null(kappa01) || is.null(nu01) || is.null(mu01) || is.null(sigma01)) {
      stop("'kappa01', 'nu01', 'mu01', and 'sigma01' must be non-NULL when prior = 'N-Inv-Chisq'")
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
    if (design %in% c("controlled", "external")) {
      for (nm in c("kappa02", "nu02", "mu02", "sigma02")) {
        val <- get(nm)
        if (is.null(val)) {
          stop(paste0("'", nm, "' must be non-NULL when prior = 'N-Inv-Chisq' and design = '",
                      design, "'"))
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
    }
  }

  # Validate mu1 (true treatment mean vector)
  if (!is.numeric(mu1) || length(mu1) < 1L || any(is.na(mu1))) {
    stop("'mu1' must be a numeric scalar or vector with no missing values")
  }

  # Set seed for reproducible results
  set.seed(seed)

  # Generate nsim x n1 matrix of standardized residuals for group 1 (mean=0, sd=sigma1).
  # Shared across all mu1 scenarios: bar.y1[scenario] = rowMeans(Z1) + mu1[scenario],
  # and s1 is independent of mu1, so it is computed once and reused.
  Z1 <- matrix(rnorm(nsim * n1, mean = 0, sd = sigma1), nrow = nsim)
  Z1_rowmean <- rowSums(Z1) / n1
  s1 <- sqrt(rowSums((Z1 - Z1_rowmean) ^ 2) / (n1 - 1))

  if (design == 'controlled' | design == 'external') {
    # Generate nsim x n2 matrix of standardized residuals for group 2 (mean=0, sd=sigma2).
    # mu2 is fixed to a single value, so bar.y2 and s2 are computed once.
    Z2 <- matrix(rnorm(nsim * n2, mean = 0, sd = sigma2), nrow = nsim)
    Z2_rowmean <- rowSums(Z2) / n2
    bar.y2 <- Z2_rowmean + mu2
    s2     <- sqrt(rowSums((Z2 - Z2_rowmean) ^ 2) / (n2 - 1))
  } else if (design == 'uncontrolled') {
    # For uncontrolled design: no group 2 data generated
    bar.y2 <- NULL
    s2     <- NULL
  }

  # External data are fixed historical sample statistics (no generation needed)

  # Ensure mu1 is a numeric vector to support multiple scenarios
  mu1 <- as.numeric(mu1)
  n_scenarios <- length(mu1)

  # Set threshold values based on probability type
  if (prob == 'posterior') {
    theta0_Go   <- theta.TV
    theta0_NoGo <- theta.MAV
  } else {
    theta0_Go   <- theta.NULL
    theta0_NoGo <- theta.NULL
  }

  # Expand bar.y1 and s1 across all mu1 scenarios into a single vector of
  # length nsim * n_scenarios. For bar.y1, shift each scenario's block by mu1[i].
  # s1 is independent of mu1 so it is simply replicated.
  bar.y1_all <- rep(Z1_rowmean, times = n_scenarios) + rep(mu1, each = nsim)
  s1_all     <- rep(s1, times = n_scenarios)

  if (design == 'controlled' | design == 'external') {
    bar.y2_all <- rep(bar.y2, times = n_scenarios)
    s2_all     <- rep(s2,     times = n_scenarios)
  } else {
    bar.y2_all <- NULL
    s2_all     <- NULL
  }

  # Call pbayespostpred1cont once with the expanded vectors.
  # Returns a vector of length nsim * n_scenarios.
  gPost_Go_all <- pbayespostpred1cont(
    prob = prob, design = design, prior = prior, CalcMethod = CalcMethod,
    theta0 = theta0_Go, nMC = nMC, n1 = n1, n2 = n2, m1 = m1, m2 = m2,
    kappa01 = kappa01, kappa02 = kappa02, nu01 = nu01, nu02 = nu02,
    mu01 = mu01, mu02 = mu02, sigma01 = sigma01, sigma02 = sigma02,
    bar.y1 = bar.y1_all, bar.y2 = bar.y2_all, s1 = s1_all, s2 = s2_all,
    r = r, ne1 = ne1, ne2 = ne2, alpha01 = alpha01, alpha02 = alpha02,
    bar.ye1 = bar.ye1, bar.ye2 = bar.ye2, se1 = se1, se2 = se2,
    lower.tail = FALSE
  )

  gPost_NoGo_all <- pbayespostpred1cont(
    prob = prob, design = design, prior = prior, CalcMethod = CalcMethod,
    theta0 = theta0_NoGo, nMC = nMC, n1 = n1, n2 = n2, m1 = m1, m2 = m2,
    kappa01 = kappa01, kappa02 = kappa02, nu01 = nu01, nu02 = nu02,
    mu01 = mu01, mu02 = mu02, sigma01 = sigma01, sigma02 = sigma02,
    bar.y1 = bar.y1_all, bar.y2 = bar.y2_all, s1 = s1_all, s2 = s2_all,
    r = r, ne1 = ne1, ne2 = ne2, alpha01 = alpha01, alpha02 = alpha02,
    bar.ye1 = bar.ye1, bar.ye2 = bar.ye2, se1 = se1, se2 = se2,
    lower.tail = TRUE
  )

  # Reshape results into nsim x n_scenarios matrices and compute decision
  # probabilities per scenario using colMeans (no R-level loop needed).
  probs_Go_mat   <- matrix((gPost_Go_all >= gamma1) & (gPost_NoGo_all <  gamma2),
                           nrow = nsim)
  probs_NoGo_mat <- matrix((gPost_Go_all <  gamma1) & (gPost_NoGo_all >= gamma2),
                           nrow = nsim)
  probs_Miss_mat <- matrix((gPost_Go_all >= gamma1) & (gPost_NoGo_all >= gamma2),
                           nrow = nsim)

  GoNogoProb <- cbind(
    colMeans(probs_Go_mat),
    colMeans(probs_NoGo_mat),
    colMeans(probs_Miss_mat)
  )  # n_scenarios x 3 matrix

  # Check for positive Miss probabilities (indicates inappropriate thresholds)
  if (error_if_Miss) {
    if (any(GoNogoProb[, 3] > 0)) {
      stop("Positive Miss probability detected. Please re-consider the chosen thresholds.")
    }
  }

  # Calculate Gray probability for each scenario (complement of Go and NoGo)
  if (Gray_inc_Miss) {
    # Include Miss in Gray probability
    GrayProb <- 1 - rowSums(GoNogoProb[, -3, drop = FALSE])
  } else {
    # Exclude Miss from Gray probability
    GrayProb <- 1 - rowSums(GoNogoProb)
  }

  # Create results data frame
  if (design == 'uncontrolled') {
    # For uncontrolled design, only include mu1
    results <- data.frame(
      mu1, Go = GoNogoProb[, 1], Gray = GrayProb, NoGo = GoNogoProb[, 2]
    )
  } else {
    # For controlled and external designs, include both mu1 and mu2
    results <- data.frame(
      mu1, mu2, Go = GoNogoProb[, 1], Gray = GrayProb, NoGo = GoNogoProb[, 2]
    )
  }

  # Add Miss column when error_if_Miss is FALSE and Gray_inc_Miss is FALSE
  if (!error_if_Miss) {
    if (!Gray_inc_Miss) {
      results$Miss <- GoNogoProb[, 3]
    }
  }

  # Address floating point error (apply only to probability columns, not mu1/mu2)
  prob_cols <- names(results)[!names(results) %in% c('mu1', 'mu2')]
  results[prob_cols] <- lapply(results[prob_cols], function(col) {
    ifelse(col < .Machine$double.eps ^ 0.25, 0, col)
  })

  # Attach metadata as attributes for use in print()
  attr(results, 'prob')          <- prob
  attr(results, 'design')        <- design
  attr(results, 'prior')         <- prior
  attr(results, 'CalcMethod')    <- CalcMethod
  attr(results, 'nsim')          <- nsim
  attr(results, 'nMC')           <- nMC
  attr(results, 'theta.TV')      <- theta.TV
  attr(results, 'theta.MAV')     <- theta.MAV
  attr(results, 'theta.NULL')    <- theta.NULL
  attr(results, 'gamma1')        <- gamma1
  attr(results, 'gamma2')        <- gamma2
  attr(results, 'n1')            <- n1
  attr(results, 'n2')            <- n2
  attr(results, 'sigma1')        <- sigma1
  attr(results, 'sigma2')        <- sigma2
  attr(results, 'r')             <- r
  attr(results, 'm1')            <- m1
  attr(results, 'm2')            <- m2
  attr(results, 'kappa01')       <- kappa01
  attr(results, 'kappa02')       <- kappa02
  attr(results, 'nu01')          <- nu01
  attr(results, 'nu02')          <- nu02
  attr(results, 'mu01')          <- mu01
  attr(results, 'mu02')          <- mu02
  attr(results, 'sigma01')       <- sigma01
  attr(results, 'sigma02')       <- sigma02
  attr(results, 'ne1')           <- ne1
  attr(results, 'ne2')           <- ne2
  attr(results, 'alpha01')       <- alpha01
  attr(results, 'alpha02')       <- alpha02
  attr(results, 'bar.ye1')       <- bar.ye1
  attr(results, 'bar.ye2')       <- bar.ye2
  attr(results, 'se1')           <- se1
  attr(results, 'se2')           <- se2
  attr(results, 'error_if_Miss') <- error_if_Miss
  attr(results, 'Gray_inc_Miss') <- Gray_inc_Miss
  attr(results, 'seed')          <- seed

  # Assign S3 class
  class(results) <- c('pbayesdecisionprob1cont', 'data.frame')

  return(results)
}

#' Print Method for pbayesdecisionprob1cont Objects
#'
#' Displays a formatted summary of Go/NoGo/Gray decision probabilities
#' for continuous endpoint results returned by \code{\link{pbayesdecisionprob1cont}}.
#'
#' @param x An object of class \code{pbayesdecisionprob1cont}.
#' @param digits A positive integer specifying the number of decimal places
#'        for probability values. Default is 4.
#' @param ... Further arguments passed to or from other methods (ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.pbayesdecisionprob1cont <- function(x, digits = 4, ...) {
  # Helper to format a value as string (NULL -> "NULL")
  fmt <- function(v) if(is.null(v)) 'NULL' else as.character(v)

  # Extract metadata from attributes
  prob          <- attr(x, 'prob')
  design        <- attr(x, 'design')
  prior         <- attr(x, 'prior')
  CalcMethod    <- attr(x, 'CalcMethod')
  nsim          <- attr(x, 'nsim')
  nMC           <- attr(x, 'nMC')
  gamma1        <- attr(x, 'gamma1')
  gamma2        <- attr(x, 'gamma2')
  n1            <- attr(x, 'n1')
  n2            <- attr(x, 'n2')
  sigma1        <- attr(x, 'sigma1')
  sigma2        <- attr(x, 'sigma2')
  r             <- attr(x, 'r')
  m1            <- attr(x, 'm1')
  m2            <- attr(x, 'm2')
  kappa01       <- attr(x, 'kappa01')
  kappa02       <- attr(x, 'kappa02')
  nu01          <- attr(x, 'nu01')
  nu02          <- attr(x, 'nu02')
  mu01          <- attr(x, 'mu01')
  mu02          <- attr(x, 'mu02')
  sigma01       <- attr(x, 'sigma01')
  sigma02       <- attr(x, 'sigma02')
  ne1           <- attr(x, 'ne1')
  ne2           <- attr(x, 'ne2')
  alpha01       <- attr(x, 'alpha01')
  alpha02       <- attr(x, 'alpha02')
  bar.ye1       <- attr(x, 'bar.ye1')
  bar.ye2       <- attr(x, 'bar.ye2')
  se1           <- attr(x, 'se1')
  se2           <- attr(x, 'se2')
  error_if_Miss <- attr(x, 'error_if_Miss')
  Gray_inc_Miss <- attr(x, 'Gray_inc_Miss')
  seed          <- attr(x, 'seed')

  # Build threshold string based on probability type
  if(prob == 'posterior') {
    theta_str <- sprintf('TV = %s, MAV = %s',
                         fmt(attr(x, 'theta.TV')), fmt(attr(x, 'theta.MAV')))
  } else {
    theta_str <- sprintf('NULL = %s', fmt(attr(x, 'theta.NULL')))
  }

  # Print header
  cat('Go/NoGo/Gray Decision Probabilities (Single Continuous Endpoint)\n')
  cat(strrep('-', 60), '\n')
  cat(sprintf('  Probability type : %s\n',   prob))
  cat(sprintf('  Design           : %s\n',   design))
  cat(sprintf('  Prior            : %s\n',   prior))
  cat(sprintf('  Calc method      : %s\n',   CalcMethod))
  cat(sprintf('  Simulations      : nsim = %s\n', fmt(nsim)))
  if (!is.null(nMC)) {
    cat(sprintf('  MC draws         : nMC = %s\n', fmt(nMC)))
  }
  cat(sprintf('  Threshold(s)     : %s\n',   theta_str))
  cat(sprintf('  Go  threshold    : gamma1 = %s\n', fmt(gamma1)))
  cat(sprintf('  NoGo threshold   : gamma2 = %s\n', fmt(gamma2)))
  cat(sprintf('  Sample size      : n1 = %s, n2 = %s\n', fmt(n1), fmt(n2)))
  cat(sprintf('  True SD          : sigma1 = %s, sigma2 = %s\n', fmt(sigma1), fmt(sigma2)))
  if (design == 'uncontrolled') {
    cat(sprintf('  Variance ratio   : r = %s\n', fmt(r)))
  }
  if (!is.null(m1) || !is.null(m2)) {
    cat(sprintf('  Future size      : m1 = %s, m2 = %s\n', fmt(m1), fmt(m2)))
  }
  if (prior == 'N-Inv-Chisq') {
    cat(sprintf('  Prior group 1    : kappa01 = %s, nu01 = %s, mu01 = %s, sigma01 = %s\n',
                fmt(kappa01), fmt(nu01), fmt(mu01), fmt(sigma01)))
    cat(sprintf('  Prior group 2    : kappa02 = %s, nu02 = %s, mu02 = %s, sigma02 = %s\n',
                fmt(kappa02), fmt(nu02), fmt(mu02), fmt(sigma02)))
  }
  if (design == 'external') {
    cat(sprintf('  External grp 1   : ne1 = %s, alpha01 = %s, bar.ye1 = %s, se1 = %s\n',
                fmt(ne1), fmt(alpha01), fmt(bar.ye1), fmt(se1)))
    cat(sprintf('  External grp 2   : ne2 = %s, alpha02 = %s, bar.ye2 = %s, se2 = %s\n',
                fmt(ne2), fmt(alpha02), fmt(bar.ye2), fmt(se2)))
  }
  cat(sprintf('  error_if_Miss    : %s\n', fmt(error_if_Miss)))
  cat(sprintf('  Gray_inc_Miss    : %s\n', fmt(Gray_inc_Miss)))
  cat(sprintf('  Seed             : %s\n', fmt(seed)))
  cat(strrep('-', 60), '\n')

  # Print results table
  df_print        <- as.data.frame(x)
  prob_cols       <- names(df_print)[!names(df_print) %in% c('mu1', 'mu2')]
  df_print[prob_cols] <- lapply(df_print[prob_cols],
                                function(col) round(col, digits))
  print.data.frame(df_print, row.names = FALSE)

  invisible(x)
}
