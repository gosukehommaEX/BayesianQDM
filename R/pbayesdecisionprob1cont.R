#' Bayesian Go/NoGo/Gray Decision Probabilities for Single Continuous Endpoint
#'
#' Evaluates Go/NoGo/Gray decision probabilities for a single continuous endpoint
#' via Monte Carlo simulation. Supports controlled (parallel control), uncontrolled
#' (single-arm with informative priors), and external (power prior borrowing)
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
#' @param theta_TV A numeric value representing the target value (TV) threshold for the
#'        Go decision. Required if \code{prob = 'posterior'}; otherwise \code{NULL}.
#'        Default is \code{NULL}.
#' @param theta_MAV A numeric value representing the minimum acceptable value (MAV) threshold
#'        for the NoGo decision. Required if \code{prob = 'posterior'}; otherwise \code{NULL}.
#'        Default is \code{NULL}.
#' @param theta_NULL A numeric value representing the null hypothesis threshold.
#'        Required if \code{prob = 'predictive'}; otherwise \code{NULL}.
#'        Default is \code{NULL}.
#' @param nMC A positive integer specifying the number of Monte Carlo draws for computing
#'        posterior probabilities. Required if \code{CalcMethod = 'MC'};
#'        otherwise \code{NULL}.
#' @param gamma_go A numeric value in (0, 1) representing the Go decision threshold.
#' @param gamma_nogo A numeric value in (0, 1) representing the NoGo decision threshold.
#'        Must be strictly less than \code{gamma_go}.
#' @param n_t A positive integer representing the sample size for the treatment group.
#' @param n_c A positive integer representing the sample size for the control group.
#'        Required if \code{design = 'controlled'} or \code{design = 'external'}.
#'        Set to \code{NULL} if \code{design = 'uncontrolled'}.
#' @param m_t A positive integer representing the future sample size for the
#'        treatment group. Required if \code{prob = 'predictive'};
#'        otherwise \code{NULL}. Default is \code{NULL}.
#' @param m_c A positive integer representing the future sample size for the control group.
#'        Required if \code{prob = 'predictive'}; otherwise \code{NULL}.
#'        Default is \code{NULL}.
#' @param kappa0_t A positive numeric value representing the prior hyperparameter kappa for
#'        the treatment group. Required if \code{prior = 'N-Inv-Chisq'}; otherwise \code{NULL}.
#'        Default is \code{NULL}.
#' @param kappa0_c A positive numeric value representing the prior hyperparameter kappa for
#'        the control group. Required if \code{prior = 'N-Inv-Chisq'} and
#'        \code{design \%in\% c('controlled', 'external')}; otherwise \code{NULL}.
#'        Default is \code{NULL}.
#' @param nu0_t A positive numeric value representing the prior hyperparameter nu for
#'        the treatment group. Required if \code{prior = 'N-Inv-Chisq'}; otherwise \code{NULL}.
#'        Default is \code{NULL}.
#' @param nu0_c A positive numeric value representing the prior hyperparameter nu for
#'        the control group. Required if \code{prior = 'N-Inv-Chisq'} and
#'        \code{design \%in\% c('controlled', 'external')}; otherwise \code{NULL}.
#'        Default is \code{NULL}.
#' @param mu0_t A numeric value representing the prior mean for the treatment group. Required if
#'        \code{prior = 'N-Inv-Chisq'}; otherwise \code{NULL}. Default is \code{NULL}.
#' @param mu0_c A numeric value representing the prior mean for the control group. For
#'        \code{design = 'uncontrolled'}, this represents the hypothetical control mean.
#'        Required if \code{prior = 'N-Inv-Chisq'} and
#'        \code{design \%in\% c('controlled', 'external')}, or if
#'        \code{design = 'uncontrolled'}; otherwise \code{NULL}. Default is \code{NULL}.
#' @param sigma0_t A positive numeric value representing the prior standard deviation for
#'        the treatment group. Required if \code{prior = 'N-Inv-Chisq'}; otherwise \code{NULL}.
#'        Default is \code{NULL}.
#' @param sigma0_c A positive numeric value representing the prior standard deviation for
#'        the control group. Required if \code{prior = 'N-Inv-Chisq'} and
#'        \code{design \%in\% c('controlled', 'external')}; otherwise \code{NULL}.
#'        Default is \code{NULL}.
#' @param mu_t A numeric value representing the true mean for the treatment group in the simulation.
#' @param mu_c A numeric value representing the true mean for the control group in the simulation.
#'        For uncontrolled design, this represents the historical control mean.
#'        Set to \code{NULL} if \code{design = 'uncontrolled'}.
#' @param sigma_t A positive numeric value representing the true standard deviation
#'        for the treatment group in the simulation.
#' @param sigma_c A positive numeric value representing the true standard deviation
#'        for the control group in the simulation. For uncontrolled design, this represents
#'        the historical control standard deviation.
#'        Set to \code{NULL} if \code{design = 'uncontrolled'}.
#' @param r A positive numeric value representing the variance scaling factor that allows
#'        the scale of hypothetical control to be different from treatment. Specifically,
#'        \code{sd.control = sqrt(r) * sd.treatment}. Required if \code{design = 'uncontrolled'}.
#'        When \code{r = 1}, the control and treatment have the same variance scale.
#' @param ne_t A positive integer representing the number of patients in the treatment group for
#'        the external data. Required if \code{design = 'external'} and external
#'        treatment data are available.
#' @param ne_c A positive integer representing the number of patients in the control group for
#'        the external data. Required if \code{design = 'external'} and external
#'        control data are available.
#' @param alpha0e_t A numeric value in (0, 1] representing the power prior scale parameter
#'        for the treatment group. Controls the degree of borrowing from external treatment data:
#'        0 = no borrowing, 1 = full borrowing. Required if \code{ne_t} is specified.
#' @param alpha0e_c A numeric value in (0, 1] representing the power prior scale parameter
#'        for the control group. Controls the degree of borrowing from external control data:
#'        0 = no borrowing, 1 = full borrowing. Required if \code{ne_c} is specified.
#' @param bar_ye_t A numeric value representing the sample mean of the external data
#'        for the treatment group. Required if \code{ne_t} is specified.
#' @param bar_ye_c A numeric value representing the sample mean of the external data
#'        for the control group. Required if \code{ne_c} is specified.
#' @param se_t A positive numeric value representing the sample standard deviation
#'        of the external data for the treatment group. Required if \code{ne_t} is specified.
#' @param se_c A positive numeric value representing the sample standard deviation
#'        of the external data for the control group. Required if \code{ne_c} is specified.
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
#'         only mu_t is included (not mu_c).
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
#' For external designs, power priors are incorporated using exact conjugate representation:
#' \itemize{
#'   \item Power priors for normal data are mathematically equivalent to Normal-Inverse-Chi-squared distributions
#'   \item This enables closed-form computation without MCMC sampling
#'   \item Alpha parameters control the degree of borrowing (0 = no borrowing, 1 = full borrowing)
#' }
#'
#' \strong{Decision rules}:
#' \itemize{
#'   \item \strong{Go}: P(treatment effect > threshold) >= gamma_go
#'   \item \strong{NoGo}: P(treatment effect > threshold) <= gamma_nogo
#'   \item \strong{Gray}: gamma_nogo < P(treatment effect > threshold) < gamma_go
#'   \item \strong{Miss}: Both Go and NoGo criteria are met simultaneously (indicates
#'                   poorly chosen thresholds)
#' }
#'
#' \strong{Handling Miss probability}:
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
#'   \item \strong{Controlled design}: Two-arm randomized trial
#'   \item \strong{Uncontrolled design}: Single-arm trial with informative priors (historical control)
#'   \item \strong{External design}: Incorporating historical data through power priors
#' }
#'
#' @examples
#' # Example 1: Controlled design with vague prior and NI method
#' # (default: error_if_Miss = TRUE, Gray_inc_Miss = FALSE)
#' pbayesdecisionprob1cont(
#'   nsim = 100, prob = 'posterior', design = 'controlled', prior = 'vague', CalcMethod = 'NI',
#'   theta_TV = 1.5, theta_MAV = -0.5, theta_NULL = NULL,
#'   nMC = NULL, gamma_go = 0.7, gamma_nogo = 0.2,
#'   n_t = 15, n_c = 15, m_t = NULL, m_c = NULL,
#'   kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
#'   mu0_t = NULL, mu0_c = NULL, sigma0_t = NULL, sigma0_c = NULL,
#'   mu_t = 3, mu_c = 1, sigma_t = 1.2, sigma_c = 1.1,
#'   r = NULL, ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
#'   bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed = 2
#' )
#'
#' # Example 2: Uncontrolled design with informative prior
#' pbayesdecisionprob1cont(
#'   nsim = 100, prob = 'posterior', design = 'uncontrolled', prior = 'N-Inv-Chisq', CalcMethod = 'NI',
#'   theta_TV = 1.0, theta_MAV = 0.0, theta_NULL = NULL,
#'   nMC = NULL, gamma_go = 0.8, gamma_nogo = 0.2,
#'   n_t = 20, n_c = NULL, m_t = NULL, m_c = NULL,
#'   kappa0_t = 2, kappa0_c = NULL, nu0_t = 5, nu0_c = NULL,
#'   mu0_t = 3.0, mu0_c = 1.5, sigma0_t = 1.5, sigma0_c = NULL,
#'   mu_t = 3.5, mu_c = NULL, sigma_t = 1.3, sigma_c = NULL,
#'   r = 1, ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
#'   bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed = 3
#' )
#'
#' # Example 3: External design with control data using MM approximation
#' \dontrun{
#' pbayesdecisionprob1cont(
#'   nsim = 100, prob = 'posterior', design = 'external', prior = 'vague', CalcMethod = 'MM',
#'   theta_TV = 1.0, theta_MAV = 0.0, theta_NULL = NULL,
#'   nMC = NULL, gamma_go = 0.8, gamma_nogo = 0.2,
#'   n_t = 12, n_c = 12, m_t = NULL, m_c = NULL,
#'   kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
#'   mu0_t = NULL, mu0_c = NULL, sigma0_t = NULL, sigma0_c = NULL,
#'   mu_t = 2, mu_c = 0, sigma_t = 1, sigma_c = 1,
#'   r = NULL, ne_t = NULL, ne_c = 20, alpha0e_t = NULL, alpha0e_c = 0.5,
#'   bar_ye_t = NULL, bar_ye_c = 0, se_t = NULL, se_c = 1,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed = 4
#' )
#' }
#'
#' # Example 4: Controlled design with predictive probability
#' pbayesdecisionprob1cont(
#'   nsim = 100, prob = 'predictive', design = 'controlled', prior = 'N-Inv-Chisq', CalcMethod = 'NI',
#'   theta_TV = NULL, theta_MAV = NULL, theta_NULL = 2.0,
#'   nMC = NULL, gamma_go = 0.75, gamma_nogo = 0.35,
#'   n_t = 15, n_c = 15, m_t = 50, m_c = 50,
#'   kappa0_t = 3, kappa0_c = 3, nu0_t = 4, nu0_c = 4,
#'   mu0_t = 3.5, mu0_c = 1.5, sigma0_t = 1.5, sigma0_c = 1.5,
#'   mu_t = 3.2, mu_c = 1.3, sigma_t = 1.4, sigma_c = 1.2,
#'   r = NULL, ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
#'   bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed = 5
#' )
#'
#' # Example 5: Uncontrolled design with predictive probability
#' pbayesdecisionprob1cont(
#'   nsim = 100, prob = 'predictive', design = 'uncontrolled', prior = 'vague', CalcMethod = 'NI',
#'   theta_TV = NULL, theta_MAV = NULL, theta_NULL = 1.0,
#'   nMC = NULL, gamma_go = 0.75, gamma_nogo = 0.35,
#'   n_t = 20, n_c = NULL, m_t = 40, m_c = 40,
#'   kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
#'   mu0_t = NULL, mu0_c = 1.5, sigma0_t = NULL, sigma0_c = NULL,
#'   mu_t = 3, mu_c = NULL, sigma_t = 1.3, sigma_c = NULL,
#'   r = 1, ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
#'   bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed = 9
#' )
#'
#' # Example 6: External design with predictive probability using MC method
#' \dontrun{
#' pbayesdecisionprob1cont(
#'   nsim = 100, prob = 'predictive', design = 'external', prior = 'vague', CalcMethod = 'MC',
#'   theta_TV = NULL, theta_MAV = NULL, theta_NULL = 1.5,
#'   nMC = 5000, gamma_go = 0.7, gamma_nogo = 0.4,
#'   n_t = 12, n_c = 12, m_t = 30, m_c = 30,
#'   kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
#'   mu0_t = NULL, mu0_c = NULL, sigma0_t = NULL, sigma0_c = NULL,
#'   mu_t = 2.5, mu_c = 1.0, sigma_t = 1.3, sigma_c = 1.1,
#'   r = NULL, ne_t = 15, ne_c = 18, alpha0e_t = 0.6, alpha0e_c = 0.7,
#'   bar_ye_t = 2.3, bar_ye_c = 0.9, se_t = 1.2, se_c = 1.0,
#'   error_if_Miss = FALSE, Gray_inc_Miss = FALSE, seed = 6
#' )
#' }
#'
#' @importFrom stats rnorm
#' @export
pbayesdecisionprob1cont <- function(nsim, prob, design, prior, CalcMethod,
                                    theta_TV = NULL, theta_MAV = NULL, theta_NULL = NULL,
                                    nMC = NULL, gamma_go, gamma_nogo, n_t, n_c = NULL,
                                    m_t = NULL, m_c = NULL,
                                    kappa0_t = NULL, kappa0_c = NULL,
                                    nu0_t = NULL, nu0_c = NULL,
                                    mu0_t = NULL, mu0_c = NULL,
                                    sigma0_t = NULL, sigma0_c = NULL,
                                    mu_t, mu_c = NULL, sigma_t, sigma_c = NULL,
                                    r = NULL, ne_t = NULL, ne_c = NULL,
                                    alpha0e_t = NULL, alpha0e_c = NULL,
                                    bar_ye_t = NULL, bar_ye_c = NULL,
                                    se_t = NULL, se_c = NULL,
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

  # Validate n_t (always required)
  if (!is.numeric(n_t) || length(n_t) != 1L || is.na(n_t) ||
      n_t != floor(n_t) || n_t < 1L) {
    stop("'n_t' must be a single positive integer")
  }

  # Validate n_c (required only for controlled and external designs)
  if (design %in% c("controlled", "external")) {
    if (!is.numeric(n_c) || length(n_c) != 1L || is.na(n_c) ||
        n_c != floor(n_c) || n_c < 1L) {
      stop("'n_c' must be a single positive integer")
    }
  }

  if (!is.numeric(gamma_go) || length(gamma_go) != 1L || is.na(gamma_go) ||
      gamma_go <= 0 || gamma_go >= 1) {
    stop("'gamma_go' must be a single numeric value in (0, 1)")
  }

  if (!is.numeric(gamma_nogo) || length(gamma_nogo) != 1L || is.na(gamma_nogo) ||
      gamma_nogo <= 0 || gamma_nogo >= 1) {
    stop("'gamma_nogo' must be a single numeric value in (0, 1)")
  }

  if (gamma_nogo >= gamma_go) {
    stop("'gamma_nogo' must be strictly less than 'gamma_go'")
  }

  if (!is.numeric(sigma_t) || length(sigma_t) != 1L || is.na(sigma_t) || sigma_t <= 0) {
    stop("'sigma_t' must be a single positive numeric value")
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
    if (is.null(theta_TV) || is.null(theta_MAV)) {
      stop("'theta_TV' and 'theta_MAV' must be non-NULL when prob = 'posterior'")
    }
    if (!is.numeric(theta_TV) || length(theta_TV) != 1L || is.na(theta_TV)) {
      stop("'theta_TV' must be a single numeric value")
    }
    if (!is.numeric(theta_MAV) || length(theta_MAV) != 1L || is.na(theta_MAV)) {
      stop("'theta_MAV' must be a single numeric value")
    }
    if (theta_TV <= theta_MAV) {
      stop("'theta_TV' must be strictly greater than 'theta_MAV'")
    }
  } else {
    if (is.null(theta_NULL)) {
      stop("'theta_NULL' must be non-NULL when prob = 'predictive'")
    }
    if (!is.numeric(theta_NULL) || length(theta_NULL) != 1L || is.na(theta_NULL)) {
      stop("'theta_NULL' must be a single numeric value")
    }
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

  # Validate prior-specific parameters
  if (prior == "N-Inv-Chisq") {
    if (is.null(kappa0_t) || is.null(nu0_t) || is.null(mu0_t) || is.null(sigma0_t)) {
      stop("'kappa0_t', 'nu0_t', 'mu0_t', and 'sigma0_t' must be non-NULL when prior = 'N-Inv-Chisq'")
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
      if (!is.numeric(mu0_c) || length(mu0_c) != 1L || is.na(mu0_c)) {
        stop("'mu0_c' must be a single numeric value")
      }
      if (!is.numeric(sigma0_c) || length(sigma0_c) != 1L || is.na(sigma0_c) || sigma0_c <= 0) {
        stop("'sigma0_c' must be a single positive numeric value")
      }
    }
  }

  # Validate design-specific parameters
  if (design == "uncontrolled") {
    if (is.null(mu0_c)) {
      stop("'mu0_c' must be non-NULL when design = 'uncontrolled'")
    }
    if (!is.numeric(mu0_c) || length(mu0_c) != 1L || is.na(mu0_c)) {
      stop("'mu0_c' must be a single numeric value")
    }
    if (is.null(r)) {
      stop("'r' must be non-NULL when design = 'uncontrolled'")
    }
    if (!is.numeric(r) || length(r) != 1L || is.na(r) || r <= 0) {
      stop("'r' must be a single positive numeric value")
    }
  }

  if (design == "external") {
    has_ext1 <- !is.null(ne_t) && !is.null(alpha0e_t) && !is.null(bar_ye_t) && !is.null(se_t)
    has_ext2 <- !is.null(ne_c) && !is.null(alpha0e_c) && !is.null(bar_ye_c) && !is.null(se_c)
    if (!has_ext1 && !has_ext2) {
      stop(paste0("For design = 'external', at least one complete set of external data ",
                  "(ne_t + alpha0e_t + bar_ye_t + se_t, or ne_c + alpha0e_c + bar_ye_c + se_c) ",
                  "must be provided"))
    }
    if (!is.null(ne_t)) {
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
  }

  # --- Simulation setup ---
  set.seed(seed)

  # Shared across all mu_t scenarios: bar_y_t[scenario] = rowMeans(Z_t) + mu_t[scenario],
  # and s_t is independent of mu_t, so it is computed once and reused.
  Z_t <- matrix(rnorm(nsim * n_t, mean = 0, sd = sigma_t), nrow = nsim)
  Z_t_rowmean <- rowSums(Z_t) / n_t
  s_t <- sqrt(rowSums((Z_t - Z_t_rowmean) ^ 2) / (n_t - 1))

  if (design == 'controlled' | design == 'external') {
    # Generate nsim x n_c matrix of standardized residuals for the control group (mean=0, sd=sigma_c).
    # mu_c is fixed to a single value, so bar_y_c and s_c are computed once.
    Z_c <- matrix(rnorm(nsim * n_c, mean = 0, sd = sigma_c), nrow = nsim)
    Z_c_rowmean <- rowSums(Z_c) / n_c
    bar_y_c <- Z_c_rowmean + mu_c
    s_c     <- sqrt(rowSums((Z_c - Z_c_rowmean) ^ 2) / (n_c - 1))
  } else if (design == 'uncontrolled') {
    # For uncontrolled design: no control group data generated
    bar_y_c <- NULL
    s_c     <- NULL
  }

  # External data are fixed historical sample statistics (no generation needed)

  # Ensure mu_t is a numeric vector to support multiple scenarios
  mu_t <- as.numeric(mu_t)
  n_scenarios <- length(mu_t)

  # Set threshold values based on probability type
  if (prob == 'posterior') {
    theta0_Go   <- theta_TV
    theta0_NoGo <- theta_MAV
  } else {
    theta0_Go   <- theta_NULL
    theta0_NoGo <- theta_NULL
  }

  # Expand bar_y_t and s_t across all mu_t scenarios into a single vector of
  # length nsim * n_scenarios. For bar_y_t, shift each scenario's block by mu_t[i].
  # s_t is independent of mu_t so it is simply replicated.
  bar_y_t_all   <- rep(Z_t_rowmean, times = n_scenarios) + rep(mu_t, each = nsim)
  s_t_all        <- rep(s_t, times = n_scenarios)

  if (design == 'controlled' | design == 'external') {
    bar_y_c_all <- rep(bar_y_c, times = n_scenarios)
    s_c_all      <- rep(s_c,     times = n_scenarios)
  } else {
    bar_y_c_all <- NULL
    s_c_all      <- NULL
  }

  # Call pbayespostpred1cont once with the expanded vectors.
  # Returns a vector of length nsim * n_scenarios.
  gPost_Go_all <- pbayespostpred1cont(
    prob = prob, design = design, prior = prior, CalcMethod = CalcMethod,
    theta0 = theta0_Go, nMC = nMC, n_t = n_t, n_c = n_c, m_t = m_t, m_c = m_c,
    kappa0_t = kappa0_t, kappa0_c = kappa0_c, nu0_t = nu0_t, nu0_c = nu0_c,
    mu0_t = mu0_t, mu0_c = mu0_c, sigma0_t = sigma0_t, sigma0_c = sigma0_c,
    bar_y_t = bar_y_t_all, bar_y_c = bar_y_c_all, s_t = s_t_all, s_c = s_c_all,
    r = r, ne_t = ne_t, ne_c = ne_c, alpha0e_t = alpha0e_t, alpha0e_c = alpha0e_c,
    bar_ye_t = bar_ye_t, bar_ye_c = bar_ye_c, se_t = se_t, se_c = se_c,
    lower.tail = FALSE
  )

  gPost_NoGo_all <- pbayespostpred1cont(
    prob = prob, design = design, prior = prior, CalcMethod = CalcMethod,
    theta0 = theta0_NoGo, nMC = nMC, n_t = n_t, n_c = n_c, m_t = m_t, m_c = m_c,
    kappa0_t = kappa0_t, kappa0_c = kappa0_c, nu0_t = nu0_t, nu0_c = nu0_c,
    mu0_t = mu0_t, mu0_c = mu0_c, sigma0_t = sigma0_t, sigma0_c = sigma0_c,
    bar_y_t = bar_y_t_all, bar_y_c = bar_y_c_all, s_t = s_t_all, s_c = s_c_all,
    r = r, ne_t = ne_t, ne_c = ne_c, alpha0e_t = alpha0e_t, alpha0e_c = alpha0e_c,
    bar_ye_t = bar_ye_t, bar_ye_c = bar_ye_c, se_t = se_t, se_c = se_c,
    lower.tail = TRUE
  )

  # Reshape results into nsim x n_scenarios matrices and compute decision
  # probabilities per scenario using colMeans (no R-level loop needed).
  probs_Go_mat   <- matrix((gPost_Go_all >= gamma_go) & (gPost_NoGo_all <  gamma_nogo),
                           nrow = nsim)
  probs_NoGo_mat <- matrix((gPost_Go_all <  gamma_go) & (gPost_NoGo_all >= gamma_nogo),
                           nrow = nsim)
  probs_Miss_mat <- matrix((gPost_Go_all >= gamma_go) & (gPost_NoGo_all >= gamma_nogo),
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
    # For uncontrolled design, only include mu_t
    results <- data.frame(
      mu_t, Go = GoNogoProb[, 1], Gray = GrayProb, NoGo = GoNogoProb[, 2]
    )
  } else {
    # For controlled and external designs, include both mu_t and mu_c
    results <- data.frame(
      mu_t, mu_c, Go = GoNogoProb[, 1], Gray = GrayProb, NoGo = GoNogoProb[, 2]
    )
  }

  # Add Miss column when error_if_Miss is FALSE and Gray_inc_Miss is FALSE
  if (!error_if_Miss) {
    if (!Gray_inc_Miss) {
      results$Miss <- GoNogoProb[, 3]
    }
  }

  # Address floating point error (apply only to probability columns, not mu_t/mu_c)
  prob_cols <- names(results)[!names(results) %in% c('mu_t', 'mu_c')]
  results[prob_cols] <- lapply(results[prob_cols], function(col) {
    ifelse(col < .Machine$double.eps ^ 0.25, 0, col)
  })

  # Attach metadata as attributes for use in print()
  attr(results, 'prob')           <- prob
  attr(results, 'design')         <- design
  attr(results, 'prior')          <- prior
  attr(results, 'CalcMethod')     <- CalcMethod
  attr(results, 'nsim')           <- nsim
  attr(results, 'nMC')            <- nMC
  attr(results, 'theta_TV')       <- theta_TV
  attr(results, 'theta_MAV')      <- theta_MAV
  attr(results, 'theta_NULL')     <- theta_NULL
  attr(results, 'gamma_go')       <- gamma_go
  attr(results, 'gamma_nogo')     <- gamma_nogo
  attr(results, 'n_t')            <- n_t
  attr(results, 'n_c')            <- n_c
  attr(results, 'sigma_t')        <- sigma_t
  attr(results, 'sigma_c')        <- sigma_c
  attr(results, 'r')              <- r
  attr(results, 'm_t')            <- m_t
  attr(results, 'm_c')            <- m_c
  attr(results, 'kappa0_t')       <- kappa0_t
  attr(results, 'kappa0_c')       <- kappa0_c
  attr(results, 'nu0_t')          <- nu0_t
  attr(results, 'nu0_c')          <- nu0_c
  attr(results, 'mu0_t')          <- mu0_t
  attr(results, 'mu0_c')          <- mu0_c
  attr(results, 'sigma0_t')       <- sigma0_t
  attr(results, 'sigma0_c')       <- sigma0_c
  attr(results, 'ne_t')           <- ne_t
  attr(results, 'ne_c')           <- ne_c
  attr(results, 'alpha0e_t')       <- alpha0e_t
  attr(results, 'alpha0e_c')       <- alpha0e_c
  attr(results, 'bar_ye_t')       <- bar_ye_t
  attr(results, 'bar_ye_c')       <- bar_ye_c
  attr(results, 'se_t')           <- se_t
  attr(results, 'se_c')           <- se_c
  attr(results, 'error_if_Miss')  <- error_if_Miss
  attr(results, 'Gray_inc_Miss')  <- Gray_inc_Miss
  attr(results, 'seed')           <- seed

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
  prob           <- attr(x, 'prob')
  design         <- attr(x, 'design')
  prior          <- attr(x, 'prior')
  CalcMethod     <- attr(x, 'CalcMethod')
  nsim           <- attr(x, 'nsim')
  nMC            <- attr(x, 'nMC')
  gamma_go       <- attr(x, 'gamma_go')
  gamma_nogo     <- attr(x, 'gamma_nogo')
  n_t            <- attr(x, 'n_t')
  n_c            <- attr(x, 'n_c')
  sigma_t        <- attr(x, 'sigma_t')
  sigma_c        <- attr(x, 'sigma_c')
  r              <- attr(x, 'r')
  m_t            <- attr(x, 'm_t')
  m_c            <- attr(x, 'm_c')
  kappa0_t       <- attr(x, 'kappa0_t')
  kappa0_c       <- attr(x, 'kappa0_c')
  nu0_t          <- attr(x, 'nu0_t')
  nu0_c          <- attr(x, 'nu0_c')
  mu0_t          <- attr(x, 'mu0_t')
  mu0_c          <- attr(x, 'mu0_c')
  sigma0_t       <- attr(x, 'sigma0_t')
  sigma0_c       <- attr(x, 'sigma0_c')
  ne_t           <- attr(x, 'ne_t')
  ne_c           <- attr(x, 'ne_c')
  alpha0e_t       <- attr(x, 'alpha0e_t')
  alpha0e_c       <- attr(x, 'alpha0e_c')
  bar_ye_t       <- attr(x, 'bar_ye_t')
  bar_ye_c       <- attr(x, 'bar_ye_c')
  se_t           <- attr(x, 'se_t')
  se_c           <- attr(x, 'se_c')
  error_if_Miss  <- attr(x, 'error_if_Miss')
  Gray_inc_Miss  <- attr(x, 'Gray_inc_Miss')
  seed           <- attr(x, 'seed')

  # Build threshold string based on probability type
  if(prob == 'posterior') {
    theta_str <- sprintf('TV = %s, MAV = %s',
                         fmt(attr(x, 'theta_TV')), fmt(attr(x, 'theta_MAV')))
  } else {
    theta_str <- sprintf('NULL = %s', fmt(attr(x, 'theta_NULL')))
  }

  # Print header
  cat('Go/NoGo/Gray Decision Probabilities (Single Continuous Endpoint)\n')
  cat(strrep('-', 60), '\n')
  cat(sprintf('  Probability type    : %s\n',   prob))
  cat(sprintf('  Design              : %s\n',   design))
  cat(sprintf('  Prior               : %s\n',   prior))
  cat(sprintf('  Calc method         : %s\n',   CalcMethod))
  cat(sprintf('  Simulations         : nsim = %s\n', fmt(nsim)))
  if (!is.null(nMC)) {
    cat(sprintf('  MC draws            : nMC = %s\n', fmt(nMC)))
  }
  cat(sprintf('  Threshold(s)        : %s\n',   theta_str))
  cat(sprintf('  Go  threshold       : gamma_go = %s\n', fmt(gamma_go)))
  cat(sprintf('  NoGo threshold      : gamma_nogo = %s\n', fmt(gamma_nogo)))
  cat(sprintf('  Sample size         : n_t = %s, n_c = %s\n', fmt(n_t), fmt(n_c)))
  cat(sprintf('  True SD             : sigma_t = %s, sigma_c = %s\n', fmt(sigma_t), fmt(sigma_c)))
  if (design == 'uncontrolled') {
    cat(sprintf('  Variance ratio      : r = %s\n', fmt(r)))
  }
  if (!is.null(m_t) || !is.null(m_c)) {
    cat(sprintf('  Future size         : m_t = %s, m_c = %s\n', fmt(m_t), fmt(m_c)))
  }
  if (prior == 'N-Inv-Chisq') {
    cat(sprintf('  Prior (traetment)   : kappa0_t = %s, nu0_t = %s, mu0_t = %s, sigma0_t = %s\n',
                fmt(kappa0_t), fmt(nu0_t), fmt(mu0_t), fmt(sigma0_t)))
    cat(sprintf('  Prior (control)     : kappa0_c = %s, nu0_c = %s, mu0_c = %s, sigma0_c = %s\n',
                fmt(kappa0_c), fmt(nu0_c), fmt(mu0_c), fmt(sigma0_c)))
  }
  if (design == 'external') {
    cat(sprintf('  External (traetment): ne_t = %s, alpha0e_t = %s, bar_ye_t = %s, se_t = %s\n',
                fmt(ne_t), fmt(alpha0e_t), fmt(bar_ye_t), fmt(se_t)))
    cat(sprintf('  External (control)  : ne_c = %s, alpha0e_c = %s, bar_ye_c = %s, se_c = %s\n',
                fmt(ne_c), fmt(alpha0e_c), fmt(bar_ye_c), fmt(se_c)))
  }
  cat(sprintf('  error_if_Miss       : %s\n', fmt(error_if_Miss)))
  cat(sprintf('  Gray_inc_Miss       : %s\n', fmt(Gray_inc_Miss)))
  cat(sprintf('  Seed                : %s\n', fmt(seed)))
  cat(strrep('-', 60), '\n')

  # Print results table
  df_print        <- as.data.frame(x)
  prob_cols       <- names(df_print)[!names(df_print) %in% c('mu_t', 'mu_c')]
  df_print[prob_cols] <- lapply(df_print[prob_cols],
                                function(col) round(col, digits))
  print.data.frame(df_print, row.names = FALSE)

  invisible(x)
}
