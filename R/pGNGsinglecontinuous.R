#' Calculate Go, NoGo, and Gray Probabilities for a Clinical Trial with a Single Continuous Endpoint
#' Under the Bayesian Framework Using Two Metrics
#'
#' This function calculates the Go, NoGo, and Gray probabilities for continuous outcome
#' clinical trials using the Bayesian framework. The function evaluates operating
#' characteristics by computing the probability of making each type of decision
#' (Go, NoGo, or Gray) across different true mean values through Monte Carlo simulation.
#' The function supports controlled, uncontrolled, and external control designs.
#'
#' @param nsim A positive integer representing the number of Monte Carlo simulations
#'        for evaluating operating characteristics.
#' @param prob A character string specifying the type of probability to use for
#'        decision-making. Options are \code{'posterior'} for posterior probability
#'        or \code{'predictive'} for posterior predictive probability.
#' @param design A character string specifying the type of trial design. Options are
#'        \code{'controlled'} for randomized controlled trials, \code{'uncontrolled'}
#'        for single-arm studies with informative priors, or \code{'external'} for
#'        designs incorporating external data through power priors.
#' @param prior A character string specifying the prior distribution to use.
#'        Options are \code{'vague'} for a vague prior or \code{'N-Inv-Chisq'} for
#'        a Normal-Inverse-Chi-squared prior.
#' @param CalcMethod A character string specifying the calculation method for
#'        computing probabilities. Options are \code{'NI'} for numerical integration,
#'        \code{'MC'} for Monte Carlo simulation, or \code{'WS'} for the
#'        Welch-Satterthwaite approximation.
#' @param theta.TV A numeric value representing the target value threshold for
#'        calculating Go probability when \code{prob = 'posterior'}. This represents
#'        the minimum clinically meaningful treatment effect.
#' @param theta.MAV A numeric value representing the minimum acceptable value threshold
#'        for calculating NoGo probability when \code{prob = 'posterior'}. This represents
#'        the minimum effect size below which the treatment is not worth pursuing.
#' @param theta.NULL A numeric value representing the null hypothesis value for
#'        calculating Go/NoGo probabilities when \code{prob = 'predictive'}.
#' @param nMC A positive integer representing the number of Monte Carlo samples for
#'        probability calculation when \code{CalcMethod = 'MC'} (required if
#'        \code{CalcMethod = 'MC'}).
#' @param gamma1 A numeric value in (0, 1) representing the threshold for Go decision.
#'        If P(treatment effect > threshold | data) ≥ gamma1, the decision is Go.
#' @param gamma2 A numeric value in (0, 1) representing the threshold for NoGo decision.
#'        If P(treatment effect > threshold | data) ≤ gamma2, the decision is NoGo.
#'        Must satisfy gamma2 < gamma1.
#' @param n1 A positive integer representing the number of patients in group 1
#'        (treatment) for the proof-of-concept (PoC) trial.
#' @param n2 A positive integer representing the number of patients in group 2
#'        (control) for the PoC trial. For uncontrolled design, this represents
#'        the effective sample size of the historical control (encoded in the prior).
#' @param m1 A positive integer representing the number of patients in group 1 for
#'        the future trial (required if \code{prob = 'predictive'}).
#' @param m2 A positive integer representing the number of patients in group 2 for
#'        the future trial (required if \code{prob = 'predictive'}).
#' @param kappa01 A positive numeric value representing the prior sample size parameter
#'        for group 1 when \code{prior = 'N-Inv-Chisq'}.
#' @param kappa02 A positive numeric value representing the prior sample size parameter
#'        for group 2 when \code{prior = 'N-Inv-Chisq'}.
#' @param nu01 A positive numeric value representing the prior degrees of freedom
#'        for group 1 when \code{prior = 'N-Inv-Chisq'}.
#' @param nu02 A positive numeric value representing the prior degrees of freedom
#'        for group 2 when \code{prior = 'N-Inv-Chisq'}.
#' @param mu01 A numeric value representing the prior mean for group 1 when
#'        \code{prior = 'N-Inv-Chisq'}.
#' @param mu02 A numeric value representing the prior mean for group 2 when
#'        \code{prior = 'N-Inv-Chisq'}.
#' @param sigma01 A positive numeric value representing the prior scale parameter
#'        for group 1 when \code{prior = 'N-Inv-Chisq'}.
#' @param sigma02 A positive numeric value representing the prior scale parameter
#'        for group 2 when \code{prior = 'N-Inv-Chisq'}.
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
#' For external control designs, power priors are incorporated using exact conjugate representation:
#' \itemize{
#'   \item Power priors for normal data are mathematically equivalent to Normal-Inverse-Chi-squared distributions
#'   \item This enables closed-form computation without MCMC sampling
#'   \item Alpha parameters control the degree of borrowing (0 = no borrowing, 1 = full borrowing)
#' }
#'
#' **Decision rules**:
#' \itemize{
#'   \item **Go**: P(treatment effect > threshold) ≥ γ₁
#'   \item **NoGo**: P(treatment effect > threshold) ≤ γ₂
#'   \item **Gray**: γ₂ < P(treatment effect > threshold) < γ₁
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
#' pGNGsinglecontinuous(
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
#' pGNGsinglecontinuous(
#'   nsim = 100, prob = 'posterior', design = 'uncontrolled', prior = 'N-Inv-Chisq', CalcMethod = 'NI',
#'   theta.TV = 1.0, theta.MAV = 0.0, theta.NULL = NULL,
#'   nMC = NULL, gamma1 = 0.8, gamma2 = 0.2,
#'   n1 = 20, n2 = 20, m1 = NULL, m2 = NULL,
#'   kappa01 = 2, kappa02 = 20, nu01 = 5, nu02 = 20,
#'   mu01 = 3.0, mu02 = 1.5, sigma01 = 1.5, sigma02 = 1.2,
#'   mu1 = 3.5, mu2 = 1.5, sigma1 = 1.3, sigma2 = 1.2,
#'   r = 1, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
#'   bar.ye1 = NULL, bar.ye2 = NULL, se1 = NULL, se2 = NULL,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed = 3
#' )
#'
#' # Example 3: External design with control data using WS approximation
#' \dontrun{
#' pGNGsinglecontinuous(
#'   nsim = 100, prob = 'posterior', design = 'external', prior = 'vague', CalcMethod = 'WS',
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
#' pGNGsinglecontinuous(
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
#' pGNGsinglecontinuous(
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
#' pGNGsinglecontinuous(
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
#' pGNGsinglecontinuous(
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
pGNGsinglecontinuous <- function(nsim, prob, design, prior, CalcMethod, theta.TV, theta.MAV, theta.NULL,
                                 nMC = NULL, gamma1, gamma2, n1, n2, m1, m2, kappa01, kappa02, nu01, nu02,
                                 mu01, mu02, sigma01, sigma02, mu1, mu2, sigma1, sigma2,
                                 r = NULL, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
                                 bar.ye1 = NULL, bar.ye2 = NULL, se1 = NULL, se2 = NULL,
                                 error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed) {

  # Set seed for reproducibility
  set.seed(seed)

  # Check parameter sets
  if (!prob %in% c("posterior", "predictive")) {
    stop("prob must be either 'posterior' or 'predictive'")
  }

  if (!design %in% c("controlled", "uncontrolled", "external")) {
    stop("design must be 'controlled', 'uncontrolled', or 'external'")
  }

  if (!prior %in% c("vague", "N-Inv-Chisq")) {
    stop("prior must be either 'vague' or 'N-Inv-Chisq'")
  }

  if (!CalcMethod %in% c("NI", "MC", "WS")) {
    stop("CalcMethod must be 'NI', 'MC', or 'WS'")
  }

  # Validate required parameters for each method
  if (CalcMethod == "MC" && is.null(nMC)) {
    stop("nMC must be specified for Monte Carlo method")
  }

  if (prob == "predictive" && (missing(m1) || missing(m2))) {
    stop("m1 and m2 must be specified for predictive probability")
  }

  # Validate uncontrolled design parameters
  if (design == "uncontrolled" && is.null(r)) {
    stop("For uncontrolled design, r (variance scaling factor) must be specified")
  }

  # Validate external design parameters
  if (design == "external") {
    if (is.null(ne1) && is.null(ne2)) {
      stop("For external design, at least one of ne1 or ne2 must be specified")
    }
    if (!is.null(ne1) && (is.null(alpha01) || is.null(bar.ye1) || is.null(se1))) {
      stop("For external treatment data, alpha01, bar.ye1, and se1 must be specified")
    }
    if (!is.null(ne2) && (is.null(alpha02) || is.null(bar.ye2) || is.null(se2))) {
      stop("For external control data, alpha02, bar.ye2, and se2 must be specified")
    }
  }

  # Set seed number for reproducible results
  set.seed(seed)

  # Generate random numbers of outcomes for group 1 in PoC study
  y1 <- matrix(rnorm(nsim * n1, mu1, sigma1), nrow = nsim)

  # Calculate sample mean for group 1
  bar.y1 <- rowSums(y1) / n1

  # Calculate sample standard deviation for group 1
  s1 <- sqrt(rowSums((y1 - bar.y1) ^ 2) / (n1 - 1))

  if(design == 'controlled' | design == 'external') {
    # Generate random numbers of outcomes for group 2 in PoC study
    y2 <- matrix(rnorm(nsim * n2, mu2, sigma2), nrow = nsim)

    # Calculate sample mean for group 2
    bar.y2 <- rowSums(y2) / n2

    # Calculate sample standard deviation for group 2
    s2 <- sqrt(rowSums((y2 - bar.y2) ^ 2) / (n2 - 1))
  } else if(design == 'uncontrolled') {
    # For uncontrolled design, group 2 data comes from informative prior (historical control)
    # No data generation needed - historical control encoded in prior parameters
    bar.y2 <- NULL
    s2 <- NULL
  }

  # External data are provided as fixed historical sample statistics
  # No data generation needed - these are already observed historical values

  # Initialize vectors to store Go and NoGo probabilities
  if(prob == 'posterior') {
    # Define thresholds: theta.TV for Go, theta.MAV for NoGo
    theta_values <- c(theta.TV, theta.MAV)
  } else {
    # For predictive probability, use theta.NULL (single threshold)
    theta_values <- c(theta.NULL, theta.NULL)
  }

  # Calculate posterior/posterior predictive probabilities for each threshold
  gPost <- lapply(seq_along(theta_values), function(i) {
    pPPsinglecontinuous(
      prob = prob, design = design, prior = prior, CalcMethod = CalcMethod,
      theta0 = theta_values[i], nMC = nMC, n1 = n1, n2 = n2, m1 = m1, m2 = m2,
      kappa01 = kappa01, kappa02 = kappa02, nu01 = nu01, nu02 = nu02,
      mu01 = mu01, mu02 = mu02, sigma01 = sigma01, sigma02 = sigma02,
      bar.y1 = bar.y1, bar.y2 = bar.y2, s1 = s1, s2 = s2,
      r = r, ne1 = ne1, ne2 = ne2, alpha01 = alpha01, alpha02 = alpha02,
      bar.ye1 = bar.ye1, bar.ye2 = bar.ye2, se1 = se1, se2 = se2, lower.tail = c(FALSE, TRUE)[i]
    )
  })

  # Calculate Go, NoGo and Miss probabilities based on decision criteria
  GoNogoProb <- sapply(seq(3), function(j) {
    # Create indicator matrix for Go (j=1), NoGo (j=2), or Miss (j=3) decisions
    if(j == 1) {
      # Go decision: high probability for Go threshold AND low probability for NoGo threshold
      I <- as.numeric((gPost[[1]] >= gamma1) & (gPost[[2]] < gamma2))
    } else if(j == 2) {
      # NoGo decision: low probability for Go threshold AND high probability for NoGo threshold
      I <- as.numeric((gPost[[1]] < gamma1) & (gPost[[2]] >= gamma2))
    } else {
      # Miss: both Go and NoGo criteria met simultaneously (should be rare)
      I <- as.numeric((gPost[[1]] >= gamma1) & (gPost[[2]] >= gamma2))
    }
    sum(I) / nsim
  })

  # Check for positive Miss probabilities (indicates inappropriate thresholds)
  if(error_if_Miss) {
    if(sum(GoNogoProb[3]) > 0) {
      stop('Because positive Miss probability(s) is obtained, re-consider appropriate thresholds')
    }
  }

  # Calculate Miss probability (both Go and NoGo criteria met simultaneously)
  Miss <- GoNogoProb[3]

  # Calculate Gray probability (complement of Go and NoGo)
  if(Gray_inc_Miss) {
    # Include Miss in Gray probability
    GrayProb <- 1 - sum(GoNogoProb[-3])
  } else {
    # Exclude Miss from Gray probability
    GrayProb <- 1 - sum(GoNogoProb)
  }

  # Create results data frame
  if(design == 'uncontrolled') {
    # For uncontrolled design, only include mu1
    results <- data.frame(
      mu1, Go = GoNogoProb[1], Gray = GrayProb, NoGo = GoNogoProb[2]
    )
  } else {
    # For controlled and external designs, include both mu1 and mu2
    results <- data.frame(
      mu1, mu2, Go = GoNogoProb[1], Gray = GrayProb, NoGo = GoNogoProb[2]
    )
  }

  # Add Miss column when error_if_Miss is FALSE and Gray_inc_Miss is FALSE
  if(!error_if_Miss) {
    if(!Gray_inc_Miss) {
      results$Miss <- Miss
    }
  }

  # Address floating point error
  results[results < .Machine$double.eps ^ 0.25] <- 0

  return(results)
}
