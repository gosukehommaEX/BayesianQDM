#' Posterior Probability Calculation using MCMC Sampling for Continuous Endpoints
#'
#' This function calculates the Bayesian posterior probability P(μt - μc > q)
#' for continuous endpoints using MCMC sampling to obtain the posterior distribution
#' of the treatment effect. The function optionally supports incorporation of external
#' (historical) data through power prior methodology.
#'
#' @param nMCMCsample A positive integer representing the number of Markov Chain Monte Carlo iterations.
#' @param q A numeric value representing the quantile threshold.
#' @param mu.n1 A numeric vector representing the mean values of the treatment group in current trial.
#' @param mu.n2 A numeric vector representing the mean values of the control group in current trial.
#' @param sd.n1 A positive numeric vector representing the standard deviations of the treatment group in current trial.
#' @param sd.n2 A positive numeric vector representing the standard deviations of the control group in current trial.
#' @param n1 A positive integer representing the sample size of the treatment group in current trial.
#' @param n2 A positive integer representing the sample size of the control group in current trial.
#' @param ne1 A positive integer representing the sample size of the treatment group in external trial (can be NULL if no external treatment data).
#' @param ne2 A positive integer representing the sample size of the control group in external trial (can be NULL if no external control data).
#' @param alpha01 A positive numeric value between 0 and 1 representing the power prior scaling parameter for external treatment data (can be NULL if no external treatment data).
#' @param alpha02 A positive numeric value between 0 and 1 representing the power prior scaling parameter for external control data (can be NULL if no external control data).
#'
#' @return A numeric vector representing P(μt - μc > q), the posterior probability
#' that the treatment effect (difference between treatment and control group means)
#' exceeds the quantile q, based on MCMC sampling of the posterior distribution.
#' The length of the output matches the maximum length of the input vectors.
#'
#' @details
#' This function performs Bayesian posterior inference using MCMC sampling to compute
#' the probability that the treatment effect exceeds a specified threshold. The core
#' methodology involves:
#' \itemize{
#'   \item MCMC sampling to obtain the posterior distribution of μt and μc
#'   \item Computing the posterior distribution of the treatment effect (μt - μc)
#'   \item Calculating the probability that this difference exceeds the threshold q
#' }
#'
#' The function supports standard two-group Bayesian analysis when no external data
#' is provided (ne1=NULL, ne2=NULL, alpha01=NULL, alpha02=NULL). In this case, the
#' analysis is based solely on the current trial data using non-informative priors.
#'
#' When external data is available, the function can incorporate historical information
#' using power prior methodology:
#' - alpha = 1: Full borrowing (external data weighted equally to current data)
#' - alpha = 0: No borrowing (external data ignored)
#' - 0 < alpha < 1: Partial borrowing (external data down-weighted)
#'
#' The resulting posterior distribution is typically non-conjugate and requires
#' MCMC sampling for computation. The function uses the bayesDP package's
#' bdpnormal function for robust MCMC-based Bayesian estimation.
#'
#' All main distribution parameters (mu.n1, mu.n2, sd.n1, sd.n2) support vectorized inputs.
#' When vectors of different lengths are provided, they are recycled to match the longest vector.
#' The function returns a vector of the same length containing the corresponding probabilities.
#'
#' @examples
#' \dontrun{
#' # Standard Bayesian analysis without external data - single values
#' pMCMCdiff(nMCMCsample = 1e+4, q = 4, mu.n1 = 5, mu.n2 = 0, sd.n1 = 1, sd.n2 = 1,
#'           n1 = 12, n2 = 12, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL)
#'
#' # Vectorized analysis with multiple parameter sets
#' pMCMCdiff(nMCMCsample = 1e+4, q = 4,
#'           mu.n1 = c(5, 6, 7), mu.n2 = c(0, 1, 2),
#'           sd.n1 = c(1, 1.2, 1.5), sd.n2 = c(1, 1.1, 1.3),
#'           n1 = 12, n2 = 12, ne1 = NULL, ne2 = NULL,
#'           alpha01 = NULL, alpha02 = NULL)
#'
#' # With external control data using power prior - vectorized
#' pMCMCdiff(nMCMCsample = 1e+4, q = 4,
#'           mu.n1 = c(5, 6), mu.n2 = c(0, 1),
#'           sd.n1 = c(1, 1.2), sd.n2 = c(1, 1.1),
#'           n1 = 12, n2 = 12, ne1 = NULL, ne2 = 24,
#'           alpha01 = NULL, alpha02 = 0.5)
#' }
#'
#' @importFrom bayesDP bdpnormal
#' @export
pMCMCdiff <- function(nMCMCsample, q, mu.n1, mu.n2, sd.n1, sd.n2, n1, n2, ne1, ne2, alpha01, alpha02) {
  # Check if bayesDP package is available
  if (!requireNamespace("bayesDP", quietly = TRUE)) {
    stop("bayesDP package is required for this function. Please install it from CRAN.")
  }

  # Determine the number of results to compute
  n <- max(length(mu.n1), length(mu.n2), length(sd.n1), length(sd.n2))

  # Ensure all main parameters have the same length
  mu.n1 <- rep(mu.n1, length.out = n)
  mu.n2 <- rep(mu.n2, length.out = n)
  sd.n1 <- rep(sd.n1, length.out = n)
  sd.n2 <- rep(sd.n2, length.out = n)

  # Calculate P(μt - μc > q) for each parameter set
  results <- sapply(seq(n), function(i) {

    # Check parameter sets for external data availability
    has_external_treatment <- (!is.null(ne1) & !is.null(alpha01))
    has_external_control <- (!is.null(ne2) & !is.null(alpha02))

    if(!has_external_treatment & !has_external_control) {
      # No external data case - use standard bdpnormal without historical data
      bdp_args <- list(
        mu_t = mu.n1[i],
        sigma_t = sd.n1[i],
        N_t = n1,
        mu_c = mu.n2[i],
        sigma_c = sd.n2[i],
        N_c = n2,
        method = "fixed",
        number_mcmc = nMCMCsample,
        compare = TRUE
      )
    } else {
      # At least one external data source is provided
      # Prepare arguments for bdpnormal function with external data
      bdp_args <- list(
        mu_t = mu.n1[i],
        sigma_t = sd.n1[i],
        N_t = n1,
        mu_c = mu.n2[i],
        sigma_c = sd.n2[i],
        N_c = n2,
        fix_alpha = TRUE,
        alpha_max = c(1, 1),  # Default values, will be updated based on available data
        method = "fixed",
        number_mcmc = nMCMCsample,
        compare = TRUE
      )
    }

    # Handle different external data scenarios
    if(has_external_treatment & has_external_control) {
      ## Both external treatment and control data are available
      if(!is.null(alpha01) & !is.null(alpha02) & alpha01 == 0 & alpha02 == 0) {
        # When both alphas are 0, completely ignore external data
        # Keep bdp_args as no-external-data version
        bdp_args <- bdp_args[!names(bdp_args) %in% c("fix_alpha", "alpha_max")]
      } else {
        # Set external data parameters
        bdp_args$mu0_t <- mu.n1[i]
        bdp_args$sigma0_t <- sd.n1[i]
        bdp_args$N0_t <- ne1
        bdp_args$mu0_c <- mu.n2[i]
        bdp_args$sigma0_c <- sd.n2[i]
        bdp_args$N0_c <- ne2
        bdp_args$alpha_max <- c(alpha01, alpha02)
      }

    } else if(has_external_treatment & !has_external_control) {
      ## External treatment data only
      if(!is.null(alpha01) & alpha01 == 0) {
        # When alpha01 is 0, ignore external treatment data
        bdp_args <- bdp_args[!names(bdp_args) %in% c("fix_alpha", "alpha_max")]
      } else {
        # Set external treatment data parameters
        bdp_args$mu0_t <- mu.n1[i]
        bdp_args$sigma0_t <- sd.n1[i]
        bdp_args$N0_t <- ne1
        bdp_args$alpha_max <- c(alpha01, 1)  # Only treatment side uses power prior
      }

    } else if(!has_external_treatment & has_external_control) {
      ## External control data only
      if(!is.null(alpha02) & alpha02 == 0) {
        # When alpha02 is 0, ignore external control data
        bdp_args <- bdp_args[!names(bdp_args) %in% c("fix_alpha", "alpha_max")]
      } else {
        # Set external control data parameters
        bdp_args$mu0_c <- mu.n2[i]
        bdp_args$sigma0_c <- sd.n2[i]
        bdp_args$N0_c <- ne2
        bdp_args$alpha_max <- c(1, alpha02)  # Only control side uses power prior
      }
    }

    # Fit Bayesian model using bayesDP
    result <- do.call(bayesDP::bdpnormal, bdp_args)

    # Extract posterior samples of the treatment effect (mu_t - mu_c)
    if (is.null(result$final$posterior)) {
      stop("Error: Could not obtain posterior samples from bdpnormal. Check input parameters.")
    }

    # Calculate the probability that the treatment effect exceeds the threshold q
    # This represents P(μt - μc > q | data, external_data, alpha)
    prob_result <- mean(result$final$posterior > q)

    return(prob_result)
  })

  return(results)
}
