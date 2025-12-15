#' Calculate Go, NoGo, and Gray Probabilities for a Clinical Trial with a Single Binary Endpoint
#' Under the Bayesian Framework Using Two Metrics
#'
#' This function calculates the Go, NoGo, and Gray probabilities for binary outcome
#' clinical trials using the Bayesian framework. The function evaluates operating
#' characteristics by computing the probability of making each type of decision
#' (Go, NoGo, or Gray) across different true response rates. The function supports
#' controlled, uncontrolled, and external control designs.
#'
#' @param prob A character string specifying the type of probability to use for
#'        decision-making. Options are \code{'posterior'} for posterior probability
#'        or \code{'predictive'} for posterior predictive probability.
#' @param design A character string specifying the type of trial design. Options are
#'        \code{'controlled'} for randomized controlled trials, \code{'uncontrolled'}
#'        for single-arm studies, or \code{'external'} for designs incorporating
#'        external data.
#' @param theta.TV A numeric value representing the target value threshold for
#'        calculating Go probability when \code{prob = 'posterior'}. This represents
#'        the minimum clinically meaningful treatment effect.
#' @param theta.MAV A numeric value representing the minimum acceptable value threshold
#'        for calculating NoGo probability when \code{prob = 'posterior'}. This represents
#'        the futility threshold.
#' @param theta.NULL A numeric value representing the null hypothesis threshold for
#'        calculating Go/NoGo probabilities when \code{prob = 'predictive'}.
#'        Set to NULL if \code{prob = 'posterior'}.
#' @param gamma1 A numeric value in (0, 1) representing the minimum probability
#'        threshold to declare success (Go decision). Typically set to values like
#'        0.8 or 0.9.
#' @param gamma2 A numeric value in (0, 1) representing the maximum probability
#'        threshold for declaring futility (NoGo decision). Typically set to values
#'        like 0.2 or 0.3, with gamma2 < gamma1.
#' @param pi1 A numeric value or vector representing the true response probability(s)
#'        for group 1 (treatment) under which to evaluate operating characteristics.
#' @param pi2 A numeric value or vector representing the true response probability(s)
#'        for group 2 (control) under which to evaluate operating characteristics.
#' @param n1 A positive integer representing the number of patients in group 1 for
#'        the proof-of-concept (PoC) trial.
#' @param n2 A positive integer representing the number of patients in group 2 for
#'        the PoC trial (used for controlled and external designs, also required for
#'        uncontrolled design for consistency).
#' @param a1 A positive numeric value representing the first shape parameter (α)
#'        of the prior beta distribution for group 1 (treatment).
#' @param a2 A positive numeric value representing the first shape parameter (α)
#'        of the prior beta distribution for group 2 (control or hypothetical control).
#' @param b1 A positive numeric value representing the second shape parameter (β)
#'        of the prior beta distribution for group 1 (treatment).
#' @param b2 A positive numeric value representing the second shape parameter (β)
#'        of the prior beta distribution for group 2 (control or hypothetical control).
#' @param z A non-negative integer representing the hypothetical control responder count
#'        (required if \code{design = 'uncontrolled'}, otherwise set to NULL).
#'        This specifies the expected number of responders in a hypothetical control
#'        group based on historical data or prior knowledge. This is mathematically
#'        equivalent to y2 in pPPsinglebinary() for uncontrolled design.
#' @param m1 A positive integer representing the number of patients in group 1 for
#'        the future trial (required if \code{prob = 'predictive'}, otherwise set to NULL).
#' @param m2 A positive integer representing the number of patients in group 2 for
#'        the future trial (required if \code{prob = 'predictive'}, otherwise set to NULL).
#' @param ne1 A positive integer representing the number of patients in group 1 for
#'        the external data (required if \code{design = 'external'}, otherwise set to NULL).
#' @param ne2 A positive integer representing the number of patients in group 2 for
#'        the external data (required if \code{design = 'external'}, otherwise set to NULL).
#' @param ye1 A non-negative integer representing the observed number of responders
#'        in group 1 for the external data (required if \code{design = 'external'},
#'        otherwise set to NULL).
#' @param ye2 A non-negative integer representing the observed number of responders
#'        in group 2 for the external data (required if \code{design = 'external'},
#'        otherwise set to NULL).
#' @param ae1 A numeric value in (0, 1] representing the power prior scale parameter
#'        for group 1 (required if \code{design = 'external'}, otherwise set to NULL).
#'        Controls the degree of borrowing: 0 = no borrowing, 1 = full borrowing.
#' @param ae2 A numeric value in (0, 1] representing the power prior scale parameter
#'        for group 2 (required if \code{design = 'external'}, otherwise set to NULL).
#'        Controls the degree of borrowing: 0 = no borrowing, 1 = full borrowing.
#' @param error_if_Miss A logical value; if TRUE (default), the function stops with
#'        an error if Miss probability > 0 (indicating poorly chosen thresholds).
#'        If FALSE, continues execution and reports Miss probability.
#' @param Gray_inc_Miss A logical value; if TRUE, Miss probability is added to Gray
#'        probability. If FALSE (default), Miss probability is reported as a separate
#'        category. This parameter is only active when \code{error_if_Miss = FALSE}.
#'
#' @return A data frame containing the true response probabilities for both groups
#'         and the corresponding Go, NoGo, and Gray probabilities. When
#'         \code{error_if_Miss = FALSE} and \code{Gray_inc_Miss = FALSE}, Miss
#'         probability is also included as a separate column.
#'
#' @details
#' The function evaluates operating characteristics by:
#' \itemize{
#'   \item Enumerating all possible trial outcomes (y1, y2) under the specified
#'         true response rates (pi1, pi2)
#'   \item Computing the posterior or predictive probability for each outcome
#'   \item Classifying each outcome as Go, NoGo, or Gray based on decision thresholds
#'   \item Weighting each outcome by its probability under the true response rates
#' }
#'
#' **Decision rules**:
#' \itemize{
#'   \item **Go**: P(treatment effect > threshold | data) ≥ γ₁ AND
#'                 P(treatment effect > threshold | data) in lower tail < γ₂
#'   \item **NoGo**: P(treatment effect > threshold | data) < γ₁ AND
#'                   P(treatment effect > threshold | data) in lower tail ≥ γ₂
#'   \item **Gray**: Neither Go nor NoGo criteria are met (γ₂ < probability < γ₁)
#'   \item **Miss**: Both Go and NoGo criteria are met simultaneously (indicates
#'                   poorly chosen thresholds)
#' }
#'
#' **Design-specific handling**:
#' \itemize{
#'   \item **Controlled design**: Evaluates all possible outcomes (y1, y2) where
#'         y1 ∈ {0,...,n1} and y2 ∈ {0,...,n2}
#'   \item **Uncontrolled design**: Uses hypothetical control specified by \code{z}.
#'         Evaluates outcomes (y1, z) where y1 ∈ {0,...,n1} and y2 is fixed at z.
#'         The parameter z represents the expected responder count in a hypothetical
#'         control group.
#'   \item **External design**: Incorporates external data through power priors for
#'         both treatment and control groups
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
#'   \item **Uncontrolled design**: Single-arm trial with historical control
#'   \item **External design**: Incorporating historical data through power priors
#' }
#'
#' @examples
#' # Example 1: Controlled design with posterior probability
#' pGNGsinglebinary(
#'   prob = 'posterior', design = 'controlled',
#'   theta.TV = 0.4, theta.MAV = 0.2, theta.NULL = NULL,
#'   gamma1 = 0.8, gamma2 = 0.2,
#'   pi1 = c(0.2, 0.4, 0.6, 0.8), pi2 = rep(0.2, 4),
#'   n1 = 12, n2 = 12,
#'   a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
#'   z = NULL, m1 = NULL, m2 = NULL,
#'   ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE
#' )
#'
#' # Example 2: Uncontrolled design with hypothetical control
#' # z = 5 means we expect 5 responders in hypothetical control (n2 = 15)
#' pGNGsinglebinary(
#'   prob = 'posterior', design = 'uncontrolled',
#'   theta.TV = 0.30, theta.MAV = 0.15, theta.NULL = NULL,
#'   gamma1 = 0.75, gamma2 = 0.25,
#'   pi1 = c(0.3, 0.5, 0.7), pi2 = rep(0.33, 3),
#'   n1 = 15, n2 = 15,
#'   a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
#'   z = 5, m1 = NULL, m2 = NULL,
#'   ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE
#' )
#'
#' # Example 3: Posterior predictive probability for controlled design
#' pGNGsinglebinary(
#'   prob = 'predictive', design = 'controlled',
#'   theta.TV = NULL, theta.MAV = NULL, theta.NULL = 0,
#'   gamma1 = 0.9, gamma2 = 0.3,
#'   pi1 = c(0.2, 0.4, 0.6, 0.8), pi2 = rep(0.2, 4),
#'   n1 = 12, n2 = 12,
#'   a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
#'   z = NULL, m1 = 30, m2 = 30,
#'   ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE
#' )
#'
#' @importFrom stats dbinom
#' @export
pGNGsinglebinary <- function(prob = 'posterior', design = 'controlled',
                             theta.TV, theta.MAV, theta.NULL = NULL,
                             gamma1, gamma2, pi1, pi2, n1, n2, a1, a2, b1, b2,
                             z = NULL, m1, m2, ne1, ne2, ye1, ye2, ae1, ae2,
                             error_if_Miss = TRUE, Gray_inc_Miss = FALSE) {
  # Validate parameter sets for posterior probability
  if((prob == 'posterior') & (sum(sapply(list(theta.TV, theta.MAV), is.null)) > 0)) {
    stop('If you calculate Go, NoGo, and Gray probabilities using posterior probability, theta.TV and theta.MAV must be non-null')
  }

  # Validate parameter sets for posterior predictive probability
  if((prob == 'predictive') & (sum(sapply(list(m1, m2), is.null)) > 0)) {
    stop('If you calculate Go, NoGo, and Gray probabilities using posterior predictive probability, m1 and m2 must be non-null')
  }

  # Validate parameter sets for posterior predictive probability threshold
  if((prob == 'predictive') & (is.null(theta.NULL))) {
    stop('If you calculate Go, NoGo, and Gray probabilities using posterior predictive probability, theta.NULL must be non-null')
  }

  # Validate parameter sets for uncontrolled design
  if((design == 'uncontrolled') & (is.null(z))) {
    stop('If you consider uncontrolled design, z (hypothetical control responder count) must be non-null')
  }

  # Validate parameter sets for external design
  if((design == 'external') & (sum(sapply(list(ne1, ne2, ye1, ye2, ae1, ae2), is.null)) > 0)) {
    stop('If you use external data, ne1, ne2, ye1, ye2, ae1, and ae2 must be non-null')
  }

  # Set threshold values based on probability type
  if(prob == 'posterior') {
    # For posterior probability: use theta.TV for Go, theta.MAV for NoGo
    theta0 <- c(theta.TV, theta.MAV)
  } else {
    # For predictive probability: use theta.NULL for both
    theta0 <- c(theta.NULL, theta.NULL)
  }

  # Define possible outcomes for group 1 (treatment)
  Y1 <- 0:n1

  # Define possible outcomes for group 2 (control or fixed hypothetical value)
  if(design == 'uncontrolled') {
    # For uncontrolled: y2 is fixed at hypothetical value z
    Y2 <- z
  } else {
    # For controlled/external: y2 ranges from 0 to n2
    Y2 <- 0:n2
  }

  # Calculate posterior/predictive probabilities for each threshold
  # gPost[[1]]: probability for Go threshold (theta.TV or theta.NULL)
  # gPost[[2]]: probability for NoGo threshold (theta.MAV or theta.NULL)
  gPost <- lapply(seq(length(theta0)), function(i) {
    sapply(Y1, function(y1) {
      sapply(Y2, function(y2) {
        pPPsinglebinary(
          prob, design, theta0[i],
          n1, n2, y1, y2, a1, a2, b1, b2,
          m1, m2, ne1, ne2, ye1, ye2, ae1, ae2, c(FALSE, TRUE)[i]
        )
      })
    })
  })

  # Calculate Go, NoGo, and Miss probabilities based on decision criteria
  GoNogoProb <- matrix(
    sapply(seq(3), function(j) {

      # Create indicator matrix for each decision type
      if(j == 1) {
        # Go decision: high probability for Go threshold AND low probability for NoGo threshold
        I <- matrix((gPost[[1]] >= gamma1) & (gPost[[2]] < gamma2), nrow = length(Y2))
      } else if(j == 2) {
        # NoGo decision: low probability for Go threshold AND high probability for NoGo threshold
        I <- matrix((gPost[[1]] < gamma1) & (gPost[[2]] >= gamma2), nrow = length(Y2))
      } else {
        # Miss: both Go and NoGo criteria met simultaneously (should be rare)
        I <- matrix((gPost[[1]] >= gamma1) & (gPost[[2]] >= gamma2), nrow = length(Y2))
      }

      if(design == 'uncontrolled') {
        # For uncontrolled design: sum over group 1 outcomes only (y2 is fixed at z)
        colSums(outer(col(I)[I] - 1, pi1, FUN = function(X, Y) dbinom(X, n1, Y)))
      } else {
        # For controlled/external design: sum over both group outcomes
        diag(crossprod(
          outer(col(I)[I] - 1, pi1, FUN = function(X, Y) dbinom(X, n1, Y)),
          outer(row(I)[I] - 1, pi2, FUN = function(X, Y) dbinom(X, n2, Y))
        ))
      }
    }),
    ncol = 3
  )

  # Check for positive Miss probabilities (indicates inappropriate thresholds)
  if(error_if_Miss) {
    if(sum(GoNogoProb[, 3, drop = FALSE]) > 0) {
      stop('Because positive Miss probability(s) is obtained, re-consider appropriate thresholds')
    }
  }

  # Calculate Miss probability (both Go and NoGo criteria met simultaneously)
  Miss <- GoNogoProb[, 3, drop = FALSE]

  # Calculate Gray probability (complement of Go and NoGo)
  if(Gray_inc_Miss) {
    # Include Miss in Gray probability
    GrayProb <- 1 - rowSums(GoNogoProb[, -3, drop = FALSE])
  } else {
    # Exclude Miss from Gray probability
    GrayProb <- 1 - rowSums(GoNogoProb[, drop = FALSE])
  }

  # Prepare results data frame
  results <- data.frame(
    pi1, pi2,
    Go = GoNogoProb[, 1],
    Gray = GrayProb,
    NoGo = GoNogoProb[, 2]
  )

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
