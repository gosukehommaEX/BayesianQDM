#' Calculate Go, NoGo, and Gray Probabilities for a Clinical Trial with Two Binary Endpoints
#' Under the Bayesian Framework
#'
#' This function calculates the Go, NoGo, and Gray probabilities for clinical trials
#' with two binary endpoints using the Bayesian framework. The function evaluates
#' operating characteristics by computing the probability of making each type of
#' decision (Go, NoGo, or Gray) across different true response rates. The function
#' supports controlled, uncontrolled, and external control designs. Decision criteria
#' can be based on joint probabilities across the nine regions defined by target values
#' (TV) and minimum acceptable values (MAV) for both endpoints.
#'
#' @param nsim A positive integer representing the number of simulation iterations
#'        to estimate Go/NoGo/Gray probabilities. Default is 10000.
#' @param prob A character string specifying the type of probability to use for
#'        decision-making. Options are \code{'posterior'} for posterior probability
#'        or \code{'predictive'} for posterior predictive probability.
#' @param design A character string specifying the type of trial design. Options are
#'        \code{'controlled'} for randomized controlled trials, \code{'uncontrolled'}
#'        for single-arm studies, or \code{'external'} for designs incorporating
#'        external data.
#' @param theta.TV1 A numeric value representing the target value threshold for
#'        Endpoint 1 when \code{prob = 'posterior'}. This defines the upper boundary
#'        for the "Gray" region.
#' @param theta.MAV1 A numeric value representing the minimum acceptable value threshold
#'        for Endpoint 1 when \code{prob = 'posterior'}. This defines the lower boundary
#'        for the "Gray" region.
#' @param theta.TV2 A numeric value representing the target value threshold for
#'        Endpoint 2 when \code{prob = 'posterior'}. This defines the upper boundary
#'        for the "Gray" region.
#' @param theta.MAV2 A numeric value representing the minimum acceptable value threshold
#'        for Endpoint 2 when \code{prob = 'posterior'}. This defines the lower boundary
#'        for the "Gray" region.
#' @param theta.NULL1 A numeric value representing the null hypothesis threshold for
#'        Endpoint 1 when \code{prob = 'predictive'}. Required when \code{prob = 'predictive'},
#'        set to NULL if \code{prob = 'posterior'}.
#' @param theta.NULL2 A numeric value representing the null hypothesis threshold for
#'        Endpoint 2 when \code{prob = 'predictive'}. Required when \code{prob = 'predictive'},
#'        set to NULL if \code{prob = 'posterior'}.
#' @param gammaG1 A numeric value in (0, 1) representing the Go threshold for Endpoint 1.
#'        Must satisfy gammaG1 > gammaN1. Typically set to 0.7-0.9.
#' @param gammaN1 A numeric value in (0, 1) representing the NoGo threshold for Endpoint 1.
#'        Must satisfy gammaN1 < gammaG1. Typically set to 0.2-0.4.
#' @param gammaG2 A numeric value in (0, 1) representing the Go threshold for Endpoint 2.
#'        Must satisfy gammaG2 > gammaN2. Typically set to 0.7-0.9.
#' @param gammaN2 A numeric value in (0, 1) representing the NoGo threshold for Endpoint 2.
#'        Must satisfy gammaN2 < gammaG2. Typically set to 0.2-0.4.
#' @param decision_criteria An integer (1, 2, 3, or 4) specifying the decision rule:
#'        \itemize{
#'          \item \strong{1}: Equations (4)-(5) - Go if BOTH endpoints meet criteria (AND),
#'                NoGo if BOTH endpoints meet criteria (AND)
#'          \item \strong{2}: Equations (6)-(7) - Go if BOTH endpoints meet criteria (AND),
#'                NoGo if EITHER endpoint meets criteria (OR)
#'          \item \strong{3}: Equations (8)-(9) - Go if EITHER endpoint meets criteria (OR),
#'                NoGo if BOTH endpoints meet criteria (AND)
#'          \item \strong{4}: Equations (10)-(11) - Go if EITHER endpoint meets criteria (OR),
#'                NoGo if EITHER endpoint meets criteria (OR)
#'        }
#' @param pi1 A numeric vector of length 4 representing the true joint response
#'        probabilities for group 1 (treatment). Elements correspond to (0,0), (0,1),
#'        (1,0), (1,1) for (Endpoint1, Endpoint2). Must sum to 1.
#' @param pi2 A numeric vector of length 4 representing the true joint response
#'        probabilities for group 2 (control). Elements correspond to (0,0), (0,1),
#'        (1,0), (1,1) for (Endpoint1, Endpoint2). Must sum to 1.
#' @param a1_00 A positive numeric value representing the Dirichlet prior parameter
#'        for response pattern (0,0) in group 1.
#' @param a1_01 A positive numeric value representing the Dirichlet prior parameter
#'        for response pattern (0,1) in group 1.
#' @param a1_10 A positive numeric value representing the Dirichlet prior parameter
#'        for response pattern (1,0) in group 1.
#' @param a1_11 A positive numeric value representing the Dirichlet prior parameter
#'        for response pattern (1,1) in group 1.
#' @param a2_00 A positive numeric value representing the Dirichlet prior parameter
#'        for response pattern (0,0) in group 2.
#' @param a2_01 A positive numeric value representing the Dirichlet prior parameter
#'        for response pattern (0,1) in group 2.
#' @param a2_10 A positive numeric value representing the Dirichlet prior parameter
#'        for response pattern (1,0) in group 2.
#' @param a2_11 A positive numeric value representing the Dirichlet prior parameter
#'        for response pattern (1,1) in group 2.
#' @param n1 A positive integer representing the number of patients in group 1 for
#'        the proof-of-concept (PoC) trial.
#' @param n2 A positive integer representing the number of patients in group 2 for
#'        the PoC trial.
#' @param m1 A positive integer representing the future sample size for group 1 in
#'        the Phase III trial. Required when \code{prob = 'predictive'}, NULL otherwise.
#' @param m2 A positive integer representing the future sample size for group 2 in
#'        the Phase III trial. Required when \code{prob = 'predictive'}, NULL otherwise.
#' @param z A numeric vector of length 4 representing the hypothetical control response
#'        counts for \code{design = 'uncontrolled'}. Elements correspond to (0,0), (0,1),
#'        (1,0), (1,1) for (Endpoint1, Endpoint2). Must sum to n2. This specifies the
#'        expected number of responses in each category for a hypothetical control group
#'        based on historical data or prior knowledge. Set to NULL for controlled and
#'        external designs.
#' @param xe1_00 A non-negative integer representing external control data count for
#'        response pattern (0,0) in group 1. Used when \code{design = 'external'}.
#' @param xe1_01 A non-negative integer representing external control data count for
#'        response pattern (0,1) in group 1. Used when \code{design = 'external'}.
#' @param xe1_10 A non-negative integer representing external control data count for
#'        response pattern (1,0) in group 1. Used when \code{design = 'external'}.
#' @param xe1_11 A non-negative integer representing external control data count for
#'        response pattern (1,1) in group 1. Used when \code{design = 'external'}.
#' @param xe2_00 A non-negative integer representing external control data count for
#'        response pattern (0,0) in group 2. Used when \code{design = 'external'}.
#' @param xe2_01 A non-negative integer representing external control data count for
#'        response pattern (0,1) in group 2. Used when \code{design = 'external'}.
#' @param xe2_10 A non-negative integer representing external control data count for
#'        response pattern (1,0) in group 2. Used when \code{design = 'external'}.
#' @param xe2_11 A non-negative integer representing external control data count for
#'        response pattern (1,1) in group 2. Used when \code{design = 'external'}.
#' @param ae1 A numeric value in \[0, 1\] representing the power prior parameter for
#'        external data in group 1. Used when \code{design = 'external'}.
#' @param ae2 A numeric value in \[0, 1\] representing the power prior parameter for
#'        external data in group 2. Used when \code{design = 'external'}.
#' @param Go1_regions A numeric vector specifying which regions count as Go regions
#'        for Endpoint 1. For \code{prob = 'posterior'}: values should be between 1-9.
#'        For \code{prob = 'predictive'}: values should be between 1-4.
#'        Default is \code{c(1,2,3)} for posterior or \code{c(1,2)} for predictive.
#' @param NoGo1_regions A numeric vector specifying which regions count as NoGo regions
#'        for Endpoint 1. For \code{prob = 'posterior'}: values should be between 1-9.
#'        For \code{prob = 'predictive'}: values should be between 1-4.
#'        Default is \code{c(7,8,9)} for posterior or \code{c(3,4)} for predictive.
#' @param Go2_regions A numeric vector specifying which regions count as Go regions
#'        for Endpoint 2. For \code{prob = 'posterior'}: values should be between 1-9.
#'        For \code{prob = 'predictive'}: values should be between 1-4.
#'        Default is \code{c(1,4,7)} for posterior or \code{c(1,3)} for predictive.
#' @param NoGo2_regions A numeric vector specifying which regions count as NoGo regions
#'        for Endpoint 2. For \code{prob = 'posterior'}: values should be between 1-9.
#'        For \code{prob = 'predictive'}: values should be between 1-4.
#'        Default is \code{c(3,6,9)} for posterior or \code{c(2,4)} for predictive.
#' @param error_if_Miss A logical value indicating whether to stop with an error if
#'        positive Miss probability is detected. Default is TRUE.
#' @param Gray_inc_Miss A logical value indicating whether to include Miss probability
#'        in Gray probability calculation. Default is TRUE.
#' @param nMC A positive integer representing the number of Monte Carlo samples for
#'        probability estimation within each simulation iteration. Default is 10000.
#' @param seed An integer value for random number generator seed to ensure reproducibility.
#'        If NULL (default), no seed is set.
#'
#' @return A data frame containing the following columns:
#' \describe{
#'   \item{pi1_00, pi1_01, pi1_10, pi1_11}{The true joint probabilities for group 1}
#'   \item{pi2_00, pi2_01, pi2_10, pi2_11}{The true joint probabilities for group 2}
#'   \item{Go}{The probability of making a Go decision}
#'   \item{Gray}{The probability of making a Gray (inconclusive) decision}
#'   \item{NoGo}{The probability of making a NoGo decision}
#'   \item{Miss}{(Optional) The probability where both Go and NoGo criteria are met
#'               simultaneously. Included when \code{error_if_Miss = FALSE} and
#'               \code{Gray_inc_Miss = FALSE}}
#' }
#'
#' @details
#' The function simulates clinical trial data from multinomial distributions with
#' specified true response probabilities, then calculates posterior or posterior
#' predictive probabilities for each simulated dataset using \code{pPPtwobinary()}.
#'
#' **Decision Regions for Posterior Probability (9 regions)**:
#' The nine regions are defined by the thresholds for both endpoints as follows
#' (regions numbered 1-9 in column-major order):
#' \tabular{lccc}{
#'   \tab \strong{θ1 > TV1} \tab \strong{TV1 ≥ θ1 > MAV1} \tab \strong{MAV1 ≥ θ1} \cr
#'   \strong{θ2 > TV2} \tab Region 1 \tab Region 4 \tab Region 7 \cr
#'   \strong{TV2 ≥ θ2 > MAV2} \tab Region 2 \tab Region 5 \tab Region 8 \cr
#'   \strong{MAV2 ≥ θ2} \tab Region 3 \tab Region 6 \tab Region 9
#' }
#'
#' **Decision Regions for Predictive Probability (4 regions)**:
#' The four regions are defined by the NULL hypothesis thresholds for both endpoints
#' (regions numbered 1-4 in column-major order):
#' \tabular{lcc}{
#'   \tab \strong{θ1 > NULL1} \tab \strong{NULL1 ≥ θ1} \cr
#'   \strong{θ2 > NULL2} \tab Region 1 \tab Region 3 \cr
#'   \strong{NULL2 ≥ θ2} \tab Region 2 \tab Region 4
#' }
#'
#' **Decision Criteria** (corresponding to Documentation.pdf equations):
#' \itemize{
#'   \item **Criteria 1** - Equations (4)-(5): AND/AND decision
#'     \itemize{
#'       \item Go: Pr(Endpoint1 success) >= gammaG1 AND Pr(Endpoint2 success) >= gammaG2
#'       \item NoGo: Pr(Endpoint1 failure) >= gammaN1 AND Pr(Endpoint2 failure) >= gammaN2
#'     }
#'   \item **Criteria 2** - Equations (6)-(7): AND/OR decision
#'     \itemize{
#'       \item Go: Pr(Endpoint1 success) >= gammaG1 AND Pr(Endpoint2 success) >= gammaG2
#'       \item NoGo: Pr(Endpoint1 failure) >= gammaN1 OR Pr(Endpoint2 failure) >= gammaN2
#'     }
#'   \item **Criteria 3** - Equations (8)-(9): OR/AND decision
#'     \itemize{
#'       \item Go: Pr(Endpoint1 success) >= gammaG1 OR Pr(Endpoint2 success) >= gammaG2
#'       \item NoGo: Pr(Endpoint1 failure) >= gammaN1 AND Pr(Endpoint2 failure) >= gammaN2
#'     }
#'   \item **Criteria 4** - Equations (10)-(11): OR/OR decision
#'     \itemize{
#'       \item Go: Pr(Endpoint1 success) >= gammaG1 OR Pr(Endpoint2 success) >= gammaG2
#'       \item NoGo: Pr(Endpoint1 failure) >= gammaN1 OR Pr(Endpoint2 failure) >= gammaN2
#'     }
#' }
#'
#' **Gray Decision**: Made when neither Go nor NoGo criteria are met (the complement).
#'
#' **Miss Decision**: Occurs when both Go and NoGo criteria are met simultaneously,
#' which should be rare with appropriate threshold selection.
#'
#' @examples
#' # Example 1: Decision criteria 1 - Equations (4)-(5): AND/AND decision
#' result1 <- pGNGtwobinary(
#'   nsim = 1000, prob = 'posterior', design = 'controlled',
#'   theta.TV1 = 0.20, theta.MAV1 = 0.10,
#'   theta.TV2 = 0.20, theta.MAV2 = 0.10,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   gammaG1 = 0.80, gammaN1 = 0.20,
#'   gammaG2 = 0.80, gammaN2 = 0.20,
#'   decision_criteria = 1,
#'   pi1 = c(0.25, 0.20, 0.20, 0.35),
#'   pi2 = c(0.40, 0.25, 0.25, 0.10),
#'   a1_00 = 0.5, a1_01 = 0.5, a1_10 = 0.5, a1_11 = 0.5,
#'   a2_00 = 0.5, a2_01 = 0.5, a2_10 = 0.5, a2_11 = 0.5,
#'   n1 = 15, n2 = 15, m1 = NULL, m2 = NULL,
#'   xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
#'   xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
#'   ae1 = NULL, ae2 = NULL,
#'   nMC = 10000, seed = 123
#' )
#' print(result1)
#'
#' # Example 2: Decision criteria 2 - Equations (6)-(7): AND/OR decision
#' result2 <- pGNGtwobinary(
#'   nsim = 1000, prob = 'posterior', design = 'controlled',
#'   theta.TV1 = 0.20, theta.MAV1 = 0.10,
#'   theta.TV2 = 0.20, theta.MAV2 = 0.10,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   gammaG1 = 0.70, gammaN1 = 0.30,
#'   gammaG2 = 0.70, gammaN2 = 0.30,
#'   decision_criteria = 2,
#'   pi1 = c(0.25, 0.20, 0.20, 0.35),
#'   pi2 = c(0.40, 0.25, 0.25, 0.10),
#'   a1_00 = 0.5, a1_01 = 0.5, a1_10 = 0.5, a1_11 = 0.5,
#'   a2_00 = 0.5, a2_01 = 0.5, a2_10 = 0.5, a2_11 = 0.5,
#'   n1 = 15, n2 = 15, m1 = NULL, m2 = NULL,
#'   xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
#'   xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
#'   ae1 = NULL, ae2 = NULL,
#'   nMC = 10000, seed = 124
#' )
#' print(result2)
#'
#' # Example 3: Decision criteria 3 - Equations (8)-(9): OR/AND decision
#' result3 <- pGNGtwobinary(
#'   nsim = 1000, prob = 'posterior', design = 'controlled',
#'   theta.TV1 = 0.20, theta.MAV1 = 0.10,
#'   theta.TV2 = 0.20, theta.MAV2 = 0.10,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   gammaG1 = 0.75, gammaN1 = 0.25,
#'   gammaG2 = 0.75, gammaN2 = 0.25,
#'   decision_criteria = 3,
#'   pi1 = c(0.25, 0.20, 0.20, 0.35),
#'   pi2 = c(0.40, 0.25, 0.25, 0.10),
#'   a1_00 = 0.5, a1_01 = 0.5, a1_10 = 0.5, a1_11 = 0.5,
#'   a2_00 = 0.5, a2_01 = 0.5, a2_10 = 0.5, a2_11 = 0.5,
#'   n1 = 15, n2 = 15, m1 = NULL, m2 = NULL,
#'   xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
#'   xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
#'   ae1 = NULL, ae2 = NULL,
#'   nMC = 10000, seed = 125
#' )
#' print(result3)
#'
#' # Example 4: Decision criteria 4 - Equations (10)-(11): OR/OR decision
#' result4 <- pGNGtwobinary(
#'   nsim = 500, prob = 'posterior', design = 'controlled',
#'   theta.TV1 = 0.15, theta.MAV1 = 0.08,
#'   theta.TV2 = 0.15, theta.MAV2 = 0.08,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   gammaG1 = 0.80, gammaN1 = 0.20,
#'   gammaG2 = 0.80, gammaN2 = 0.20,
#'   decision_criteria = 4,
#'   pi1 = c(0.30, 0.18, 0.18, 0.34),
#'   pi2 = c(0.45, 0.22, 0.22, 0.11),
#'   a1_00 = 0.5, a1_01 = 0.5, a1_10 = 0.5, a1_11 = 0.5,
#'   a2_00 = 0.5, a2_01 = 0.5, a2_10 = 0.5, a2_11 = 0.5,
#'   n1 = 15, n2 = 15, m1 = NULL, m2 = NULL,
#'   xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
#'   xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
#'   ae1 = NULL, ae2 = NULL,
#'   error_if_Miss = FALSE, Gray_inc_Miss = FALSE,
#'   nMC = 10000, seed = 126
#' )
#' print(result4)
#'
#' # Example 5: Posterior predictive probability with decision criteria 1
#' result5 <- pGNGtwobinary(
#'   nsim = 1000, prob = 'predictive', design = 'controlled',
#'   theta.TV1 = 0.15, theta.MAV1 = 0.15,
#'   theta.TV2 = 0.15, theta.MAV2 = 0.15,
#'   theta.NULL1 = 0.15, theta.NULL2 = 0.15,
#'   gammaG1 = 0.75, gammaN1 = 0.25,
#'   gammaG2 = 0.75, gammaN2 = 0.25,
#'   decision_criteria = 1,
#'   pi1 = c(0.30, 0.15, 0.15, 0.40),
#'   pi2 = c(0.45, 0.20, 0.20, 0.15),
#'   a1_00 = 1, a1_01 = 1, a1_10 = 1, a1_11 = 1,
#'   a2_00 = 1, a2_01 = 1, a2_10 = 1, a2_11 = 1,
#'   n1 = 20, n2 = 20, m1 = 50, m2 = 50,
#'   xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
#'   xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
#'   ae1 = NULL, ae2 = NULL,
#'   nMC = 10000, seed = 127
#' )
#' print(result5)
#'
#' # Example 6: External control design with decision criteria 2
#' result6 <- pGNGtwobinary(
#'   nsim = 500, prob = 'posterior', design = 'external',
#'   theta.TV1 = 0.15, theta.MAV1 = 0.08,
#'   theta.TV2 = 0.15, theta.MAV2 = 0.08,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   gammaG1 = 0.80, gammaN1 = 0.20,
#'   gammaG2 = 0.80, gammaN2 = 0.20,
#'   decision_criteria = 2,
#'   pi1 = c(0.30, 0.18, 0.18, 0.34),
#'   pi2 = c(0.45, 0.22, 0.22, 0.11),
#'   a1_00 = 0.5, a1_01 = 0.5, a1_10 = 0.5, a1_11 = 0.5,
#'   a2_00 = 0.5, a2_01 = 0.5, a2_10 = 0.5, a2_11 = 0.5,
#'   n1 = 15, n2 = 15, m1 = NULL, m2 = NULL,
#'   xe1_00 = 5, xe1_01 = 3, xe1_10 = 3, xe1_11 = 9,
#'   xe2_00 = 10, xe2_01 = 5, xe2_10 = 5, xe2_11 = 5,
#'   ae1 = 0.5, ae2 = 0.5,
#'   nMC = 10000, seed = 128
#' )
#' print(result6)
#'
#' # Example 7: Uncontrolled design with hypothetical control and decision criteria 3
#' result7 <- pGNGtwobinary(
#'   nsim = 1000, prob = 'posterior', design = 'uncontrolled',
#'   theta.TV1 = 0.20, theta.MAV1 = 0.10,
#'   theta.TV2 = 0.20, theta.MAV2 = 0.10,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   gammaG1 = 0.80, gammaN1 = 0.20,
#'   gammaG2 = 0.80, gammaN2 = 0.20,
#'   decision_criteria = 3,
#'   pi1 = c(0.25, 0.20, 0.20, 0.35),
#'   pi2 = c(0.40, 0.25, 0.25, 0.10),  # Not used for uncontrolled
#'   a1_00 = 0.5, a1_01 = 0.5, a1_10 = 0.5, a1_11 = 0.5,
#'   a2_00 = 0.5, a2_01 = 0.5, a2_10 = 0.5, a2_11 = 0.5,
#'   n1 = 20, n2 = 20, m1 = NULL, m2 = NULL,
#'   z = c(8, 5, 5, 2),  # Hypothetical control: 8 for (0,0), 5 for (0,1), etc.
#'   xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
#'   xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
#'   ae1 = NULL, ae2 = NULL,
#'   nMC = 10000, seed = 129
#' )
#' print(result7)
#'
#' @importFrom stats rmultinom
#' @export
pGNGtwobinary <- function(nsim = 10000, prob = 'posterior', design = 'controlled',
                          theta.TV1, theta.MAV1, theta.TV2, theta.MAV2,
                          theta.NULL1 = NULL, theta.NULL2 = NULL,
                          gammaG1 = NULL, gammaN1 = NULL,
                          gammaG2 = NULL, gammaN2 = NULL,
                          decision_criteria = 1,
                          pi1, pi2,
                          a1_00, a1_01, a1_10, a1_11,
                          a2_00, a2_01, a2_10, a2_11,
                          n1, n2, m1 = NULL, m2 = NULL, z = NULL,
                          xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
                          xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
                          ae1 = NULL, ae2 = NULL,
                          Go1_regions = NULL, NoGo1_regions = NULL,
                          Go2_regions = NULL, NoGo2_regions = NULL,
                          error_if_Miss = TRUE, Gray_inc_Miss = TRUE,
                          nMC = 10000, seed = NULL) {

  # Set seed for reproducibility if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Set default regions if not specified by user
  if (is.null(Go1_regions) || is.null(NoGo1_regions) || is.null(Go2_regions) || is.null(NoGo2_regions)) {
    if (prob == 'posterior') {
      # Default regions for posterior probability (9 regions)
      Go1_regions <- c(1, 2, 3)
      NoGo1_regions <- c(7, 8, 9)
      Go2_regions <- c(1, 4, 7)
      NoGo2_regions <- c(3, 6, 9)
    } else {
      # Default regions for predictive probability (4 regions)
      Go1_regions <- c(1, 2)
      NoGo1_regions <- c(3, 4)
      Go2_regions <- c(1, 3)
      NoGo2_regions <- c(2, 4)
    }
  }

  # Validate decision_criteria parameter
  if (!decision_criteria %in% c(1, 2, 3, 4)) {
    stop("decision_criteria must be 1, 2, 3, or 4")
  }

  # Validate gamma parameters - all decision criteria use endpoint-specific gammas
  if (is.null(gammaG1) || is.null(gammaN1) || is.null(gammaG2) || is.null(gammaN2)) {
    stop("gammaG1, gammaN1, gammaG2, and gammaN2 must be specified")
  }
  if (gammaG1 <= gammaN1) {
    stop("gammaG1 must be greater than gammaN1")
  }
  if (gammaG2 <= gammaN2) {
    stop("gammaG2 must be greater than gammaN2")
  }

  # Input validation
  if (!is.numeric(nsim) || nsim <= 0 || nsim != as.integer(nsim)) {
    stop("nsim must be a positive integer")
  }

  if (!prob %in% c('posterior', 'predictive')) {
    stop("prob must be either 'posterior' or 'predictive'")
  }

  if (!design %in% c('controlled', 'external', 'uncontrolled')) {
    stop("design must be 'controlled', 'external', or 'uncontrolled'")
  }

  if (length(pi1) != 4 || abs(sum(pi1) - 1) > 1e-10) {
    stop("pi1 must be a numeric vector of length 4 that sums to 1")
  }

  if (length(pi2) != 4 || abs(sum(pi2) - 1) > 1e-10) {
    stop("pi2 must be a numeric vector of length 4 that sums to 1")
  }

  # Validate region parameters based on prob type
  if (prob == 'posterior') {
    if (any(Go1_regions < 1 | Go1_regions > 9)) {
      stop("For posterior probability, Go1_regions must contain values between 1 and 9")
    }
    if (any(NoGo1_regions < 1 | NoGo1_regions > 9)) {
      stop("For posterior probability, NoGo1_regions must contain values between 1 and 9")
    }
    if (any(Go2_regions < 1 | Go2_regions > 9)) {
      stop("For posterior probability, Go2_regions must contain values between 1 and 9")
    }
    if (any(NoGo2_regions < 1 | NoGo2_regions > 9)) {
      stop("For posterior probability, NoGo2_regions must contain values between 1 and 9")
    }
  } else {
    if (any(Go1_regions < 1 | Go1_regions > 4)) {
      stop("For predictive probability, Go1_regions must contain values between 1 and 4")
    }
    if (any(NoGo1_regions < 1 | NoGo1_regions > 4)) {
      stop("For predictive probability, NoGo1_regions must contain values between 1 and 4")
    }
    if (any(Go2_regions < 1 | Go2_regions > 4)) {
      stop("For predictive probability, Go2_regions must contain values between 1 and 4")
    }
    if (any(NoGo2_regions < 1 | NoGo2_regions > 4)) {
      stop("For predictive probability, NoGo2_regions must contain values between 1 and 4")
    }
  }

  # Validate theta.NULL parameters for predictive probability
  if (prob == 'predictive') {
    if (is.null(theta.NULL1) || is.null(theta.NULL2)) {
      stop("theta.NULL1 and theta.NULL2 must be specified when prob = 'predictive'")
    }
  }

  # Validate parameter sets for uncontrolled design
  if (design == 'uncontrolled') {
    if (is.null(z)) {
      stop("For uncontrolled design, z (hypothetical control response counts) must be specified")
    }
    if (length(z) != 4) {
      stop("z must be a numeric vector of length 4 for (0,0), (0,1), (1,0), (1,1) responses")
    }
    if (sum(z) != n2) {
      stop("Sum of z must equal n2")
    }
  }

  # Set threshold values based on probability type
  if (prob == 'posterior') {
    theta01 <- c(theta.TV1, theta.MAV1)
    theta02 <- c(theta.TV2, theta.MAV2)
  } else {
    # For predictive probability, TV and MAV should be equal to NULL
    theta01 <- c(theta.NULL1, theta.NULL1)
    theta02 <- c(theta.NULL2, theta.NULL2)
  }

  # Generate simulated data
  x1 <- t(rmultinom(nsim, n1, pi1))
  if (design != 'uncontrolled') {
    x2 <- t(rmultinom(nsim, n2, pi2))
  } else {
    # For uncontrolled design, use fixed hypothetical control counts
    x2 <- matrix(rep(z, nsim), nrow = nsim, byrow = TRUE)
  }

  # Calculate posterior/predictive probabilities for each simulation
  pPost_and_Pred <- do.call(
    rbind,
    lapply(seq(nsim), function(i) {
      pPPtwobinary(
        prob = prob, design = design,
        theta.TV1 = theta01[1], theta.MAV1 = theta01[2],
        theta.TV2 = theta02[1], theta.MAV2 = theta02[2],
        x1_00 = x1[i, 1], x1_01 = x1[i, 2], x1_10 = x1[i, 3], x1_11 = x1[i, 4],
        x2_00 = x2[i, 1], x2_01 = x2[i, 2], x2_10 = x2[i, 3], x2_11 = x2[i, 4],
        a1_00 = a1_00, a1_01 = a1_01, a1_10 = a1_10, a1_11 = a1_11,
        a2_00 = a2_00, a2_01 = a2_01, a2_10 = a2_10, a2_11 = a2_11,
        m1 = m1, m2 = m2,
        xe1_00 = xe1_00, xe1_01 = xe1_01, xe1_10 = xe1_10, xe1_11 = xe1_11,
        xe2_00 = xe2_00, xe2_01 = xe2_01, xe2_10 = xe2_10, xe2_11 = xe2_11,
        ae1 = ae1, ae2 = ae2,
        nMC = nMC
      )
    })
  )

  # Calculate Go, NoGo and Miss probabilities based on decision criteria
  # Use the regions defined at the beginning of the function (either user-specified or defaults)

  Go_prob1 <- rowSums(pPost_and_Pred[, Go1_regions, drop = FALSE])
  NoGo_prob1 <- rowSums(pPost_and_Pred[, NoGo1_regions, drop = FALSE])
  Go_prob2 <- rowSums(pPost_and_Pred[, Go2_regions, drop = FALSE])
  NoGo_prob2 <- rowSums(pPost_and_Pred[, NoGo2_regions, drop = FALSE])

  if (decision_criteria == 1) {
    # Criteria 1 (4)-(5): Go if BOTH endpoints meet criteria (AND)
    #                     NoGo if BOTH endpoints meet criteria (AND)
    GoNogoProb <- sapply(seq(3), function(j) {
      if (j == 1) {
        # Go decision: (Go_prob1 >= gammaG1) AND (Go_prob2 >= gammaG2)
        I <- as.numeric((Go_prob1 >= gammaG1) & (Go_prob2 >= gammaG2) &
                          !((NoGo_prob1 >= gammaN1) & (NoGo_prob2 >= gammaN2)))
      } else if (j == 2) {
        # NoGo decision: (NoGo_prob1 >= gammaN1) AND (NoGo_prob2 >= gammaN2)
        I <- as.numeric((NoGo_prob1 >= gammaN1) & (NoGo_prob2 >= gammaN2) &
                          !((Go_prob1 >= gammaG1) & (Go_prob2 >= gammaG2)))
      } else {
        # Miss: Go and NoGo both satisfied
        I <- as.numeric((Go_prob1 >= gammaG1) & (Go_prob2 >= gammaG2) &
                          (NoGo_prob1 >= gammaN1) & (NoGo_prob2 >= gammaN2))
      }
      sum(I) / nsim
    })

  } else if (decision_criteria == 2) {
    # Criteria 2 (6)-(7): Go if BOTH endpoints meet criteria (AND)
    #                     NoGo if EITHER endpoint meets criteria (OR)
    GoNogoProb <- sapply(seq(3), function(j) {
      if (j == 1) {
        # Go decision: (Go_prob1 >= gammaG1) AND (Go_prob2 >= gammaG2)
        I <- as.numeric((Go_prob1 >= gammaG1) & (Go_prob2 >= gammaG2) &
                          !((NoGo_prob1 >= gammaN1) | (NoGo_prob2 >= gammaN2)))
      } else if (j == 2) {
        # NoGo decision: (NoGo_prob1 >= gammaN1) OR (NoGo_prob2 >= gammaN2)
        I <- as.numeric(((NoGo_prob1 >= gammaN1) | (NoGo_prob2 >= gammaN2)) &
                          !((Go_prob1 >= gammaG1) & (Go_prob2 >= gammaG2)))
      } else {
        # Miss: Go and NoGo both satisfied
        I <- as.numeric((Go_prob1 >= gammaG1) & (Go_prob2 >= gammaG2) &
                          ((NoGo_prob1 >= gammaN1) | (NoGo_prob2 >= gammaN2)))
      }
      sum(I) / nsim
    })

  } else if (decision_criteria == 3) {
    # Criteria 3 (8)-(9): Go if EITHER endpoint meets criteria (OR)
    #                     NoGo if BOTH endpoints meet criteria (AND)
    GoNogoProb <- sapply(seq(3), function(j) {
      if (j == 1) {
        # Go decision: (Go_prob1 >= gammaG1) OR (Go_prob2 >= gammaG2)
        I <- as.numeric(((Go_prob1 >= gammaG1) | (Go_prob2 >= gammaG2)) &
                          !((NoGo_prob1 >= gammaN1) & (NoGo_prob2 >= gammaN2)))
      } else if (j == 2) {
        # NoGo decision: (NoGo_prob1 >= gammaN1) AND (NoGo_prob2 >= gammaN2)
        I <- as.numeric((NoGo_prob1 >= gammaN1) & (NoGo_prob2 >= gammaN2) &
                          !((Go_prob1 >= gammaG1) | (Go_prob2 >= gammaG2)))
      } else {
        # Miss: Go and NoGo both satisfied
        I <- as.numeric(((Go_prob1 >= gammaG1) | (Go_prob2 >= gammaG2)) &
                          (NoGo_prob1 >= gammaN1) & (NoGo_prob2 >= gammaN2))
      }
      sum(I) / nsim
    })

  } else if (decision_criteria == 4) {
    # Criteria 4 (10)-(11): Go if EITHER endpoint meets criteria (OR)
    #                       NoGo if EITHER endpoint meets criteria (OR)
    GoNogoProb <- sapply(seq(3), function(j) {
      if (j == 1) {
        # Go decision: (Go_prob1 >= gammaG1) OR (Go_prob2 >= gammaG2)
        I <- as.numeric(((Go_prob1 >= gammaG1) | (Go_prob2 >= gammaG2)) &
                          !((NoGo_prob1 >= gammaN1) | (NoGo_prob2 >= gammaN2)))
      } else if (j == 2) {
        # NoGo decision: (NoGo_prob1 >= gammaN1) OR (NoGo_prob2 >= gammaN2)
        I <- as.numeric(((NoGo_prob1 >= gammaN1) | (NoGo_prob2 >= gammaN2)) &
                          !((Go_prob1 >= gammaG1) | (Go_prob2 >= gammaG2)))
      } else {
        # Miss: Go and NoGo both satisfied
        I <- as.numeric(((Go_prob1 >= gammaG1) | (Go_prob2 >= gammaG2)) &
                          ((NoGo_prob1 >= gammaN1) | (NoGo_prob2 >= gammaN2)))
      }
      sum(I) / nsim
    })
  }

  # Check for positive Miss probabilities (indicates inappropriate thresholds)
  if (error_if_Miss) {
    if (sum(GoNogoProb[3]) > 0) {
      stop('Because positive Miss probability(s) is obtained, re-consider appropriate thresholds')
    }
  }

  # Calculate Miss probability (both Go and NoGo criteria met simultaneously)
  Miss <- GoNogoProb[3]

  # Calculate Gray probability (complement of Go and NoGo)
  if (Gray_inc_Miss) {
    # Include Miss in Gray probability
    GrayProb <- 1 - sum(GoNogoProb[-3])
  } else {
    # Exclude Miss from Gray probability
    GrayProb <- 1 - sum(GoNogoProb)
  }

  # Prepare results data frame
  results <- data.frame(
    pi1_00 = pi1[1], pi1_01 = pi1[2], pi1_10 = pi1[3], pi1_11 = pi1[4],
    pi2_00 = pi2[1], pi2_01 = pi2[2], pi2_10 = pi2[3], pi2_11 = pi2[4],
    Go = GoNogoProb[1],
    Gray = GrayProb,
    NoGo = GoNogoProb[2]
  )

  # Add Miss column when error_if_Miss is FALSE and Gray_inc_Miss is FALSE
  if (!error_if_Miss) {
    if (!Gray_inc_Miss) {
      results$Miss <- Miss
    }
  }

  # Address floating point error
  results[results < .Machine$double.eps ^ 0.25] <- 0

  return(results)
}
