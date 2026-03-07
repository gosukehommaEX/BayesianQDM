#' Go/NoGo/Gray Decision Probabilities for a Clinical Trial with a Single Binary Endpoint
#'
#' Evaluates operating characteristics (Go, NoGo, Gray probabilities) for binary-outcome
#' clinical trials under the Bayesian framework by enumerating all possible trial
#' outcomes. The function supports controlled, uncontrolled, and external designs.
#'
#' @param prob A character string specifying the probability type for decision-making.
#'        Must be \code{'posterior'} or \code{'predictive'}.
#' @param design A character string specifying the trial design.
#'        Must be \code{'controlled'}, \code{'uncontrolled'}, or \code{'external'}.
#' @param theta_TV A numeric scalar giving the target value (TV) threshold used for
#'        the Go decision when \code{prob = 'posterior'}. Set to \code{NULL} when
#'        \code{prob = 'predictive'}.
#' @param theta_MAV A numeric scalar giving the minimum acceptable value (MAV)
#'        threshold used for the NoGo decision when \code{prob = 'posterior'}.
#'        Must satisfy \code{theta_TV > theta_MAV}. Set to \code{NULL} when
#'        \code{prob = 'predictive'}.
#' @param theta_NULL A numeric scalar giving the null hypothesis threshold used for
#'        both Go and NoGo decisions when \code{prob = 'predictive'}. Set to
#'        \code{NULL} when \code{prob = 'posterior'}.
#' @param gamma_go A numeric scalar in \code{(0, 1)} giving the minimum posterior or
#'        predictive probability required for a Go decision. Typically 0.8 or higher.
#' @param gamma_nogo A numeric scalar in \code{(0, 1)} giving the maximum posterior or
#'        predictive probability that triggers a NoGo decision. Must satisfy
#'        \code{gamma_nogo < gamma_go}.
#' @param pi_t A numeric value or vector giving the true response probability(s) for
#'        the treatment group used to evaluate operating characteristics. Each element
#'        must be in \code{(0, 1)}.
#' @param pi_c A numeric value or vector giving the true response probability(s) for
#'        the control group. For \code{design = 'uncontrolled'}, this parameter is
#'        not used in calculations but must be supplied; it is excluded from the output.
#'        When supplied as a vector, must have the same length as \code{pi_t}.
#' @param n_t A positive integer giving the number of patients in the treatment group in the
#'        proof-of-concept (PoC) trial.
#' @param n_c A positive integer giving the number of patients in the control group in the
#'        PoC trial (also used as the hypothetical control size for uncontrolled design).
#' @param a_t A positive numeric scalar giving the first shape parameter (alpha) of the
#'        prior Beta distribution for the treatment group.
#' @param a_c A positive numeric scalar giving the first shape parameter (alpha) of the
#'        prior Beta distribution for the control group.
#' @param b_t A positive numeric scalar giving the second shape parameter (beta) of the
#'        prior Beta distribution for the treatment group.
#' @param b_c A positive numeric scalar giving the second shape parameter (beta) of the
#'        prior Beta distribution for the control group.
#' @param z A non-negative integer giving the hypothetical control responder count.
#'        Required when \code{design = 'uncontrolled'}; otherwise set to \code{NULL}.
#' @param m_t A positive integer giving the future sample size for the treatment group.
#'        Required when \code{prob = 'predictive'}; otherwise set to \code{NULL}.
#' @param m_c A positive integer giving the future sample size for the control group.
#'        Required when \code{prob = 'predictive'}; otherwise set to \code{NULL}.
#' @param ne_t A positive integer giving the number of patients in the treatment group
#'        of the external data set. Required when \code{design = 'external'}; otherwise
#'        \code{NULL}.
#' @param ne_c A positive integer giving the number of patients in the control group
#'        of the external data set. Required when \code{design = 'external'}; otherwise
#'        \code{NULL}.
#' @param ye_t A non-negative integer giving the number of responders in the
#'        treatment group of the external data set. Required when
#'        \code{design = 'external'}; otherwise \code{NULL}.
#' @param ye_c A non-negative integer giving the number of responders in the
#'        control group of the external data set. Required when
#'        \code{design = 'external'}; otherwise \code{NULL}.
#' @param ae_t A numeric scalar in \code{(0, 1]} giving the power prior weight for
#'        the treatment group. Required when \code{design = 'external'};
#'        otherwise \code{NULL}.
#' @param ae_c A numeric scalar in \code{(0, 1]} giving the power prior weight for
#'        the control group. Required when \code{design = 'external'};
#'        otherwise \code{NULL}.
#' @param error_if_Miss A logical scalar; if \code{TRUE} (default), the function stops
#'        with an error if Miss probability > 0, prompting reconsideration of thresholds.
#' @param Gray_inc_Miss A logical scalar; if \code{TRUE}, Miss probability is added
#'        to Gray probability. If \code{FALSE} (default), Miss is reported separately.
#'        Active only when \code{error_if_Miss = FALSE}.
#'
#' @return A data frame with one row per \code{pi_t} scenario and columns:
#' \describe{
#'   \item{pi_t}{True treatment response probability.}
#'   \item{pi_c}{True control response probability (omitted for uncontrolled design).}
#'   \item{Go}{Probability of making a Go decision.}
#'   \item{Gray}{Probability of making a Gray (inconclusive) decision.}
#'   \item{NoGo}{Probability of making a NoGo decision.}
#'   \item{Miss}{(Optional) Probability where Go and NoGo criteria are simultaneously
#'               met. Included when \code{error_if_Miss = FALSE} and
#'               \code{Gray_inc_Miss = FALSE}.}
#' }
#' The returned object has S3 class \code{pbayesdecisionprob1bin} with an associated
#' \code{print} method.
#'
#' @details
#' Operating characteristics are computed by exact enumeration:
#' \enumerate{
#'   \item All possible outcome pairs \eqn{(y_t, y_c)} with \eqn{y_t \in \{0,\ldots,n_t\}}
#'         and \eqn{y_c \in \{0,\ldots,n_c\}} (or fixed at \eqn{z} for uncontrolled) are
#'         evaluated.
#'   \item For each pair, \code{\link{pbayespostpred1bin}} computes the posterior or predictive
#'         probability at both thresholds (TV/MAV or NULL).
#'   \item Outcomes are classified into Go, NoGo, Miss, or Gray:
#'         \itemize{
#'           \item \strong{Go}: \eqn{P(\mathrm{Go}) \ge \gamma_1} AND
#'                 \eqn{P(\mathrm{NoGo}) < \gamma_2}
#'           \item \strong{NoGo}: \eqn{P(\mathrm{Go}) < \gamma_1} AND
#'                 \eqn{P(\mathrm{NoGo}) \ge \gamma_2}
#'           \item \strong{Miss}: both Go and NoGo criteria met simultaneously
#'           \item \strong{Gray}: neither Go nor NoGo criteria met
#'         }
#'   \item Each outcome is weighted by its binomial probability under the true rates.
#' }
#'
#' @examples
#' # Example 1: Controlled design with posterior probability
#' pbayesdecisionprob1bin(
#'   prob = 'posterior', design = 'controlled',
#'   theta_TV = 0.4, theta_MAV = 0.2, theta_NULL = NULL,
#'   gamma_go = 0.8, gamma_nogo = 0.2,
#'   pi_t = c(0.2, 0.4, 0.6, 0.8), pi_c = rep(0.2, 4),
#'   n_t = 12, n_c = 12,
#'   a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
#'   z = NULL, m_t = NULL, m_c = NULL,
#'   ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL, ae_t = NULL, ae_c = NULL,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE
#' )
#'
#' # Example 2: Uncontrolled design with hypothetical control
#' pbayesdecisionprob1bin(
#'   prob = 'posterior', design = 'uncontrolled',
#'   theta_TV = 0.30, theta_MAV = 0.15, theta_NULL = NULL,
#'   gamma_go = 0.75, gamma_nogo = 0.25,
#'   pi_t = c(0.3, 0.5, 0.7), pi_c = NULL,
#'   n_t = 15, n_c = 15,
#'   a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
#'   z = 5, m_t = NULL, m_c = NULL,
#'   ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL, ae_t = NULL, ae_c = NULL,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE
#' )
#'
#' # Example 3: External design with 50 percent power prior borrowing
#' pbayesdecisionprob1bin(
#'   prob = 'posterior', design = 'external',
#'   theta_TV = 0.4, theta_MAV = 0.2, theta_NULL = NULL,
#'   gamma_go = 0.8, gamma_nogo = 0.2,
#'   pi_t = c(0.2, 0.4, 0.6, 0.8), pi_c = rep(0.2, 4),
#'   n_t = 12, n_c = 12,
#'   a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
#'   z = NULL, m_t = NULL, m_c = NULL,
#'   ne_t = 15, ne_c = 15, ye_t = 6, ye_c = 4, ae_t = 0.5, ae_c = 0.5,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE
#' )
#'
#' # Example 4: Posterior predictive probability for controlled design
#' pbayesdecisionprob1bin(
#'   prob = 'predictive', design = 'controlled',
#'   theta_TV = NULL, theta_MAV = NULL, theta_NULL = 0,
#'   gamma_go = 0.9, gamma_nogo = 0.3,
#'   pi_t = c(0.2, 0.4, 0.6, 0.8), pi_c = rep(0.2, 4),
#'   n_t = 12, n_c = 12,
#'   a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
#'   z = NULL, m_t = 30, m_c = 30,
#'   ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL, ae_t = NULL, ae_c = NULL,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE
#' )
#'
#' # Example 5: Uncontrolled design with posterior predictive probability
#' pbayesdecisionprob1bin(
#'   prob = 'predictive', design = 'uncontrolled',
#'   theta_TV = NULL, theta_MAV = NULL, theta_NULL = 0,
#'   gamma_go = 0.75, gamma_nogo = 0.25,
#'   pi_t = c(0.3, 0.5, 0.7), pi_c = NULL,
#'   n_t = 15, n_c = 15,
#'   a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
#'   z = 5, m_t = 30, m_c = 30,
#'   ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL, ae_t = NULL, ae_c = NULL,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE
#' )
#'
#' # Example 6: External design with posterior predictive probability
#' pbayesdecisionprob1bin(
#'   prob = 'predictive', design = 'external',
#'   theta_TV = NULL, theta_MAV = NULL, theta_NULL = 0,
#'   gamma_go = 0.9, gamma_nogo = 0.3,
#'   pi_t = c(0.2, 0.4, 0.6, 0.8), pi_c = rep(0.2, 4),
#'   n_t = 12, n_c = 12,
#'   a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
#'   z = NULL, m_t = 30, m_c = 30,
#'   ne_t = 15, ne_c = 15, ye_t = 6, ye_c = 4, ae_t = 0.5, ae_c = 0.5,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE
#' )
#'
#' @importFrom stats dbinom
#' @export
pbayesdecisionprob1bin <- function(prob = 'posterior', design = 'controlled',
                                   theta_TV = NULL, theta_MAV = NULL,
                                   theta_NULL = NULL,
                                   gamma_go, gamma_nogo, pi_t, pi_c = NULL,
                                   n_t, n_c, a_t, a_c, b_t, b_c,
                                   z = NULL, m_t = NULL, m_c = NULL,
                                   ne_t = NULL, ne_c = NULL,
                                   ye_t = NULL, ye_c = NULL,
                                   ae_t = NULL, ae_c = NULL,
                                   error_if_Miss = TRUE,
                                   Gray_inc_Miss = FALSE) {

  # --- Input validation ---
  if (!is.character(prob) || length(prob) != 1L ||
      !prob %in% c('posterior', 'predictive')) {
    stop("'prob' must be either 'posterior' or 'predictive'")
  }

  if (!is.character(design) || length(design) != 1L ||
      !design %in% c('controlled', 'uncontrolled', 'external')) {
    stop("'design' must be 'controlled', 'uncontrolled', or 'external'")
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

  for (nm in c("n_t", "n_c")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
        val != floor(val) || val < 1L) {
      stop(paste0("'", nm, "' must be a single positive integer"))
    }
  }

  for (nm in c("a_t", "a_c", "b_t", "b_c")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) || val <= 0) {
      stop(paste0("'", nm, "' must be a single positive numeric value"))
    }
  }

  if (!is.logical(error_if_Miss) || length(error_if_Miss) != 1L ||
      is.na(error_if_Miss)) {
    stop("'error_if_Miss' must be a single logical value (TRUE or FALSE)")
  }

  if (!is.logical(Gray_inc_Miss) || length(Gray_inc_Miss) != 1L ||
      is.na(Gray_inc_Miss)) {
    stop("'Gray_inc_Miss' must be a single logical value (TRUE or FALSE)")
  }

  # Validate pi vectors
  pi_t <- as.numeric(pi_t)
  pi_c <- as.numeric(pi_c)

  if (any(is.na(pi_t)) || any(pi_t <= 0) || any(pi_t >= 1)) {
    stop("All elements of 'pi_t' must be in (0, 1)")
  }

  if (design != 'uncontrolled') {
    if (any(is.na(pi_c)) || any(pi_c <= 0) || any(pi_c >= 1)) {
      stop("All elements of 'pi_c' must be in (0, 1)")
    }
    if (length(pi_t) != length(pi_c)) {
      stop("'pi_t' and 'pi_c' must have the same length")
    }
  }

  # Validate prob-specific parameters
  if (prob == 'posterior') {
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

  # Validate design-specific parameters
  if (design == 'uncontrolled') {
    if (is.null(z)) {
      stop("'z' must be non-NULL when design = 'uncontrolled'")
    }
    if (!is.numeric(z) || length(z) != 1L || is.na(z) ||
        z != floor(z) || z < 0L || z > n_c) {
      stop("'z' must be a single non-negative integer not exceeding n_c")
    }
  }

  if (design == 'external') {
    if (any(sapply(list(ne_t, ne_c, ye_t, ye_c, ae_t, ae_c), is.null))) {
      stop("'ne_t', 'ne_c', 'ye_t', 'ye_c', 'ae_t', and 'ae_c' must all be non-NULL when design = 'external'")
    }
    for (nm in c("ne_t", "ne_c")) {
      val <- get(nm)
      if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
          val != floor(val) || val < 1L) {
        stop(paste0("'", nm, "' must be a single positive integer"))
      }
    }
    for (nm in c("ye_t", "ye_c")) {
      val  <- get(nm)
      ne_v <- if (nm == "ye_t") ne_t else ne_c
      if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
          val != floor(val) || val < 0L || val > ne_v) {
        stop(paste0("'", nm, "' must be a single non-negative integer not exceeding the corresponding ne"))
      }
    }
    for (nm in c("ae_t", "ae_c")) {
      val <- get(nm)
      if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
          val <= 0 || val > 1) {
        stop(paste0("'", nm, "' must be a single numeric value in (0, 1]"))
      }
    }
  }

  # --- Set decision thresholds ---
  if (prob == 'posterior') {
    # Go threshold: TV, NoGo threshold: MAV
    theta0 <- c(theta_TV, theta_MAV)
  } else {
    # Both thresholds equal to theta_NULL for predictive probability
    theta0 <- c(theta_NULL, theta_NULL)
  }

  # --- Enumerate all possible outcome combinations ---
  Y_t <- 0:n_t
  if (design == 'uncontrolled') {
    all_y_t <- Y_t
    all_y_c <- rep(z, length(Y_t))
  } else {
    grid   <- expand.grid(y_t = Y_t, y_c = 0:n_c)
    all_y_t <- grid$y_t
    all_y_c <- grid$y_c
  }

  # --- Compute posterior/predictive probabilities for all outcome pairs ---
  # Use z argument only for uncontrolled design; pass NULL otherwise
  z_arg <- if (design == 'uncontrolled') z else NULL

  gPost_Go <- pbayespostpred1bin(
    prob = prob, design = design, theta0 = theta0[1],
    n_t = n_t, n_c = n_c, y_t = all_y_t, y_c = if (design == 'uncontrolled') NULL else all_y_c,
    a_t = a_t, a_c = a_c, b_t = b_t, b_c = b_c,
    m_t = m_t, m_c = m_c, z = z_arg,
    ne_t = ne_t, ne_c = ne_c, ye_t = ye_t, ye_c = ye_c, ae_t = ae_t, ae_c = ae_c,
    lower.tail = FALSE
  )

  gPost_NoGo <- pbayespostpred1bin(
    prob = prob, design = design, theta0 = theta0[2],
    n_t = n_t, n_c = n_c, y_t = all_y_t, y_c = if (design == 'uncontrolled') NULL else all_y_c,
    a_t = a_t, a_c = a_c, b_t = b_t, b_c = b_c,
    m_t = m_t, m_c = m_c, z = z_arg,
    ne_t = ne_t, ne_c = ne_c, ye_t = ye_t, ye_c = ye_c, ae_t = ae_t, ae_c = ae_c,
    lower.tail = TRUE
  )

  # --- Decision indicators (mutually exclusive: Go, NoGo, Miss; Gray = complement) ---
  probs_Go   <- (gPost_Go  >= gamma_go) & (gPost_NoGo <  gamma_nogo)
  probs_NoGo <- (gPost_Go  <  gamma_go) & (gPost_NoGo >= gamma_nogo)
  probs_Miss <- (gPost_Go  >= gamma_go) & (gPost_NoGo >= gamma_nogo)

  # --- Weight outcomes by binomial probability under each scenario ---
  n_scenarios <- length(pi_t)
  GoNogoProb  <- matrix(0, nrow = n_scenarios, ncol = 3L)

  for (scenario in seq_len(n_scenarios)) {
    if (design == 'uncontrolled') {
      # y_c is fixed at z; sum over y_t only
      binom_w <- dbinom(all_y_t, n_t, pi_t[scenario])
      GoNogoProb[scenario, 1L] <- sum(probs_Go   * binom_w)
      GoNogoProb[scenario, 2L] <- sum(probs_NoGo * binom_w)
      GoNogoProb[scenario, 3L] <- sum(probs_Miss * binom_w)
    } else {
      # Sum over all (y_t, y_c) combinations
      w <- dbinom(all_y_t, n_t, pi_t[scenario]) * dbinom(all_y_c, n_c, pi_c[scenario])
      GoNogoProb[scenario, 1L] <- sum(probs_Go   * w)
      GoNogoProb[scenario, 2L] <- sum(probs_NoGo * w)
      GoNogoProb[scenario, 3L] <- sum(probs_Miss * w)
    }
  }

  # --- Check for positive Miss probability ---
  if (error_if_Miss && sum(GoNogoProb[, 3L]) > 0) {
    stop("Positive Miss probability detected. Please re-consider the chosen thresholds.")
  }

  # --- Gray probability (complement of Go and NoGo) ---
  if (Gray_inc_Miss) {
    # Include Miss in Gray: Gray = 1 - Go - NoGo
    GrayProb <- 1 - rowSums(GoNogoProb[, -3L, drop = FALSE])
  } else {
    # Exclude Miss from Gray: Gray = 1 - Go - NoGo - Miss
    GrayProb <- 1 - rowSums(GoNogoProb[, , drop = FALSE])
  }

  # --- Build results data frame ---
  if (design == 'uncontrolled') {
    results <- data.frame(
      pi_t  = pi_t,
      Go   = GoNogoProb[, 1L],
      Gray = GrayProb,
      NoGo = GoNogoProb[, 2L]
    )
  } else {
    results <- data.frame(
      pi_t  = pi_t,
      pi_c  = pi_c,
      Go   = GoNogoProb[, 1L],
      Gray = GrayProb,
      NoGo = GoNogoProb[, 2L]
    )
  }

  if (!error_if_Miss && !Gray_inc_Miss) {
    results$Miss <- GoNogoProb[, 3L]
  }

  # Address floating-point rounding artefacts
  results[results < .Machine$double.eps ^ 0.25] <- 0

  # Attach metadata as attributes for use in print()
  attr(results, 'prob')           <- prob
  attr(results, 'design')         <- design
  attr(results, 'theta_TV')       <- theta_TV
  attr(results, 'theta_MAV')      <- theta_MAV
  attr(results, 'theta_NULL')     <- theta_NULL
  attr(results, 'gamma_go')       <- gamma_go
  attr(results, 'gamma_nogo')     <- gamma_nogo
  attr(results, 'n_t')            <- n_t
  attr(results, 'n_c')            <- n_c
  attr(results, 'a_t')            <- a_t
  attr(results, 'a_c')            <- a_c
  attr(results, 'b_t')            <- b_t
  attr(results, 'b_c')            <- b_c
  attr(results, 'z')              <- z
  attr(results, 'm_t')            <- m_t
  attr(results, 'm_c')            <- m_c
  attr(results, 'ne_t')           <- ne_t
  attr(results, 'ne_c')           <- ne_c
  attr(results, 'ye_t')           <- ye_t
  attr(results, 'ye_c')           <- ye_c
  attr(results, 'ae_t')           <- ae_t
  attr(results, 'ae_c')           <- ae_c
  attr(results, 'error_if_Miss')  <- error_if_Miss
  attr(results, 'Gray_inc_Miss')  <- Gray_inc_Miss

  # Assign S3 class
  class(results) <- c('pbayesdecisionprob1bin', 'data.frame')

  return(results)
}

#' Print Method for pbayesdecisionprob1bin Objects
#'
#' Displays a formatted summary of Go/NoGo/Gray decision probabilities
#' for binary endpoint results returned by \code{\link{pbayesdecisionprob1bin}}.
#'
#' @param x An object of class \code{pbayesdecisionprob1bin}.
#' @param digits A positive integer specifying the number of decimal places
#'        for probability values. Default is 4.
#' @param ... Further arguments passed to or from other methods (ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.pbayesdecisionprob1bin <- function(x, digits = 4, ...) {
  # Helper to format a value as string (NULL -> "NULL")
  fmt <- function(v) if(is.null(v)) 'NULL' else as.character(v)

  # Extract metadata from attributes
  prob           <- attr(x, 'prob')
  design         <- attr(x, 'design')
  gamma_go       <- attr(x, 'gamma_go')
  gamma_nogo     <- attr(x, 'gamma_nogo')
  n_t            <- attr(x, 'n_t')
  n_c            <- attr(x, 'n_c')
  a_t            <- attr(x, 'a_t')
  a_c            <- attr(x, 'a_c')
  b_t            <- attr(x, 'b_t')
  b_c            <- attr(x, 'b_c')
  z              <- attr(x, 'z')
  m_t            <- attr(x, 'm_t')
  m_c            <- attr(x, 'm_c')
  ne_t           <- attr(x, 'ne_t')
  ne_c           <- attr(x, 'ne_c')
  ye_t           <- attr(x, 'ye_t')
  ye_c           <- attr(x, 'ye_c')
  ae_t           <- attr(x, 'ae_t')
  ae_c           <- attr(x, 'ae_c')
  error_if_Miss  <- attr(x, 'error_if_Miss')
  Gray_inc_Miss  <- attr(x, 'Gray_inc_Miss')

  # Build threshold string based on probability type
  if(prob == 'posterior') {
    theta_str <- sprintf('TV = %s, MAV = %s',
                         fmt(attr(x, 'theta_TV')), fmt(attr(x, 'theta_MAV')))
  } else {
    theta_str <- sprintf('NULL = %s', fmt(attr(x, 'theta_NULL')))
  }

  # Prior label depends on probability type
  prior_label <- if(prob == 'posterior') 'Prior (Beta)     ' else 'Prior (Beta-bin) '

  # Print header
  cat('Go/NoGo/Gray Decision Probabilities (Single Binary Endpoint)\n')
  cat(strrep('-', 60), '\n')
  cat(sprintf('  Probability type : %s\n',   prob))
  cat(sprintf('  Design           : %s\n',   design))
  cat(sprintf('  Threshold(s)     : %s\n',   theta_str))
  cat(sprintf('  Go  threshold    : gamma_go = %s\n', fmt(gamma_go)))
  cat(sprintf('  NoGo threshold   : gamma_nogo = %s\n', fmt(gamma_nogo)))
  cat(sprintf('  Sample size      : n_t = %s, n_c = %s\n', fmt(n_t), fmt(n_c)))
  cat(sprintf('  %s: a_t = %s, a_c = %s, b_t = %s, b_c = %s\n',
              prior_label, fmt(a_t), fmt(a_c), fmt(b_t), fmt(b_c)))
  if(design == 'uncontrolled') {
    cat(sprintf('  Uncontrolled     : z = %s\n', fmt(z)))
  }
  if(prob == 'predictive') {
    cat(sprintf('  Future trial     : m_t = %s, m_c = %s\n', fmt(m_t), fmt(m_c)))
  }
  if(design == 'external') {
    cat(sprintf('  External data    : ne_t = %s, ne_c = %s, ye_t = %s, ye_c = %s, ae_t = %s, ae_c = %s\n',
                fmt(ne_t), fmt(ne_c), fmt(ye_t), fmt(ye_c), fmt(ae_t), fmt(ae_c)))
  }
  cat(sprintf('  Miss handling    : error_if_Miss = %s, Gray_inc_Miss = %s\n',
              fmt(error_if_Miss), fmt(Gray_inc_Miss)))
  cat(strrep('-', 60), '\n')

  # Format numeric columns (probability columns only, not pi_t/pi_c)
  prob_cols <- names(x)[!names(x) %in% c('pi_t', 'pi_c')]
  x_print <- x
  x_print[prob_cols] <- lapply(x[prob_cols], function(col) {
    formatC(col, digits = digits, format = 'f')
  })

  # Print table without row names (call print.data.frame explicitly to avoid recursion)
  print.data.frame(x_print, row.names = FALSE, quote = FALSE)
  cat(strrep('-', 60), '\n')

  invisible(x)
}
