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
#' @param theta.TV A numeric scalar giving the target value (TV) threshold used for
#'        the Go decision when \code{prob = 'posterior'}. Set to \code{NULL} when
#'        \code{prob = 'predictive'}.
#' @param theta.MAV A numeric scalar giving the minimum acceptable value (MAV)
#'        threshold used for the NoGo decision when \code{prob = 'posterior'}.
#'        Must satisfy \code{theta.TV > theta.MAV}. Set to \code{NULL} when
#'        \code{prob = 'predictive'}.
#' @param theta.NULL A numeric scalar giving the null hypothesis threshold used for
#'        both Go and NoGo decisions when \code{prob = 'predictive'}. Set to
#'        \code{NULL} when \code{prob = 'posterior'}.
#' @param gamma1 A numeric scalar in \code{(0, 1)} giving the minimum posterior/predictive
#'        probability required for a Go decision. Typically 0.8 or higher.
#' @param gamma2 A numeric scalar in \code{(0, 1)} giving the maximum posterior/predictive
#'        probability that triggers a NoGo decision. Must satisfy \code{gamma2 < gamma1}.
#' @param pi1 A numeric value or vector giving the true response probability(s) for
#'        group 1 (treatment) used to evaluate operating characteristics. Each element
#'        must be in \code{(0, 1)}.
#' @param pi2 A numeric value or vector giving the true response probability(s) for
#'        group 2 (control). For \code{design = 'uncontrolled'}, this parameter is
#'        not used in calculations but must be supplied; it is excluded from the output.
#'        When supplied as a vector, must have the same length as \code{pi1}.
#' @param n1 A positive integer giving the number of patients in group 1 in the
#'        proof-of-concept (PoC) trial.
#' @param n2 A positive integer giving the number of patients in group 2 in the
#'        PoC trial (also used as the hypothetical control size for uncontrolled design).
#' @param a1 A positive numeric scalar giving the first shape parameter (alpha) of the
#'        prior Beta distribution for group 1.
#' @param a2 A positive numeric scalar giving the first shape parameter (alpha) of the
#'        prior Beta distribution for group 2.
#' @param b1 A positive numeric scalar giving the second shape parameter (beta) of the
#'        prior Beta distribution for group 1.
#' @param b2 A positive numeric scalar giving the second shape parameter (beta) of the
#'        prior Beta distribution for group 2.
#' @param z A non-negative integer giving the hypothetical control responder count.
#'        Required when \code{design = 'uncontrolled'}; otherwise set to \code{NULL}.
#' @param m1 A positive integer giving the future sample size for group 1.
#'        Required when \code{prob = 'predictive'}; otherwise set to \code{NULL}.
#' @param m2 A positive integer giving the future sample size for group 2.
#'        Required when \code{prob = 'predictive'}; otherwise set to \code{NULL}.
#' @param ne1 A positive integer giving the number of patients in group 1 of the
#'        external data set. Required when \code{design = 'external'}; otherwise
#'        \code{NULL}.
#' @param ne2 A positive integer giving the number of patients in group 2 of the
#'        external data set. Required when \code{design = 'external'}; otherwise
#'        \code{NULL}.
#' @param ye1 A non-negative integer giving the number of responders in group 1 of
#'        the external data set. Required when \code{design = 'external'}; otherwise
#'        \code{NULL}.
#' @param ye2 A non-negative integer giving the number of responders in group 2 of
#'        the external data set. Required when \code{design = 'external'}; otherwise
#'        \code{NULL}.
#' @param ae1 A numeric scalar in \code{(0, 1]} giving the power prior weight for
#'        group 1. Required when \code{design = 'external'}; otherwise \code{NULL}.
#' @param ae2 A numeric scalar in \code{(0, 1]} giving the power prior weight for
#'        group 2. Required when \code{design = 'external'}; otherwise \code{NULL}.
#' @param error_if_Miss A logical scalar; if \code{TRUE} (default), the function stops
#'        with an error if Miss probability > 0, prompting reconsideration of thresholds.
#' @param Gray_inc_Miss A logical scalar; if \code{TRUE}, Miss probability is added
#'        to Gray probability. If \code{FALSE} (default), Miss is reported separately.
#'        Active only when \code{error_if_Miss = FALSE}.
#'
#' @return A data frame with one row per \code{pi1} scenario and columns:
#' \describe{
#'   \item{pi1}{True treatment response probability.}
#'   \item{pi2}{True control response probability (omitted for uncontrolled design).}
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
#'   \item All possible outcome pairs \eqn{(y_1, y_2)} with \eqn{y_1 \in \{0,\ldots,n_1\}}
#'         and \eqn{y_2 \in \{0,\ldots,n_2\}} (or fixed at \eqn{z} for uncontrolled) are
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
#' pbayesdecisionprob1bin(
#'   prob = 'posterior', design = 'uncontrolled',
#'   theta.TV = 0.30, theta.MAV = 0.15, theta.NULL = NULL,
#'   gamma1 = 0.75, gamma2 = 0.25,
#'   pi1 = c(0.3, 0.5, 0.7), pi2 = NULL,
#'   n1 = 15, n2 = 15,
#'   a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
#'   z = 5, m1 = NULL, m2 = NULL,
#'   ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE
#' )
#'
#' # Example 3: Posterior predictive probability for controlled design
#' pbayesdecisionprob1bin(
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
#' # Example 4: External design with 50 percent power prior borrowing
#' pbayesdecisionprob1bin(
#'   prob = 'posterior', design = 'external',
#'   theta.TV = 0.4, theta.MAV = 0.2, theta.NULL = NULL,
#'   gamma1 = 0.8, gamma2 = 0.2,
#'   pi1 = c(0.2, 0.4, 0.6, 0.8), pi2 = rep(0.2, 4),
#'   n1 = 12, n2 = 12,
#'   a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
#'   z = NULL, m1 = NULL, m2 = NULL,
#'   ne1 = 15, ne2 = 15, ye1 = 6, ye2 = 4, ae1 = 0.5, ae2 = 0.5,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE
#' )
#'
#' @importFrom stats dbinom
#' @export
pbayesdecisionprob1bin <- function(prob = 'posterior', design = 'controlled',
                                   theta.TV = NULL, theta.MAV = NULL,
                                   theta.NULL = NULL,
                                   gamma1, gamma2, pi1, pi2 = NULL,
                                   n1, n2, a1, a2, b1, b2,
                                   z = NULL, m1 = NULL, m2 = NULL,
                                   ne1 = NULL, ne2 = NULL,
                                   ye1 = NULL, ye2 = NULL,
                                   ae1 = NULL, ae2 = NULL,
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

  for (nm in c("n1", "n2")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
        val != floor(val) || val < 1L) {
      stop(paste0("'", nm, "' must be a single positive integer"))
    }
  }

  for (nm in c("a1", "a2", "b1", "b2")) {
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
  pi1 <- as.numeric(pi1)
  pi2 <- as.numeric(pi2)

  if (any(is.na(pi1)) || any(pi1 <= 0) || any(pi1 >= 1)) {
    stop("All elements of 'pi1' must be in (0, 1)")
  }

  if (design != 'uncontrolled') {
    if (any(is.na(pi2)) || any(pi2 <= 0) || any(pi2 >= 1)) {
      stop("All elements of 'pi2' must be in (0, 1)")
    }
    if (length(pi1) != length(pi2)) {
      stop("'pi1' and 'pi2' must have the same length")
    }
  }

  # Validate prob-specific parameters
  if (prob == 'posterior') {
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

  # Validate design-specific parameters
  if (design == 'uncontrolled') {
    if (is.null(z)) {
      stop("'z' must be non-NULL when design = 'uncontrolled'")
    }
    if (!is.numeric(z) || length(z) != 1L || is.na(z) ||
        z != floor(z) || z < 0L || z > n2) {
      stop("'z' must be a single non-negative integer not exceeding n2")
    }
  }

  if (design == 'external') {
    if (any(sapply(list(ne1, ne2, ye1, ye2, ae1, ae2), is.null))) {
      stop("'ne1', 'ne2', 'ye1', 'ye2', 'ae1', and 'ae2' must all be non-NULL when design = 'external'")
    }
    for (nm in c("ne1", "ne2")) {
      val <- get(nm)
      if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
          val != floor(val) || val < 1L) {
        stop(paste0("'", nm, "' must be a single positive integer"))
      }
    }
    for (nm in c("ye1", "ye2")) {
      val  <- get(nm)
      ne_v <- if (nm == "ye1") ne1 else ne2
      if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
          val != floor(val) || val < 0L || val > ne_v) {
        stop(paste0("'", nm, "' must be a single non-negative integer not exceeding the corresponding ne"))
      }
    }
    for (nm in c("ae1", "ae2")) {
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
    theta0 <- c(theta.TV, theta.MAV)
  } else {
    # Both thresholds equal to theta.NULL for predictive probability
    theta0 <- c(theta.NULL, theta.NULL)
  }

  # --- Enumerate all possible outcome combinations ---
  Y1 <- 0:n1
  if (design == 'uncontrolled') {
    all_y1 <- Y1
    all_y2 <- rep(z, length(Y1))
  } else {
    grid   <- expand.grid(y1 = Y1, y2 = 0:n2)
    all_y1 <- grid$y1
    all_y2 <- grid$y2
  }

  # --- Compute posterior/predictive probabilities for all outcome pairs ---
  # Use z argument only for uncontrolled design; pass NULL otherwise
  z_arg <- if (design == 'uncontrolled') z else NULL

  gPost_Go <- pbayespostpred1bin(
    prob = prob, design = design, theta0 = theta0[1],
    n1 = n1, n2 = n2, y1 = all_y1, y2 = if (design == 'uncontrolled') NULL else all_y2,
    a1 = a1, a2 = a2, b1 = b1, b2 = b2,
    m1 = m1, m2 = m2, z = z_arg,
    ne1 = ne1, ne2 = ne2, ye1 = ye1, ye2 = ye2, ae1 = ae1, ae2 = ae2,
    lower.tail = FALSE
  )

  gPost_NoGo <- pbayespostpred1bin(
    prob = prob, design = design, theta0 = theta0[2],
    n1 = n1, n2 = n2, y1 = all_y1, y2 = if (design == 'uncontrolled') NULL else all_y2,
    a1 = a1, a2 = a2, b1 = b1, b2 = b2,
    m1 = m1, m2 = m2, z = z_arg,
    ne1 = ne1, ne2 = ne2, ye1 = ye1, ye2 = ye2, ae1 = ae1, ae2 = ae2,
    lower.tail = TRUE
  )

  # --- Decision indicators (mutually exclusive: Go, NoGo, Miss; Gray = complement) ---
  probs_Go   <- (gPost_Go  >= gamma1) & (gPost_NoGo <  gamma2)
  probs_NoGo <- (gPost_Go  <  gamma1) & (gPost_NoGo >= gamma2)
  probs_Miss <- (gPost_Go  >= gamma1) & (gPost_NoGo >= gamma2)

  # --- Weight outcomes by binomial probability under each scenario ---
  n_scenarios <- length(pi1)
  GoNogoProb  <- matrix(0, nrow = n_scenarios, ncol = 3L)

  for (scenario in seq_len(n_scenarios)) {
    if (design == 'uncontrolled') {
      # y2 is fixed at z; sum over y1 only
      binom_w <- dbinom(all_y1, n1, pi1[scenario])
      GoNogoProb[scenario, 1L] <- sum(probs_Go   * binom_w)
      GoNogoProb[scenario, 2L] <- sum(probs_NoGo * binom_w)
      GoNogoProb[scenario, 3L] <- sum(probs_Miss * binom_w)
    } else {
      # Sum over all (y1, y2) combinations
      w <- dbinom(all_y1, n1, pi1[scenario]) * dbinom(all_y2, n2, pi2[scenario])
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
      pi1  = pi1,
      Go   = GoNogoProb[, 1L],
      Gray = GrayProb,
      NoGo = GoNogoProb[, 2L]
    )
  } else {
    results <- data.frame(
      pi1  = pi1,
      pi2  = pi2,
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
  attr(results, 'prob')          <- prob
  attr(results, 'design')        <- design
  attr(results, 'theta.TV')      <- theta.TV
  attr(results, 'theta.MAV')     <- theta.MAV
  attr(results, 'theta.NULL')    <- theta.NULL
  attr(results, 'gamma1')        <- gamma1
  attr(results, 'gamma2')        <- gamma2
  attr(results, 'n1')            <- n1
  attr(results, 'n2')            <- n2
  attr(results, 'a1')            <- a1
  attr(results, 'a2')            <- a2
  attr(results, 'b1')            <- b1
  attr(results, 'b2')            <- b2
  attr(results, 'z')             <- z
  attr(results, 'm1')            <- m1
  attr(results, 'm2')            <- m2
  attr(results, 'ne1')           <- ne1
  attr(results, 'ne2')           <- ne2
  attr(results, 'ye1')           <- ye1
  attr(results, 'ye2')           <- ye2
  attr(results, 'ae1')           <- ae1
  attr(results, 'ae2')           <- ae2
  attr(results, 'error_if_Miss') <- error_if_Miss
  attr(results, 'Gray_inc_Miss') <- Gray_inc_Miss

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
  prob          <- attr(x, 'prob')
  design        <- attr(x, 'design')
  gamma1        <- attr(x, 'gamma1')
  gamma2        <- attr(x, 'gamma2')
  n1            <- attr(x, 'n1')
  n2            <- attr(x, 'n2')
  a1            <- attr(x, 'a1')
  a2            <- attr(x, 'a2')
  b1            <- attr(x, 'b1')
  b2            <- attr(x, 'b2')
  z             <- attr(x, 'z')
  m1            <- attr(x, 'm1')
  m2            <- attr(x, 'm2')
  ne1           <- attr(x, 'ne1')
  ne2           <- attr(x, 'ne2')
  ye1           <- attr(x, 'ye1')
  ye2           <- attr(x, 'ye2')
  ae1           <- attr(x, 'ae1')
  ae2           <- attr(x, 'ae2')
  error_if_Miss <- attr(x, 'error_if_Miss')
  Gray_inc_Miss <- attr(x, 'Gray_inc_Miss')

  # Build threshold string based on probability type
  if(prob == 'posterior') {
    theta_str <- sprintf('TV = %s, MAV = %s',
                         fmt(attr(x, 'theta.TV')), fmt(attr(x, 'theta.MAV')))
  } else {
    theta_str <- sprintf('NULL = %s', fmt(attr(x, 'theta.NULL')))
  }

  # Prior label depends on probability type
  prior_label <- if(prob == 'posterior') 'Prior (Beta)     ' else 'Prior (Beta-bin) '

  # Print header
  cat('Go/NoGo/Gray Decision Probabilities (Single Binary Endpoint)\n')
  cat(strrep('-', 60), '\n')
  cat(sprintf('  Probability type : %s\n',   prob))
  cat(sprintf('  Design           : %s\n',   design))
  cat(sprintf('  Threshold(s)     : %s\n',   theta_str))
  cat(sprintf('  Go  threshold    : gamma1 = %s\n', fmt(gamma1)))
  cat(sprintf('  NoGo threshold   : gamma2 = %s\n', fmt(gamma2)))
  cat(sprintf('  Sample size      : n1 = %s, n2 = %s\n', fmt(n1), fmt(n2)))
  cat(sprintf('  %s: a1 = %s, a2 = %s, b1 = %s, b2 = %s\n',
              prior_label, fmt(a1), fmt(a2), fmt(b1), fmt(b2)))
  if(design == 'uncontrolled') {
    cat(sprintf('  Uncontrolled     : z = %s\n', fmt(z)))
  }
  if(prob == 'predictive') {
    cat(sprintf('  Future trial     : m1 = %s, m2 = %s\n', fmt(m1), fmt(m2)))
  }
  if(design == 'external') {
    cat(sprintf('  External data    : ne1 = %s, ne2 = %s, ye1 = %s, ye2 = %s, ae1 = %s, ae2 = %s\n',
                fmt(ne1), fmt(ne2), fmt(ye1), fmt(ye2), fmt(ae1), fmt(ae2)))
  }
  cat(sprintf('  Miss handling    : error_if_Miss = %s, Gray_inc_Miss = %s\n',
              fmt(error_if_Miss), fmt(Gray_inc_Miss)))
  cat(strrep('-', 60), '\n')

  # Format numeric columns (probability columns only, not pi1/pi2)
  prob_cols <- names(x)[!names(x) %in% c('pi1', 'pi2')]
  x_print <- x
  x_print[prob_cols] <- lapply(x[prob_cols], function(col) {
    formatC(col, digits = digits, format = 'f')
  })

  # Print table without row names (call print.data.frame explicitly to avoid recursion)
  print.data.frame(x_print, row.names = FALSE, quote = FALSE)
  cat(strrep('-', 60), '\n')

  invisible(x)
}
