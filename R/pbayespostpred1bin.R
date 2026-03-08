#' Bayesian Posterior or Posterior Predictive Probability for a Single
#' Binary Endpoint
#'
#' Computes the Bayesian posterior probability or posterior predictive
#' probability for binary-outcome clinical trials under a beta-binomial
#' conjugate model. The function supports controlled, uncontrolled, and
#' external designs, with optional incorporation of external data through
#' power priors. Vector inputs for \code{y_t} and \code{y_c} are supported
#' for efficient batch processing (e.g., across all possible trial outcomes
#' in \code{\link{pbayesdecisionprob1bin}}).
#'
#' @param prob A character string specifying the probability type.
#'        Must be \code{'posterior'} or \code{'predictive'}.
#' @param design A character string specifying the trial design.
#'        Must be \code{'controlled'}, \code{'uncontrolled'}, or
#'        \code{'external'}.
#' @param theta0 A numeric scalar in \code{(-1, 1)} giving the pre-specified
#'        threshold for the treatment effect (difference in response rates).
#' @param n_t A positive integer giving the number of patients in the
#'        treatment group in the proof-of-concept (PoC) trial.
#' @param n_c A positive integer giving the number of patients in the
#'        control group in the PoC trial. For \code{design = 'uncontrolled'},
#'        this is the hypothetical control sample size (required for
#'        consistency with other designs).
#' @param y_t A non-negative integer or integer vector giving the observed
#'        number of responders in the treatment group (must satisfy
#'        \code{0 <= y_t <= n_t}).
#' @param y_c A non-negative integer or integer vector giving the number of
#'        responders in the control group (must satisfy
#'        \code{0 <= y_c <= n_c}). Set to \code{NULL} for
#'        \code{design = 'uncontrolled'} and use \code{z} instead. When
#'        provided as a vector, must have the same length as \code{y_t}.
#' @param a_t A positive numeric scalar giving the first shape parameter
#'        (alpha) of the prior Beta distribution for the treatment group.
#' @param a_c A positive numeric scalar giving the first shape parameter
#'        (alpha) of the prior Beta distribution for the control group.
#' @param b_t A positive numeric scalar giving the second shape parameter
#'        (beta) of the prior Beta distribution for the treatment group.
#' @param b_c A positive numeric scalar giving the second shape parameter
#'        (beta) of the prior Beta distribution for the control group.
#' @param m_t A positive integer giving the number of patients in the
#'        treatment group for the future trial. Required when
#'        \code{prob = 'predictive'}; otherwise set to \code{NULL}.
#' @param m_c A positive integer giving the number of patients in the
#'        control group for the future trial. Required when
#'        \code{prob = 'predictive'}; otherwise set to \code{NULL}.
#' @param z A non-negative integer giving the hypothetical number of responders
#'        in the control group. Required when \code{design = 'uncontrolled'};
#'        otherwise set to \code{NULL}. When used, \code{y_c} should be
#'        \code{NULL}.
#' @param ne_t A positive integer giving the number of patients in the
#'        treatment group of the external data set. Required when
#'        \code{design = 'external'} and external treatment data are
#'        available; otherwise set to \code{NULL}.
#' @param ne_c A positive integer giving the number of patients in the
#'        control group of the external data set. Required when
#'        \code{design = 'external'} and external control data are available;
#'        otherwise set to \code{NULL}.
#' @param ye_t A non-negative integer giving the number of responders in the
#'        treatment group of the external data set. Required when
#'        \code{design = 'external'}; otherwise set to \code{NULL}.
#' @param ye_c A non-negative integer giving the number of responders in the
#'        control group of the external data set. Required when
#'        \code{design = 'external'}; otherwise set to \code{NULL}.
#' @param alpha0e_t A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for the external treatment data. Required when external
#'        treatment data are used; otherwise set to \code{NULL}.
#' @param alpha0e_c A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for the external control data. Required when external
#'        control data are used; otherwise set to \code{NULL}.
#' @param lower.tail A logical scalar; if \code{TRUE} (default), the function
#'        returns \eqn{P(\mathrm{effect} \le \theta_0)}, otherwise
#'        \eqn{P(\mathrm{effect} > \theta_0)}.
#'
#' @return A numeric scalar or vector in \code{[0, 1]}. When \code{y_t} and
#'         \code{y_c} are vectors of length \eqn{n}, a vector of length
#'         \eqn{n} is returned.
#'
#' @details
#' Posterior shape parameters are computed from the beta-binomial conjugate
#' model:
#' \itemize{
#'   \item Prior: \eqn{\pi_j \sim \mathrm{Beta}(a_j, b_j)}.
#'   \item Posterior:
#'         \eqn{\pi_j \mid y_j \sim \mathrm{Beta}(a_j + y_j, b_j + n_j - y_j)}.
#' }
#'
#' For \code{design = 'external'}, external data are incorporated through
#' the power prior, which inflates the prior by the weighted external
#' sufficient statistics:
#' \deqn{(\alpha_j^*, \beta_j^*) = (a_j + ae_j \cdot ye_j, b_j + ae_j \cdot
#'       (ne_j - ye_j)).}
#'
#' For \code{design = 'uncontrolled'}, the hypothetical control value
#' \code{z} is used in place of \code{y_c}, and \code{n_c} acts as the
#' hypothetical control sample size.
#'
#' The final probability is obtained by calling \code{\link{pbetadiff}} for
#' \code{prob = 'posterior'} or \code{\link{pbetabinomdiff}} for
#' \code{prob = 'predictive'}, applied element-wise via \code{mapply}.
#'
#' @examples
#' # Example 1: Controlled design - posterior probability
#' pbayespostpred1bin(
#'   prob = 'posterior', design = 'controlled', theta0 = 0.15,
#'   n_t = 12, n_c = 15, y_t = 7, y_c = 5,
#'   a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
#'   m_t = NULL, m_c = NULL, z = NULL,
#'   ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL,
#'   alpha0e_t = NULL, alpha0e_c = NULL, lower.tail = FALSE
#' )
#'
#' # Example 2: Uncontrolled design - posterior probability
#' pbayespostpred1bin(
#'   prob = 'posterior', design = 'uncontrolled', theta0 = 0.20,
#'   n_t = 20, n_c = 20, y_t = 12, y_c = NULL,
#'   a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
#'   m_t = NULL, m_c = NULL, z = 3,
#'   ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL,
#'   alpha0e_t = NULL, alpha0e_c = NULL, lower.tail = FALSE
#' )
#'
#' # Example 3: External design - posterior probability
#' pbayespostpred1bin(
#'   prob = 'posterior', design = 'external', theta0 = 0.15,
#'   n_t = 12, n_c = 15, y_t = 7, y_c = 9,
#'   a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
#'   m_t = NULL, m_c = NULL, z = NULL,
#'   ne_t = 12, ne_c = 12, ye_t = 6, ye_c = 6, alpha0e_t = 0.5, alpha0e_c = 0.5,
#'   lower.tail = FALSE
#' )
#'
#' # Example 4: Controlled design - posterior predictive probability
#' pbayespostpred1bin(
#'   prob = 'predictive', design = 'controlled', theta0 = 0.10,
#'   n_t = 12, n_c = 15, y_t = 7, y_c = 5,
#'   a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
#'   m_t = 30, m_c = 30, z = NULL,
#'   ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL,
#'   alpha0e_t = NULL, alpha0e_c = NULL, lower.tail = FALSE
#' )
#'
#' # Example 5: Uncontrolled design - posterior predictive probability
#' pbayespostpred1bin(
#'   prob = 'predictive', design = 'uncontrolled', theta0 = 0.20,
#'   n_t = 20, n_c = 20, y_t = 12, y_c = NULL,
#'   a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
#'   m_t = 30, m_c = 30, z = 3,
#'   ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL,
#'   alpha0e_t = NULL, alpha0e_c = NULL, lower.tail = FALSE
#' )
#'
#' # Example 6: External design - posterior predictive probability
#' pbayespostpred1bin(
#'   prob = 'predictive', design = 'external', theta0 = 0.15,
#'   n_t = 12, n_c = 15, y_t = 7, y_c = 9,
#'   a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
#'   m_t = 30, m_c = 30, z = NULL,
#'   ne_t = 12, ne_c = 12, ye_t = 6, ye_c = 6, alpha0e_t = 0.5, alpha0e_c = 0.5,
#'   lower.tail = FALSE
#' )
#'
#' @importFrom stats integrate
#' @export
pbayespostpred1bin <- function(prob = 'posterior', design = 'controlled', theta0,
                               n_t, n_c, y_t, y_c = NULL, a_t, a_c, b_t, b_c,
                               m_t = NULL, m_c = NULL, z = NULL,
                               ne_t = NULL, ne_c = NULL,
                               ye_t = NULL, ye_c = NULL,
                               alpha0e_t = NULL, alpha0e_c = NULL,
                               lower.tail = TRUE) {

  # --- Input validation ---
  if (!is.character(prob) || length(prob) != 1L ||
      !prob %in% c('posterior', 'predictive')) {
    stop("'prob' must be either 'posterior' or 'predictive'")
  }

  if (!is.character(design) || length(design) != 1L ||
      !design %in% c('controlled', 'uncontrolled', 'external')) {
    stop("'design' must be 'controlled', 'uncontrolled', or 'external'")
  }

  if (!is.numeric(theta0) || length(theta0) != 1L || is.na(theta0) ||
      theta0 <= -1 || theta0 >= 1) {
    stop("'theta0' must be a single numeric value in (-1, 1)")
  }

  if (!is.numeric(n_t) || length(n_t) != 1L || is.na(n_t) ||
      n_t != floor(n_t) || n_t < 1L) {
    stop("'n_t' must be a single positive integer")
  }

  if (!is.numeric(n_c) || length(n_c) != 1L || is.na(n_c) ||
      n_c != floor(n_c) || n_c < 1L) {
    stop("'n_c' must be a single positive integer")
  }

  for (nm in c("a_t", "a_c", "b_t", "b_c")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) || val <= 0) {
      stop(paste0("'", nm, "' must be a single positive numeric value"))
    }
  }

  if (!is.logical(lower.tail) || length(lower.tail) != 1L || is.na(lower.tail)) {
    stop("'lower.tail' must be a single logical value (TRUE or FALSE)")
  }

  # Validate prob-specific parameters
  if (prob == 'predictive') {
    if (is.null(m_t) || is.null(m_c)) {
      stop("'m_t' and 'm_c' must be non-NULL when prob = 'predictive'")
    }
    for (nm in c("m_t", "m_c")) {
      val <- get(nm)
      if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
          val != floor(val) || val < 1L) {
        stop(paste0("'", nm, "' must be a single positive integer"))
      }
    }
  }

  # Validate design-specific parameters
  if (design == 'uncontrolled') {
    if (is.null(z)) {
      stop("'z' must be non-NULL when design = 'uncontrolled'")
    }
    if (!is.numeric(z) || length(z) != 1L || is.na(z) ||
        z != floor(z) || z < 0 || z > n_c) {
      stop("'z' must be a single non-negative integer not exceeding n_c")
    }
  }

  if (design == 'external') {
    for (nm in c("ne_t", "ne_c", "ye_t", "ye_c", "alpha0e_t", "alpha0e_c")) {
      val <- get(nm)
      if (is.null(val)) {
        stop(paste0("'", nm, "' must be non-NULL when design = 'external'"))
      }
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
      ne_val <- get(sub("ye", "ne", nm))
      if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
          val != floor(val) || val < 0 || val > ne_val) {
        stop(paste0("'", nm, "' must be a single non-negative integer not exceeding the corresponding ne"))
      }
    }
    for (nm in c("alpha0e_t", "alpha0e_c")) {
      val <- get(nm)
      if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
          val <= 0 || val > 1) {
        stop(paste0("'", nm, "' must be a single numeric value in (0, 1]"))
      }
    }
  }

  # --- Coerce y_t to numeric vector ---
  y_t <- as.numeric(y_t)

  if (any(is.na(y_t)) || any(y_t < 0) || any(y_t > n_t) ||
      any(y_t != floor(y_t))) {
    stop("All elements of 'y_t' must be non-negative integers not exceeding n_t")
  }

  # --- Resolve y_c: fixed at z for uncontrolled design ---
  if (is.null(y_c)) {
    if (design == 'uncontrolled') {
      y_c <- rep(z, length(y_t))
    } else {
      stop("'y_c' must be provided for controlled and external designs")
    }
  }
  y_c <- as.numeric(y_c)

  if (any(is.na(y_c)) || any(y_c < 0) || any(y_c > n_c) ||
      any(y_c != floor(y_c))) {
    stop("All elements of 'y_c' must be non-negative integers not exceeding n_c")
  }

  if (length(y_t) != length(y_c)) {
    stop("'y_t' and 'y_c' must have the same length")
  }

  # --- Compute posterior shape parameters (vectorised) ---
  if (design == 'external') {
    # Power prior: augment the prior with weighted external sufficient statistics
    s_t1 <- y_t + a_t + alpha0e_t * ye_t
    s_t2 <- n_t - y_t + b_t + alpha0e_t * (ne_t - ye_t)
    s_c1 <- y_c + a_c + alpha0e_c * ye_c
    s_c2 <- n_c - y_c + b_c + alpha0e_c * (ne_c - ye_c)
  } else {
    # Controlled or uncontrolled (y_c already set to z for uncontrolled)
    s_t1 <- y_t + a_t
    s_t2 <- n_t - y_t + b_t
    s_c1 <- y_c + a_c
    s_c2 <- n_c - y_c + b_c
  }

  # --- Compute probabilities element-wise via mapply ---
  if (prob == 'posterior') {
    # P(pi_t - pi_c > theta0 | data) using Beta posterior distributions
    results <- mapply(
      pbetadiff,
      alpha_t = s_t1, alpha_c = s_c1, beta_t = s_t2, beta_c = s_c2,
      MoreArgs = list(q = theta0, lower.tail = lower.tail)
    )
  } else {
    # P(future trial success | data) using Beta-Binomial predictive distributions
    results <- mapply(
      pbetabinomdiff,
      alpha_t = s_t1, alpha_c = s_c1, beta_t = s_t2, beta_c = s_c2,
      MoreArgs = list(q = theta0, m_t = m_t, m_c = m_c, lower.tail = lower.tail)
    )
  }

  return(results)
}
