#' Bayesian Posterior or Posterior Predictive Probability for a Single Binary Endpoint
#'
#' Computes the Bayesian posterior probability or posterior predictive probability
#' for binary-outcome clinical trials under a beta-binomial conjugate model.
#' The function supports controlled, uncontrolled, and external-control designs,
#' with optional incorporation of external data through power priors.
#' Vector inputs for \code{y1} and \code{y2} are supported for efficient batch
#' processing (e.g., across all possible trial outcomes in
#' \code{\link{pGNGsinglebinary}}).
#'
#' @param prob A character string specifying the probability type.
#'        Must be \code{'posterior'} or \code{'predictive'}.
#' @param design A character string specifying the trial design.
#'        Must be \code{'controlled'}, \code{'uncontrolled'}, or \code{'external'}.
#' @param theta0 A numeric scalar in \code{(-1, 1)} giving the pre-specified
#'        threshold for the treatment effect (difference in response rates).
#' @param n1 A positive integer giving the number of patients in group 1
#'        (treatment) in the proof-of-concept (PoC) trial.
#' @param n2 A positive integer giving the number of patients in group 2 in the
#'        PoC trial. For \code{design = 'uncontrolled'}, this is the hypothetical
#'        control sample size (required for consistency with other designs).
#' @param y1 A non-negative integer or integer vector giving the observed number
#'        of responders in group 1 (must satisfy \code{0 <= y1 <= n1}).
#' @param y2 A non-negative integer or integer vector giving the number of
#'        responders in group 2 (must satisfy \code{0 <= y2 <= n2}).
#'        Set to \code{NULL} for \code{design = 'uncontrolled'} and use \code{z}
#'        instead. When provided as a vector, must have the same length as
#'        \code{y1}.
#' @param a1 A positive numeric scalar giving the first shape parameter (alpha)
#'        of the prior Beta distribution for group 1.
#' @param a2 A positive numeric scalar giving the first shape parameter (alpha)
#'        of the prior Beta distribution for group 2.
#' @param b1 A positive numeric scalar giving the second shape parameter (beta)
#'        of the prior Beta distribution for group 1.
#' @param b2 A positive numeric scalar giving the second shape parameter (beta)
#'        of the prior Beta distribution for group 2.
#' @param m1 A positive integer giving the number of patients in group 1 for
#'        the future trial. Required when \code{prob = 'predictive'};
#'        otherwise set to \code{NULL}.
#' @param m2 A positive integer giving the number of patients in group 2 for
#'        the future trial. Required when \code{prob = 'predictive'};
#'        otherwise set to \code{NULL}.
#' @param z A non-negative integer giving the hypothetical number of responders
#'        in the control group. Required when \code{design = 'uncontrolled'};
#'        otherwise set to \code{NULL}. When used, \code{y2} should be
#'        \code{NULL}.
#' @param ne1 A positive integer giving the number of patients in group 1 of
#'        the external data set. Required when \code{design = 'external'} and
#'        external treatment data are available; otherwise set to \code{NULL}.
#' @param ne2 A positive integer giving the number of patients in group 2 of
#'        the external data set. Required when \code{design = 'external'} and
#'        external control data are available; otherwise set to \code{NULL}.
#' @param ye1 A non-negative integer giving the number of responders in group 1
#'        of the external data set. Required when \code{design = 'external'};
#'        otherwise set to \code{NULL}.
#' @param ye2 A non-negative integer giving the number of responders in group 2
#'        of the external data set. Required when \code{design = 'external'};
#'        otherwise set to \code{NULL}.
#' @param ae1 A numeric scalar in \code{(0, 1]} giving the power prior weight
#'        for group 1. Controls the degree of borrowing: 1 = full borrowing.
#'        Required when \code{design = 'external'}; otherwise set to \code{NULL}.
#' @param ae2 A numeric scalar in \code{(0, 1]} giving the power prior weight
#'        for group 2. Required when \code{design = 'external'};
#'        otherwise set to \code{NULL}.
#' @param lower.tail A logical scalar; if \code{TRUE} (default), the function
#'        returns \eqn{P(\mathrm{effect} \le \theta_0)}, otherwise
#'        \eqn{P(\mathrm{effect} > \theta_0)}.
#'
#' @return A numeric scalar or vector in \code{[0, 1]}.  When \code{y1} and
#'         \code{y2} are vectors of length \eqn{n}, a vector of length \eqn{n}
#'         is returned.
#'
#' @details
#' Posterior shape parameters are computed from the beta-binomial conjugate model:
#' \itemize{
#'   \item Prior: \eqn{\pi_j \sim \mathrm{Beta}(a_j,\, b_j)}.
#'   \item Posterior: \eqn{\pi_j \mid y_j \sim
#'         \mathrm{Beta}(a_j + y_j,\; b_j + n_j - y_j)}.
#' }
#'
#' For \code{design = 'external'}, external data are incorporated through the
#' power prior, which inflates the prior by the weighted external sufficient
#' statistics:
#' \deqn{(\alpha_j^*, \beta_j^*) =
#'   (a_j + ae_j \cdot ye_j,\; b_j + ae_j \cdot (ne_j - ye_j)).}
#'
#' For \code{design = 'uncontrolled'}, the hypothetical control value \code{z}
#' is used in place of \code{y2}, and \code{n2} acts as the hypothetical control
#' sample size.
#'
#' The final probability is obtained by calling \code{\link{p2betadiff}} for
#' \code{prob = 'posterior'} or \code{\link{p2betabinomdiff}} for
#' \code{prob = 'predictive'}, applied element-wise via \code{mapply}.
#'
#' @examples
#' # Example 1: Controlled design - posterior probability (scalar)
#' pPPsinglebinary(
#'   prob = 'posterior', design = 'controlled', theta0 = 0.15,
#'   n1 = 12, n2 = 15, y1 = 7, y2 = 5,
#'   a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
#'   m1 = NULL, m2 = NULL, z = NULL,
#'   ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL,
#'   lower.tail = FALSE
#' )
#'
#' # Example 2: Controlled design - vectorised over outcomes
#' pPPsinglebinary(
#'   prob = 'posterior', design = 'controlled', theta0 = 0.15,
#'   n1 = 12, n2 = 15, y1 = c(5, 6, 7, 8), y2 = c(3, 4, 5, 6),
#'   a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
#'   m1 = NULL, m2 = NULL, z = NULL,
#'   ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL,
#'   lower.tail = FALSE
#' )
#'
#' # Example 3: Uncontrolled design - hypothetical control via z
#' pPPsinglebinary(
#'   prob = 'posterior', design = 'uncontrolled', theta0 = 0.20,
#'   n1 = 20, n2 = 20, y1 = 12, y2 = NULL,
#'   a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
#'   m1 = NULL, m2 = NULL, z = 3,
#'   ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL,
#'   lower.tail = FALSE
#' )
#'
#' # Example 4: Controlled design - posterior predictive probability
#' pPPsinglebinary(
#'   prob = 'predictive', design = 'controlled', theta0 = 0.10,
#'   n1 = 12, n2 = 15, y1 = 7, y2 = 5,
#'   a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
#'   m1 = 30, m2 = 30, z = NULL,
#'   ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL,
#'   lower.tail = FALSE
#' )
#'
#' # Example 5: External design - power prior with 50 percent borrowing
#' pPPsinglebinary(
#'   prob = 'posterior', design = 'external', theta0 = 0.15,
#'   n1 = 12, n2 = 15, y1 = 7, y2 = 9,
#'   a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
#'   m1 = NULL, m2 = NULL, z = NULL,
#'   ne1 = 12, ne2 = 12, ye1 = 6, ye2 = 6, ae1 = 0.5, ae2 = 0.5,
#'   lower.tail = FALSE
#' )
#'
#' @importFrom stats integrate
#' @export
pPPsinglebinary <- function(prob = 'posterior', design = 'controlled', theta0,
                            n1, n2, y1, y2, a1, a2, b1, b2,
                            m1, m2, z = NULL,
                            ne1, ne2, ye1, ye2, ae1, ae2,
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

  if (!is.numeric(n1) || length(n1) != 1L || is.na(n1) ||
      n1 != floor(n1) || n1 < 1L) {
    stop("'n1' must be a single positive integer")
  }

  if (!is.numeric(n2) || length(n2) != 1L || is.na(n2) ||
      n2 != floor(n2) || n2 < 1L) {
    stop("'n2' must be a single positive integer")
  }

  for (nm in c("a1", "a2", "b1", "b2")) {
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

  # --- Coerce y1 to a numeric vector ---
  y1 <- as.numeric(y1)

  if (any(is.na(y1)) || any(y1 < 0) || any(y1 > n1) ||
      any(y1 != floor(y1))) {
    stop("All elements of 'y1' must be non-negative integers not exceeding n1")
  }

  # --- Resolve y2: fixed at z for uncontrolled design ---
  if (is.null(y2)) {
    if (design == 'uncontrolled') {
      y2 <- rep(z, length(y1))
    } else {
      stop("'y2' must be provided for controlled and external designs")
    }
  }
  y2 <- as.numeric(y2)

  if (any(is.na(y2)) || any(y2 < 0) || any(y2 > n2) ||
      any(y2 != floor(y2))) {
    stop("All elements of 'y2' must be non-negative integers not exceeding n2")
  }

  if (length(y1) != length(y2)) {
    stop("'y1' and 'y2' must have the same length")
  }

  # --- Compute posterior shape parameters (vectorised) ---
  if (design == 'external') {
    # Power prior: augment the prior with weighted external sufficient statistics
    s11 <- y1 + a1 + ae1 * ye1
    s21 <- n1 - y1 + b1 + ae1 * (ne1 - ye1)
    s12 <- y2 + a2 + ae2 * ye2
    s22 <- n2 - y2 + b2 + ae2 * (ne2 - ye2)
  } else {
    # Controlled or uncontrolled (y2 already set to z for uncontrolled)
    s11 <- y1 + a1
    s21 <- n1 - y1 + b1
    s12 <- y2 + a2
    s22 <- n2 - y2 + b2
  }

  # --- Compute probabilities element-wise via mapply ---
  if (prob == 'posterior') {
    # P(pi1 - pi2 > theta0 | data) using Beta posterior distributions
    results <- mapply(
      p2betadiff,
      alpha1 = s11, alpha2 = s12, beta1 = s21, beta2 = s22,
      MoreArgs = list(q = theta0, lower.tail = lower.tail)
    )
  } else {
    # P(future trial success | data) using Beta-Binomial predictive distributions
    results <- mapply(
      p2betabinomdiff,
      alpha1 = s11, alpha2 = s12, beta1 = s21, beta2 = s22,
      MoreArgs = list(q = theta0, m1 = m1, m2 = m2, lower.tail = lower.tail)
    )
  }

  return(results)
}
