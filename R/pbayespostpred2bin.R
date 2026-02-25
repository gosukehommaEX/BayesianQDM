#' Bayesian Posterior or Posterior Predictive Probability for a Clinical Trial
#' with Two Binary Endpoints
#'
#' Computes the Bayesian posterior probability or posterior predictive
#' probability for clinical trials with two binary endpoints under a
#' Dirichlet-multinomial conjugate model.  The function returns probabilities
#' for nine decision regions (posterior) or four decision regions (predictive)
#' defined by target values (TV) and minimum acceptable values (MAV) for both
#' endpoints.  Three study designs are supported: controlled, uncontrolled
#' (hypothetical control), and external-control (power prior).
#'
#' @param prob A character string specifying the probability type.
#'        Must be \code{'posterior'} or \code{'predictive'}.
#' @param design A character string specifying the trial design.
#'        Must be \code{'controlled'}, \code{'uncontrolled'}, or
#'        \code{'external'}.
#' @param theta.TV1 A numeric scalar giving the target value (TV) threshold
#'        for Endpoint 1.  Required when \code{prob = 'posterior'}; must
#'        satisfy \code{theta.TV1 > theta.MAV1}.  Set to \code{NULL} when
#'        \code{prob = 'predictive'}.
#' @param theta.MAV1 A numeric scalar giving the minimum acceptable value
#'        (MAV) threshold for Endpoint 1.  Required when
#'        \code{prob = 'posterior'}; must satisfy
#'        \code{theta.TV1 > theta.MAV1}.  Set to \code{NULL} when
#'        \code{prob = 'predictive'}.
#' @param theta.TV2 A numeric scalar giving the target value (TV) threshold
#'        for Endpoint 2.  Required when \code{prob = 'posterior'}; must
#'        satisfy \code{theta.TV2 > theta.MAV2}.  Set to \code{NULL} when
#'        \code{prob = 'predictive'}.
#' @param theta.MAV2 A numeric scalar giving the minimum acceptable value
#'        (MAV) threshold for Endpoint 2.  Required when
#'        \code{prob = 'posterior'}; must satisfy
#'        \code{theta.TV2 > theta.MAV2}.  Set to \code{NULL} when
#'        \code{prob = 'predictive'}.
#' @param theta.NULL1 A numeric scalar giving the null hypothesis threshold
#'        for Endpoint 1.  Required when \code{prob = 'predictive'};
#'        set to \code{NULL} when \code{prob = 'posterior'}.
#' @param theta.NULL2 A numeric scalar giving the null hypothesis threshold
#'        for Endpoint 2.  Required when \code{prob = 'predictive'};
#'        set to \code{NULL} when \code{prob = 'posterior'}.
#' @param x1_00 A non-negative integer giving the count of (0,0) responses
#'        in group 1 (Endpoint 1 = 0, Endpoint 2 = 0).
#' @param x1_01 A non-negative integer giving the count of (0,1) responses
#'        in group 1.
#' @param x1_10 A non-negative integer giving the count of (1,0) responses
#'        in group 1.
#' @param x1_11 A non-negative integer giving the count of (1,1) responses
#'        in group 1.
#' @param x2_00 A non-negative integer giving the count of (0,0) responses
#'        in group 2.  Not used when \code{design = 'uncontrolled'}; set to
#'        \code{NULL} in that case.
#' @param x2_01 A non-negative integer giving the count of (0,1) responses
#'        in group 2.  Not used when \code{design = 'uncontrolled'}; set to
#'        \code{NULL} in that case.
#' @param x2_10 A non-negative integer giving the count of (1,0) responses
#'        in group 2.  Not used when \code{design = 'uncontrolled'}; set to
#'        \code{NULL} in that case.
#' @param x2_11 A non-negative integer giving the count of (1,1) responses
#'        in group 2.  Not used when \code{design = 'uncontrolled'}; set to
#'        \code{NULL} in that case.
#' @param a1_00 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the (0,0) response pattern in group 1.
#' @param a1_01 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the (0,1) response pattern in group 1.
#' @param a1_10 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the (1,0) response pattern in group 1.
#' @param a1_11 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the (1,1) response pattern in group 1.
#' @param a2_00 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the (0,0) response pattern in group 2.  For
#'        \code{design = 'uncontrolled'}, this serves as a hyperparameter
#'        of the hypothetical control distribution (see Details).
#' @param a2_01 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the (0,1) response pattern in group 2.  For
#'        \code{design = 'uncontrolled'}, serves as a hyperparameter of
#'        the hypothetical control distribution.
#' @param a2_10 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the (1,0) response pattern in group 2.  For
#'        \code{design = 'uncontrolled'}, serves as a hyperparameter of
#'        the hypothetical control distribution.
#' @param a2_11 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the (1,1) response pattern in group 2.  For
#'        \code{design = 'uncontrolled'}, serves as a hyperparameter of
#'        the hypothetical control distribution.
#' @param m1 A positive integer giving the future sample size for group 1.
#'        Required when \code{prob = 'predictive'}; otherwise set to
#'        \code{NULL}.
#' @param m2 A positive integer giving the future sample size for group 2.
#'        Required when \code{prob = 'predictive'}; otherwise set to
#'        \code{NULL}.
#' @param z00 A non-negative integer giving the hypothetical control count
#'        for pattern (0,0).  Required when \code{design = 'uncontrolled'};
#'        otherwise set to \code{NULL}.
#' @param z01 A non-negative integer giving the hypothetical control count
#'        for pattern (0,1).  Required when \code{design = 'uncontrolled'};
#'        otherwise set to \code{NULL}.
#' @param z10 A non-negative integer giving the hypothetical control count
#'        for pattern (1,0).  Required when \code{design = 'uncontrolled'};
#'        otherwise set to \code{NULL}.
#' @param z11 A non-negative integer giving the hypothetical control count
#'        for pattern (1,1).  Required when \code{design = 'uncontrolled'};
#'        otherwise set to \code{NULL}.
#' @param xe1_00 A non-negative integer giving the external group 1 count
#'        for pattern (0,0).  Required when \code{design = 'external'} and
#'        external treatment data are used; otherwise \code{NULL}.
#' @param xe1_01 A non-negative integer giving the external group 1 count
#'        for pattern (0,1).  Required for external treatment data;
#'        otherwise \code{NULL}.
#' @param xe1_10 A non-negative integer giving the external group 1 count
#'        for pattern (1,0).  Required for external treatment data;
#'        otherwise \code{NULL}.
#' @param xe1_11 A non-negative integer giving the external group 1 count
#'        for pattern (1,1).  Required for external treatment data;
#'        otherwise \code{NULL}.
#' @param xe2_00 A non-negative integer giving the external group 2 count
#'        for pattern (0,0).  Required when \code{design = 'external'} and
#'        external control data are used; otherwise \code{NULL}.
#' @param xe2_01 A non-negative integer giving the external group 2 count
#'        for pattern (0,1).  Required for external control data; otherwise
#'        \code{NULL}.
#' @param xe2_10 A non-negative integer giving the external group 2 count
#'        for pattern (1,0).  Required for external control data; otherwise
#'        \code{NULL}.
#' @param xe2_11 A non-negative integer giving the external group 2 count
#'        for pattern (1,1).  Required for external control data; otherwise
#'        \code{NULL}.
#' @param ae1 A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for group 1.  Required when external treatment data are
#'        used; otherwise \code{NULL}.
#' @param ae2 A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for group 2.  Required when external control data are
#'        used; otherwise \code{NULL}.
#' @param nMC A positive integer giving the number of Monte Carlo draws
#'        used to estimate region probabilities.  Default is \code{10000}.
#'        Larger values reduce Monte Carlo error at the cost of
#'        computation time.
#'
#' @return A named numeric vector of region probabilities.
#'         For \code{prob = 'posterior'}: length 9, named \code{R1}, ...,
#'         \code{R9}, corresponding to regions defined by TV and MAV
#'         thresholds for both endpoints (row-major order: Endpoint 1
#'         varies slowest).
#'         For \code{prob = 'predictive'}: length 4, named \code{R1}, ...,
#'         \code{R4}.
#'         All elements are non-negative and sum to 1.
#'
#' @details
#' \strong{Model.}
#' The four response categories are ordered as (0,0), (0,1), (1,0), (1,1).
#' For each arm \eqn{k \in \{1, 2\}}, the observed count vector follows a
#' multinomial distribution, and a conjugate Dirichlet prior is placed on
#' the cell probability vector \eqn{p_k}:
#' \deqn{p_k \sim \mathrm{Dir}(\alpha_{k,00},\, \alpha_{k,01},\,
#'       \alpha_{k,10},\, \alpha_{k,11}).}
#' The posterior is
#' \deqn{p_k \mid x_k \sim
#'       \mathrm{Dir}(\alpha_{k,00} + x_{k,00},\; \alpha_{k,01} + x_{k,01},\;
#'       \alpha_{k,10} + x_{k,10},\; \alpha_{k,11} + x_{k,11}).}
#' Marginal response rates are
#' \eqn{\pi_{k1} = p_{k,10} + p_{k,11}} (Endpoint 1) and
#' \eqn{\pi_{k2} = p_{k,01} + p_{k,11}} (Endpoint 2).
#' Treatment effects are
#' \eqn{\theta_1 = \pi_{11} - \pi_{21}} and
#' \eqn{\theta_2 = \pi_{12} - \pi_{22}}.
#'
#' \strong{Uncontrolled design.}
#' When \code{design = 'uncontrolled'}, no concurrent control data are
#' available.  The hypothetical control distribution is specified as a
#' Dirichlet distribution directly via the prior hyperparameters
#' \code{a2_00}, \code{a2_01}, \code{a2_10}, \code{a2_11}, augmented by
#' hypothetical control counts \code{z00}, \code{z01}, \code{z10},
#' \code{z11}:
#' \deqn{p_2 \sim \mathrm{Dir}(\alpha_{2,00} + z_{00},\; \alpha_{2,01} +
#'       z_{01},\; \alpha_{2,10} + z_{10},\; \alpha_{2,11} + z_{11}).}
#' Monte Carlo samples from this distribution replace data-derived posterior
#' samples for group 2 in all subsequent computations.
#'
#' \strong{External design.}
#' When \code{design = 'external'}, the power prior augments the Dirichlet
#' prior parameters:
#' \deqn{\alpha_{k,lm}^{*} = \alpha_{k,lm} + a_{0k} \cdot x_{e k,lm}.}
#'
#' \strong{Posterior probability regions (prob = 'posterior').}
#' The nine regions follow Table 1 of the reference (row-major order,
#' Endpoint 1 varies slowest):
#' \itemize{
#'   \item R1: \eqn{\theta_1 > TV_1} AND \eqn{\theta_2 > TV_2}
#'   \item R2: \eqn{\theta_1 > TV_1} AND \eqn{TV_2 \ge \theta_2 > MAV_2}
#'   \item R3: \eqn{\theta_1 > TV_1} AND \eqn{\theta_2 \le MAV_2}
#'   \item R4: \eqn{TV_1 \ge \theta_1 > MAV_1} AND \eqn{\theta_2 > TV_2}
#'   \item R5: \eqn{TV_1 \ge \theta_1 > MAV_1} AND
#'             \eqn{TV_2 \ge \theta_2 > MAV_2}
#'   \item R6: \eqn{TV_1 \ge \theta_1 > MAV_1} AND \eqn{\theta_2 \le MAV_2}
#'   \item R7: \eqn{\theta_1 \le MAV_1} AND \eqn{\theta_2 > TV_2}
#'   \item R8: \eqn{\theta_1 \le MAV_1} AND
#'             \eqn{TV_2 \ge \theta_2 > MAV_2}
#'   \item R9: \eqn{\theta_1 \le MAV_1} AND \eqn{\theta_2 \le MAV_2}
#' }
#'
#' \strong{Predictive probability regions (prob = 'predictive').}
#' The four regions follow Table 2 of the reference:
#' \itemize{
#'   \item R1: \eqn{\tilde\theta_1 > NULL_1} AND
#'             \eqn{\tilde\theta_2 > NULL_2}
#'   \item R2: \eqn{\tilde\theta_1 > NULL_1} AND
#'             \eqn{\tilde\theta_2 \le NULL_2}
#'   \item R3: \eqn{\tilde\theta_1 \le NULL_1} AND
#'             \eqn{\tilde\theta_2 > NULL_2}
#'   \item R4: \eqn{\tilde\theta_1 \le NULL_1} AND
#'             \eqn{\tilde\theta_2 \le NULL_2}
#' }
#'
#' @examples
#' # Example 1: Posterior probability, controlled design
#' pbayespostpred2bin(
#'   prob = 'posterior', design = 'controlled',
#'   theta.TV1 = 0.20, theta.MAV1 = 0.10,
#'   theta.TV2 = 0.20, theta.MAV2 = 0.10,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   x1_00 = 5, x1_01 = 3, x1_10 = 4, x1_11 = 8,
#'   x2_00 = 8, x2_01 = 4, x2_10 = 5, x2_11 = 3,
#'   a1_00 = 0.5, a1_01 = 0.5, a1_10 = 0.5, a1_11 = 0.5,
#'   a2_00 = 0.5, a2_01 = 0.5, a2_10 = 0.5, a2_11 = 0.5,
#'   m1 = NULL, m2 = NULL,
#'   z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
#'   xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
#'   xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
#'   ae1 = NULL, ae2 = NULL,
#'   nMC = 10000
#' )
#'
#' # Example 2: Posterior predictive probability, controlled design
#' pbayespostpred2bin(
#'   prob = 'predictive', design = 'controlled',
#'   theta.TV1 = NULL, theta.MAV1 = NULL,
#'   theta.TV2 = NULL, theta.MAV2 = NULL,
#'   theta.NULL1 = 0.15, theta.NULL2 = 0.15,
#'   x1_00 = 3, x1_01 = 2, x1_10 = 3, x1_11 = 7,
#'   x2_00 = 6, x2_01 = 3, x2_10 = 4, x2_11 = 2,
#'   a1_00 = 1, a1_01 = 1, a1_10 = 1, a1_11 = 1,
#'   a2_00 = 1, a2_01 = 1, a2_10 = 1, a2_11 = 1,
#'   m1 = 50, m2 = 50,
#'   z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
#'   xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
#'   xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
#'   ae1 = NULL, ae2 = NULL,
#'   nMC = 10000
#' )
#'
#' # Example 3: Posterior probability, external control design
#' pbayespostpred2bin(
#'   prob = 'posterior', design = 'external',
#'   theta.TV1 = 0.20, theta.MAV1 = 0.10,
#'   theta.TV2 = 0.20, theta.MAV2 = 0.10,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   x1_00 = 5, x1_01 = 3, x1_10 = 4, x1_11 = 8,
#'   x2_00 = 8, x2_01 = 4, x2_10 = 5, x2_11 = 3,
#'   a1_00 = 0.5, a1_01 = 0.5, a1_10 = 0.5, a1_11 = 0.5,
#'   a2_00 = 0.5, a2_01 = 0.5, a2_10 = 0.5, a2_11 = 0.5,
#'   m1 = NULL, m2 = NULL,
#'   z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
#'   xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
#'   xe2_00 = 3L, xe2_01 = 2L, xe2_10 = 3L, xe2_11 = 2L,
#'   ae1 = NULL, ae2 = 0.5,
#'   nMC = 10000
#' )
#'
#' @importFrom stats rbinom
#' @export
pbayespostpred2bin <- function(prob    = 'posterior',
                               design  = 'controlled',
                               theta.TV1   = NULL, theta.MAV1  = NULL,
                               theta.TV2   = NULL, theta.MAV2  = NULL,
                               theta.NULL1 = NULL, theta.NULL2 = NULL,
                               x1_00, x1_01, x1_10, x1_11,
                               x2_00  = NULL, x2_01  = NULL,
                               x2_10  = NULL, x2_11  = NULL,
                               a1_00, a1_01, a1_10, a1_11,
                               a2_00, a2_01, a2_10, a2_11,
                               m1 = NULL, m2 = NULL,
                               z00 = NULL, z01 = NULL,
                               z10 = NULL, z11 = NULL,
                               xe1_00 = NULL, xe1_01 = NULL,
                               xe1_10 = NULL, xe1_11 = NULL,
                               xe2_00 = NULL, xe2_01 = NULL,
                               xe2_10 = NULL, xe2_11 = NULL,
                               ae1 = NULL, ae2 = NULL,
                               nMC = 10000L) {

  # ---------------------------------------------------------------------------
  # Section 1: Input validation
  # ---------------------------------------------------------------------------

  if (!is.character(prob) || length(prob) != 1L ||
      !prob %in% c('posterior', 'predictive')) {
    stop("'prob' must be either 'posterior' or 'predictive'")
  }

  if (!is.character(design) || length(design) != 1L ||
      !design %in% c('controlled', 'uncontrolled', 'external')) {
    stop("'design' must be 'controlled', 'uncontrolled', or 'external'")
  }

  if (!is.numeric(nMC) || length(nMC) != 1L || is.na(nMC) ||
      nMC != floor(nMC) || nMC < 1L) {
    stop("'nMC' must be a single positive integer")
  }
  nMC <- as.integer(nMC)

  # --- Threshold parameters ---
  if (prob == 'posterior') {
    for (nm in c("theta.TV1", "theta.MAV1", "theta.TV2", "theta.MAV2")) {
      val <- get(nm)
      if (is.null(val)) {
        stop(paste0("'", nm, "' must be non-NULL when prob = 'posterior'"))
      }
      if (!is.numeric(val) || length(val) != 1L || is.na(val)) {
        stop(paste0("'", nm, "' must be a single numeric value"))
      }
    }
    if (theta.TV1 <= theta.MAV1) {
      stop("'theta.TV1' must be strictly greater than 'theta.MAV1'")
    }
    if (theta.TV2 <= theta.MAV2) {
      stop("'theta.TV2' must be strictly greater than 'theta.MAV2'")
    }
  } else {
    for (nm in c("theta.NULL1", "theta.NULL2")) {
      val <- get(nm)
      if (is.null(val)) {
        stop(paste0("'", nm, "' must be non-NULL when prob = 'predictive'"))
      }
      if (!is.numeric(val) || length(val) != 1L || is.na(val)) {
        stop(paste0("'", nm, "' must be a single numeric value"))
      }
    }
    # Convert null thresholds: TV = MAV = NULL for predictive probability
    theta.TV1  <- theta.NULL1;  theta.MAV1 <- theta.NULL1
    theta.TV2  <- theta.NULL2;  theta.MAV2 <- theta.NULL2
  }

  # --- Group 1 observed counts (always required) ---
  for (nm in c("x1_00", "x1_01", "x1_10", "x1_11")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
        val != floor(val) || val < 0L) {
      stop(paste0("'", nm, "' must be a single non-negative integer"))
    }
  }

  # --- Group 2 observed counts (required unless design = 'uncontrolled') ---
  if (design != 'uncontrolled') {
    for (nm in c("x2_00", "x2_01", "x2_10", "x2_11")) {
      val <- get(nm)
      if (is.null(val)) {
        stop(paste0("'", nm, "' must be non-NULL when design != 'uncontrolled'"))
      }
      if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
          val != floor(val) || val < 0L) {
        stop(paste0("'", nm, "' must be a single non-negative integer"))
      }
    }
  } else {
    # Warn if group 2 data are accidentally supplied for uncontrolled design
    if (!is.null(x2_00) || !is.null(x2_01) ||
        !is.null(x2_10) || !is.null(x2_11)) {
      warning(paste0("design = 'uncontrolled': 'x2_00', 'x2_01', 'x2_10', ",
                     "'x2_11' are ignored. Use 'z00', 'z01', 'z10', 'z11' ",
                     "to specify the hypothetical control counts."))
    }
  }

  # --- Dirichlet prior parameters (always required) ---
  for (nm in c("a1_00", "a1_01", "a1_10", "a1_11",
               "a2_00", "a2_01", "a2_10", "a2_11")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) || val <= 0) {
      stop(paste0("'", nm, "' must be a single positive numeric value"))
    }
  }

  # --- Future sample sizes (required for predictive) ---
  if (prob == 'predictive') {
    if (is.null(m1) || is.null(m2)) {
      stop("'m1' and 'm2' must be non-NULL when prob = 'predictive'")
    }
    for (nm in c("m1", "m2")) {
      val <- get(nm)
      if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
          val != floor(val) || val < 1L) {
        stop(paste0("'", nm, "' must be a single positive integer"))
      }
    }
    m1 <- as.integer(m1)
    m2 <- as.integer(m2)
  }

  # --- Hypothetical control counts (required for uncontrolled) ---
  if (design == 'uncontrolled') {
    for (nm in c("z00", "z01", "z10", "z11")) {
      val <- get(nm)
      if (is.null(val)) {
        stop(paste0("'", nm, "' must be non-NULL when design = 'uncontrolled'"))
      }
      if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
          val != floor(val) || val < 0L) {
        stop(paste0("'", nm, "' must be a single non-negative integer"))
      }
    }
  }

  # --- External data parameters (required for external) ---
  if (design == 'external') {
    has_ext1 <- !is.null(xe1_00) && !is.null(xe1_01) &&
      !is.null(xe1_10) && !is.null(xe1_11) && !is.null(ae1)
    has_ext2 <- !is.null(xe2_00) && !is.null(xe2_01) &&
      !is.null(xe2_10) && !is.null(xe2_11) && !is.null(ae2)

    if (!has_ext1 && !has_ext2) {
      stop(paste0("For design = 'external', at least one complete set of ",
                  "external data (xe1_* + ae1 or xe2_* + ae2) must be provided"))
    }

    for (nm in c("ae1", "ae2")) {
      val <- get(nm)
      if (!is.null(val)) {
        if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
            val <= 0 || val > 1) {
          stop(paste0("'", nm, "' must be a single numeric value in (0, 1]"))
        }
      }
    }

    ext1_nms <- c("xe1_00", "xe1_01", "xe1_10", "xe1_11")
    ext2_nms <- c("xe2_00", "xe2_01", "xe2_10", "xe2_11")

    if (has_ext1) {
      for (nm in ext1_nms) {
        val <- get(nm)
        if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
            val != floor(val) || val < 0L) {
          stop(paste0("'", nm, "' must be a single non-negative integer"))
        }
      }
    }
    if (has_ext2) {
      for (nm in ext2_nms) {
        val <- get(nm)
        if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
            val != floor(val) || val < 0L) {
          stop(paste0("'", nm, "' must be a single non-negative integer"))
        }
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Section 2: Build posterior Dirichlet parameters
  # ---------------------------------------------------------------------------

  # Group 1: base prior + observed data + optional external data (power prior)
  xe1_weight <- if (!is.null(ae1) && design == 'external') ae1 else 0
  alpha1 <- c(
    a1_00 + x1_00 + xe1_weight * ifelse(!is.null(xe1_00), xe1_00, 0),
    a1_01 + x1_01 + xe1_weight * ifelse(!is.null(xe1_01), xe1_01, 0),
    a1_10 + x1_10 + xe1_weight * ifelse(!is.null(xe1_10), xe1_10, 0),
    a1_11 + x1_11 + xe1_weight * ifelse(!is.null(xe1_11), xe1_11, 0)
  )

  # Group 2: parameterisation depends on design
  if (design == 'uncontrolled') {
    # Hypothetical control distribution: prior augmented by pseudo-counts z*
    # (Appendix A1.1.6: P_hyp(p_c) = Dir(beta_c00, ..., beta_c11) where
    #  beta_c** = a2_** + z**)
    alpha2 <- c(
      a2_00 + z00,
      a2_01 + z01,
      a2_10 + z10,
      a2_11 + z11
    )
  } else {
    # Controlled or external: posterior from observed control data
    xe2_weight <- if (!is.null(ae2) && design == 'external') ae2 else 0
    alpha2 <- c(
      a2_00 + x2_00 + xe2_weight * ifelse(!is.null(xe2_00), xe2_00, 0),
      a2_01 + x2_01 + xe2_weight * ifelse(!is.null(xe2_01), xe2_01, 0),
      a2_10 + x2_10 + xe2_weight * ifelse(!is.null(xe2_10), xe2_10, 0),
      a2_11 + x2_11 + xe2_weight * ifelse(!is.null(xe2_11), xe2_11, 0)
    )
  }

  # ---------------------------------------------------------------------------
  # Section 3: Monte Carlo estimation
  # ---------------------------------------------------------------------------

  # Draw nMC samples from the posterior (or hypothetical) Dirichlet for each arm
  # Columns of p1/p2: (p_00, p_01, p_10, p_11)
  p1 <- rdirichlet(nMC, alpha1)
  p2 <- rdirichlet(nMC, alpha2)

  if (prob == 'posterior') {

    # Marginal response rates:
    #   pi_k1 = p_k10 + p_k11  (Endpoint 1 success)
    #   pi_k2 = p_k01 + p_k11  (Endpoint 2 success)
    theta1 <- (p1[, 3L] + p1[, 4L]) - (p2[, 3L] + p2[, 4L])
    theta2 <- (p1[, 2L] + p1[, 4L]) - (p2[, 2L] + p2[, 4L])

    # Classify each draw into one of 9 regions (row-major: Endpoint 1 slowest)
    # r1: 1 = theta1 > TV1, 2 = MAV1 < theta1 <= TV1, 3 = theta1 <= MAV1
    # r2: 1 = theta2 > TV2, 2 = MAV2 < theta2 <= TV2, 3 = theta2 <= MAV2
    r1 <- 3L -
      as.integer(theta1 > theta.MAV1) -
      as.integer(theta1 > theta.TV1)
    r2 <- 3L -
      as.integer(theta2 > theta.MAV2) -
      as.integer(theta2 > theta.TV2)

    # Region index: 1 to 9, row-major (r1 varies slowest)
    region <- (r1 - 1L) * 3L + r2

    Pr_R        <- tabulate(region, nbins = 9L) / nMC
    names(Pr_R) <- paste0("R", 1:9)

  } else {

    # -------------------------------------------------------------------------
    # Posterior predictive: simulate future multinomial data from the
    # Dirichlet-Multinomial distribution by sequential binomial draws.
    # For each MC draw u:
    #   1. Draw p_k^[u] from Dirichlet (done above).
    #   2. Draw x_ek^[u] from Multinomial(m_k, p_k^[u]) via sequential
    #      conditional binomials:
    #        x_k1 ~ Bin(m_k,        p_k1)
    #        x_k2 ~ Bin(m_k - x_k1, p_k2 / (1 - p_k1))
    #        x_k3 ~ Bin(m_k - x_k1 - x_k2, p_k3 / (1 - p_k1 - p_k2))
    #        x_k4 = m_k - x_k1 - x_k2 - x_k3
    # -------------------------------------------------------------------------

    x1_fut <- matrix(0L, nrow = nMC, ncol = 4L)
    x2_fut <- matrix(0L, nrow = nMC, ncol = 4L)
    rem1   <- rep(m1, nMC)
    rem2   <- rep(m2, nMC)
    used1  <- rep(0, nMC)
    used2  <- rep(0, nMC)

    for (j in seq_len(3L)) {
      # Conditional probability of category j given remaining categories
      denom1 <- pmax(1 - used1, 0)
      denom2 <- pmax(1 - used2, 0)
      prob1_j <- pmin(pmax(ifelse(denom1 > 0, p1[, j] / denom1, 0), 0), 1)
      prob2_j <- pmin(pmax(ifelse(denom2 > 0, p2[, j] / denom2, 0), 0), 1)

      draw1 <- rbinom(nMC, rem1, prob1_j)
      draw2 <- rbinom(nMC, rem2, prob2_j)

      x1_fut[, j] <- draw1
      x2_fut[, j] <- draw2
      rem1         <- rem1 - draw1
      rem2         <- rem2 - draw2
      used1        <- used1 + p1[, j]
      used2        <- used2 + p2[, j]
    }
    # Remaining count goes to category 4
    x1_fut[, 4L] <- rem1
    x2_fut[, 4L] <- rem2

    # Treatment effects based on future sample proportions
    # theta_e1 = (x_et10 + x_et11) / mt - (x_ec10 + x_ec11) / mc
    # theta_e2 = (x_et01 + x_et11) / mt - (x_ec01 + x_ec11) / mc
    theta1 <- (x1_fut[, 3L] + x1_fut[, 4L]) / m1 -
      (x2_fut[, 3L] + x2_fut[, 4L]) / m2
    theta2 <- (x1_fut[, 2L] + x1_fut[, 4L]) / m1 -
      (x2_fut[, 2L] + x2_fut[, 4L]) / m2

    # Classify each draw into one of 4 regions (Table 2 of reference)
    # theta.TV1 = theta.MAV1 = NULL threshold for Endpoint 1
    # r1: 1 = theta1 > NULL1, 2 = theta1 <= NULL1
    # r2: 1 = theta2 > NULL2, 2 = theta2 <= NULL2
    r1 <- 2L - as.integer(theta1 > theta.TV1)
    r2 <- 2L - as.integer(theta2 > theta.TV2)

    # Region index: row-major (r1 varies slowest)
    # R1 = (1,1), R2 = (1,2), R3 = (2,1), R4 = (2,2)
    region <- (r1 - 1L) * 2L + r2

    Pr_R        <- tabulate(region, nbins = 4L) / nMC
    names(Pr_R) <- paste0("R", 1:4)
  }

  return(Pr_R)
}
