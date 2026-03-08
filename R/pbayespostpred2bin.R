#' Bayesian Posterior or Posterior Predictive Probability for Two Binary
#' Endpoints
#'
#' Computes the Bayesian posterior probability or posterior predictive
#' probability for clinical trials with two binary endpoints under a
#' Dirichlet-multinomial conjugate model. The function returns probabilities
#' for nine decision regions (posterior) or four decision regions (predictive)
#' defined by target values (TV) and minimum acceptable values (MAV) for both
#' endpoints. Three study designs are supported: controlled, uncontrolled
#' (hypothetical control), and external (power prior).
#'
#' @param prob A character string specifying the probability type.
#'        Must be \code{'posterior'} or \code{'predictive'}.
#' @param design A character string specifying the trial design.
#'        Must be \code{'controlled'}, \code{'uncontrolled'}, or
#'        \code{'external'}.
#' @param theta_TV1 A numeric scalar giving the target value (TV) threshold
#'        for Endpoint 1. Required when \code{prob = 'posterior'}; must
#'        satisfy \code{theta_TV1 > theta_MAV1}. Set to \code{NULL} when
#'        \code{prob = 'predictive'}.
#' @param theta_MAV1 A numeric scalar giving the minimum acceptable value
#'        (MAV) threshold for Endpoint 1. Required when
#'        \code{prob = 'posterior'}; must satisfy
#'        \code{theta_TV1 > theta_MAV1}. Set to \code{NULL} when
#'        \code{prob = 'predictive'}.
#' @param theta_TV2 A numeric scalar giving the target value (TV) threshold
#'        for Endpoint 2. Required when \code{prob = 'posterior'}; must
#'        satisfy \code{theta_TV2 > theta_MAV2}. Set to \code{NULL} when
#'        \code{prob = 'predictive'}.
#' @param theta_MAV2 A numeric scalar giving the minimum acceptable value
#'        (MAV) threshold for Endpoint 2. Required when
#'        \code{prob = 'posterior'}; must satisfy
#'        \code{theta_TV2 > theta_MAV2}. Set to \code{NULL} when
#'        \code{prob = 'predictive'}.
#' @param theta_NULL1 A numeric scalar giving the null hypothesis threshold
#'        for Endpoint 1. Required when \code{prob = 'predictive'}; set to
#'        \code{NULL} when \code{prob = 'posterior'}.
#' @param theta_NULL2 A numeric scalar giving the null hypothesis threshold
#'        for Endpoint 2. Required when \code{prob = 'predictive'}; set to
#'        \code{NULL} when \code{prob = 'posterior'}.
#' @param x_t_00 A non-negative integer giving the count of
#'        \code{(0, 0)} responses in the treatment group (Endpoint 1 = 0,
#'        Endpoint 2 = 0).
#' @param x_t_01 A non-negative integer giving the count of \code{(0, 1)}
#'        responses in the treatment group.
#' @param x_t_10 A non-negative integer giving the count of \code{(1, 0)}
#'        responses in the treatment group.
#' @param x_t_11 A non-negative integer giving the count of \code{(1, 1)}
#'        responses in the treatment group.
#' @param x_c_00 A non-negative integer giving the count of \code{(0, 0)}
#'        responses in the control group. Not used when
#'        \code{design = 'uncontrolled'}; set to \code{NULL} in that case.
#' @param x_c_01 A non-negative integer giving the count of \code{(0, 1)}
#'        responses in the control group. Not used when
#'        \code{design = 'uncontrolled'}; set to \code{NULL} in that case.
#' @param x_c_10 A non-negative integer giving the count of \code{(1, 0)}
#'        responses in the control group. Not used when
#'        \code{design = 'uncontrolled'}; set to \code{NULL} in that case.
#' @param x_c_11 A non-negative integer giving the count of \code{(1, 1)}
#'        responses in the control group. Not used when
#'        \code{design = 'uncontrolled'}; set to \code{NULL} in that case.
#' @param a_t_00 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the \code{(0, 0)} response pattern in the treatment
#'        group.
#' @param a_t_01 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the \code{(0, 1)} response pattern in the treatment
#'        group.
#' @param a_t_10 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the \code{(1, 0)} response pattern in the treatment
#'        group.
#' @param a_t_11 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the \code{(1, 1)} response pattern in the treatment
#'        group.
#' @param a_c_00 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the \code{(0, 0)} response pattern in the control
#'        group. For \code{design = 'uncontrolled'}, serves as a
#'        hyperparameter of the hypothetical control distribution.
#' @param a_c_01 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the \code{(0, 1)} response pattern in the control
#'        group. For \code{design = 'uncontrolled'}, serves as a
#'        hyperparameter of the hypothetical control distribution.
#' @param a_c_10 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the \code{(1, 0)} response pattern in the control
#'        group. For \code{design = 'uncontrolled'}, serves as a
#'        hyperparameter of the hypothetical control distribution.
#' @param a_c_11 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the \code{(1, 1)} response pattern in the control
#'        group. For \code{design = 'uncontrolled'}, serves as a
#'        hyperparameter of the hypothetical control distribution.
#' @param m_t A positive integer giving the number of patients in the
#'        treatment group for the future trial. Required when
#'        \code{prob = 'predictive'}; otherwise set to \code{NULL}.
#' @param m_c A positive integer giving the number of patients in the
#'        control group for the future trial. Required when
#'        \code{prob = 'predictive'}; otherwise set to \code{NULL}.
#' @param z00 A non-negative integer giving the hypothetical control count
#'        for pattern \code{(0, 0)}. Required when
#'        \code{design = 'uncontrolled'}; otherwise set to \code{NULL}.
#' @param z01 A non-negative integer giving the hypothetical control count
#'        for pattern \code{(0, 1)}. Required when
#'        \code{design = 'uncontrolled'}; otherwise set to \code{NULL}.
#' @param z10 A non-negative integer giving the hypothetical control count
#'        for pattern \code{(1, 0)}. Required when
#'        \code{design = 'uncontrolled'}; otherwise set to \code{NULL}.
#' @param z11 A non-negative integer giving the hypothetical control count
#'        for pattern \code{(1, 1)}. Required when
#'        \code{design = 'uncontrolled'}; otherwise set to \code{NULL}.
#' @param xe_t_00 A non-negative integer giving the external treatment group
#'        count for pattern \code{(0, 0)}. Required when
#'        \code{design = 'external'} and external treatment data are used;
#'        otherwise \code{NULL}.
#' @param xe_t_01 A non-negative integer giving the external treatment group
#'        count for pattern \code{(0, 1)}. Required for external treatment
#'        data; otherwise \code{NULL}.
#' @param xe_t_10 A non-negative integer giving the external treatment group
#'        count for pattern \code{(1, 0)}. Required for external treatment
#'        data; otherwise \code{NULL}.
#' @param xe_t_11 A non-negative integer giving the external treatment group
#'        count for pattern \code{(1, 1)}. Required for external treatment
#'        data; otherwise \code{NULL}.
#' @param xe_c_00 A non-negative integer giving the external control group
#'        count for pattern \code{(0, 0)}. Required when
#'        \code{design = 'external'} and external control data are used;
#'        otherwise \code{NULL}.
#' @param xe_c_01 A non-negative integer giving the external control group
#'        count for pattern \code{(0, 1)}. Required for external control
#'        data; otherwise \code{NULL}.
#' @param xe_c_10 A non-negative integer giving the external control group
#'        count for pattern \code{(1, 0)}. Required for external control
#'        data; otherwise \code{NULL}.
#' @param xe_c_11 A non-negative integer giving the external control group
#'        count for pattern \code{(1, 1)}. Required for external control
#'        data; otherwise \code{NULL}.
#' @param alpha0e_t A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for the external treatment data. Required when external
#'        treatment data are used; otherwise set to \code{NULL}.
#' @param alpha0e_c A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for the external control data. Required when external
#'        control data are used; otherwise set to \code{NULL}.
#' @param nMC A positive integer giving the number of Monte Carlo draws used
#'        to estimate region probabilities. Default is \code{10000}. Larger
#'        values reduce Monte Carlo error at the cost of computation time.
#'
#' @return A named numeric vector of region probabilities. For
#'         \code{prob = 'posterior'}: length 9, named \code{R1}, ...,
#'         \code{R9}, corresponding to regions defined by TV and MAV
#'         thresholds for both endpoints (row-major order: Endpoint 1 varies
#'         slowest). For \code{prob = 'predictive'}: length 4, named
#'         \code{R1}, ..., \code{R4}. All elements are non-negative and sum
#'         to 1.
#'
#' @details
#' \strong{Model.} The four response categories are ordered as
#' \code{(0, 0)}, \code{(0, 1)}, \code{(1, 0)}, \code{(1, 1)}. For each
#' group \eqn{j}, the observed count vector follows a multinomial
#' distribution, and a conjugate Dirichlet prior is placed on the cell
#' probability vector \eqn{p_j}:
#' \deqn{p_j \sim \mathrm{Dir}(\alpha_{j,00}, \alpha_{j,01}, \alpha_{j,10},
#'       \alpha_{j,11}).}
#' The posterior is
#' \deqn{p_j \mid x_j \sim \mathrm{Dir}(\alpha_{j,00} + x_{j,00},
#'       \alpha_{j,01} + x_{j,01}, \alpha_{j,10} + x_{j,10},
#'       \alpha_{j,11} + x_{j,11}).}
#'
#' Marginal response rates are \eqn{\pi_{j1} = p_{j,10} + p_{j,11}}
#' (Endpoint 1) and \eqn{\pi_{j2} = p_{j,01} + p_{j,11}} (Endpoint 2).
#' Treatment effects are \eqn{\theta_1 = \pi_{t1} - \pi_{c1}} and
#' \eqn{\theta_2 = \pi_{t2} - \pi_{c2}}.
#'
#' \strong{Uncontrolled design.} When \code{design = 'uncontrolled'}, the
#' hypothetical control distribution is specified as a Dirichlet
#' distribution directly via the prior hyperparameters and hypothetical
#' control counts \code{z00}, \code{z01}, \code{z10}, \code{z11}.
#'
#' \strong{External design.} When \code{design = 'external'}, the power
#' prior augments the Dirichlet prior parameters.
#'
#' @examples
#' # Example 1: Controlled design - posterior probability
#' pbayespostpred2bin(
#'   prob = 'posterior', design = 'controlled',
#'   theta_TV1 = 0.15, theta_MAV1 = 0.05, theta_TV2 = 0.10, theta_MAV2 = 0.0,
#'   theta_NULL1 = NULL, theta_NULL2 = NULL,
#'   x_t_00 = 2, x_t_01 = 3, x_t_10 = 5, x_t_11 = 2,
#'   x_c_00 = 3, x_c_01 = 2, x_c_10 = 4, x_c_11 = 1,
#'   a_t_00 = 1, a_t_01 = 1, a_t_10 = 1, a_t_11 = 1,
#'   a_c_00 = 1, a_c_01 = 1, a_c_10 = 1, a_c_11 = 1,
#'   m_t = NULL, m_c = NULL,
#'   z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
#'   xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
#'   xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
#'   alpha0e_t = NULL, alpha0e_c = NULL
#' )
#'
#' # Example 2: Uncontrolled design - posterior probability
#' pbayespostpred2bin(
#'   prob = 'posterior', design = 'uncontrolled',
#'   theta_TV1 = 0.15, theta_MAV1 = 0.05, theta_TV2 = 0.10, theta_MAV2 = 0.0,
#'   theta_NULL1 = NULL, theta_NULL2 = NULL,
#'   x_t_00 = 2, x_t_01 = 3, x_t_10 = 5, x_t_11 = 2,
#'   x_c_00 = NULL, x_c_01 = NULL, x_c_10 = NULL, x_c_11 = NULL,
#'   a_t_00 = 1, a_t_01 = 1, a_t_10 = 1, a_t_11 = 1,
#'   a_c_00 = 1, a_c_01 = 1, a_c_10 = 1, a_c_11 = 1,
#'   m_t = NULL, m_c = NULL,
#'   z00 = 1, z01 = 2, z10 = 2, z11 = 1,
#'   xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
#'   xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
#'   alpha0e_t = NULL, alpha0e_c = NULL
#' )
#'
#' # Example 3: External design - posterior probability
#' pbayespostpred2bin(
#'   prob = 'posterior', design = 'external',
#'   theta_TV1 = 0.15, theta_MAV1 = 0.05, theta_TV2 = 0.10, theta_MAV2 = 0.0,
#'   theta_NULL1 = NULL, theta_NULL2 = NULL,
#'   x_t_00 = 2, x_t_01 = 3, x_t_10 = 5, x_t_11 = 2,
#'   x_c_00 = 3, x_c_01 = 2, x_c_10 = 4, x_c_11 = 1,
#'   a_t_00 = 1, a_t_01 = 1, a_t_10 = 1, a_t_11 = 1,
#'   a_c_00 = 1, a_c_01 = 1, a_c_10 = 1, a_c_11 = 1,
#'   m_t = NULL, m_c = NULL,
#'   z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
#'   xe_t_00 = 2, xe_t_01 = 2, xe_t_10 = 3, xe_t_11 = 1,
#'   xe_c_00 = 2, xe_c_01 = 1, xe_c_10 = 2, xe_c_11 = 1,
#'   alpha0e_t = 0.5, alpha0e_c = 0.5
#' )
#'
#' # Example 4: Controlled design - posterior predictive probability
#' pbayespostpred2bin(
#'   prob = 'predictive', design = 'controlled',
#'   theta_TV1 = NULL, theta_MAV1 = NULL, theta_TV2 = NULL, theta_MAV2 = NULL,
#'   theta_NULL1 = 0.05, theta_NULL2 = 0.0,
#'   x_t_00 = 2, x_t_01 = 3, x_t_10 = 5, x_t_11 = 2,
#'   x_c_00 = 3, x_c_01 = 2, x_c_10 = 4, x_c_11 = 1,
#'   a_t_00 = 1, a_t_01 = 1, a_t_10 = 1, a_t_11 = 1,
#'   a_c_00 = 1, a_c_01 = 1, a_c_10 = 1, a_c_11 = 1,
#'   m_t = 20, m_c = 20,
#'   z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
#'   xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
#'   xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
#'   alpha0e_t = NULL, alpha0e_c = NULL
#' )
#'
#' # Example 5: Uncontrolled design - posterior predictive probability
#' pbayespostpred2bin(
#'   prob = 'predictive', design = 'uncontrolled',
#'   theta_TV1 = NULL, theta_MAV1 = NULL, theta_TV2 = NULL, theta_MAV2 = NULL,
#'   theta_NULL1 = 0.05, theta_NULL2 = 0.0,
#'   x_t_00 = 2, x_t_01 = 3, x_t_10 = 5, x_t_11 = 2,
#'   x_c_00 = NULL, x_c_01 = NULL, x_c_10 = NULL, x_c_11 = NULL,
#'   a_t_00 = 1, a_t_01 = 1, a_t_10 = 1, a_t_11 = 1,
#'   a_c_00 = 1, a_c_01 = 1, a_c_10 = 1, a_c_11 = 1,
#'   m_t = 20, m_c = 20,
#'   z00 = 1, z01 = 2, z10 = 2, z11 = 1,
#'   xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
#'   xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
#'   alpha0e_t = NULL, alpha0e_c = NULL
#' )
#'
#' # Example 6: External design - posterior predictive probability
#' pbayespostpred2bin(
#'   prob = 'predictive', design = 'external',
#'   theta_TV1 = NULL, theta_MAV1 = NULL, theta_TV2 = NULL, theta_MAV2 = NULL,
#'   theta_NULL1 = 0.05, theta_NULL2 = 0.0,
#'   x_t_00 = 2, x_t_01 = 3, x_t_10 = 5, x_t_11 = 2,
#'   x_c_00 = 3, x_c_01 = 2, x_c_10 = 4, x_c_11 = 1,
#'   a_t_00 = 1, a_t_01 = 1, a_t_10 = 1, a_t_11 = 1,
#'   a_c_00 = 1, a_c_01 = 1, a_c_10 = 1, a_c_11 = 1,
#'   m_t = 20, m_c = 20,
#'   z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
#'   xe_t_00 = 2, xe_t_01 = 2, xe_t_10 = 3, xe_t_11 = 1,
#'   xe_c_00 = 2, xe_c_01 = 1, xe_c_10 = 2, xe_c_11 = 1,
#'   alpha0e_t = 0.5, alpha0e_c = 0.5
#' )
#'
#' @export
pbayespostpred2bin <- function(prob    = 'posterior',
                               design  = 'controlled',
                               theta_TV1   = NULL, theta_MAV1  = NULL,
                               theta_TV2   = NULL, theta_MAV2  = NULL,
                               theta_NULL1 = NULL, theta_NULL2 = NULL,
                               x_t_00, x_t_01, x_t_10, x_t_11,
                               x_c_00  = NULL, x_c_01  = NULL,
                               x_c_10  = NULL, x_c_11  = NULL,
                               a_t_00, a_t_01, a_t_10, a_t_11,
                               a_c_00, a_c_01, a_c_10, a_c_11,
                               m_t = NULL, m_c = NULL,
                               z00 = NULL, z01 = NULL,
                               z10 = NULL, z11 = NULL,
                               xe_t_00 = NULL, xe_t_01 = NULL,
                               xe_t_10 = NULL, xe_t_11 = NULL,
                               xe_c_00 = NULL, xe_c_01 = NULL,
                               xe_c_10 = NULL, xe_c_11 = NULL,
                               alpha0e_t = NULL, alpha0e_c = NULL,
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
    for (nm in c("theta_TV1", "theta_MAV1", "theta_TV2", "theta_MAV2")) {
      val <- get(nm)
      if (is.null(val)) {
        stop(paste0("'", nm, "' must be non-NULL when prob = 'posterior'"))
      }
      if (!is.numeric(val) || length(val) != 1L || is.na(val)) {
        stop(paste0("'", nm, "' must be a single numeric value"))
      }
    }
    if (theta_TV1 <= theta_MAV1) {
      stop("'theta_TV1' must be strictly greater than 'theta_MAV1'")
    }
    if (theta_TV2 <= theta_MAV2) {
      stop("'theta_TV2' must be strictly greater than 'theta_MAV2'")
    }
  } else {
    for (nm in c("theta_NULL1", "theta_NULL2")) {
      val <- get(nm)
      if (is.null(val)) {
        stop(paste0("'", nm, "' must be non-NULL when prob = 'predictive'"))
      }
      if (!is.numeric(val) || length(val) != 1L || is.na(val)) {
        stop(paste0("'", nm, "' must be a single numeric value"))
      }
    }
    # Convert null thresholds: TV = MAV = NULL for predictive probability
    theta_TV1  <- theta_NULL1;  theta_MAV1 <- theta_NULL1
    theta_TV2  <- theta_NULL2;  theta_MAV2 <- theta_NULL2
  }

  # --- Group 1 observed counts (always required) ---
  for (nm in c("x_t_00", "x_t_01", "x_t_10", "x_t_11")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
        val != floor(val) || val < 0L) {
      stop(paste0("'", nm, "' must be a single non-negative integer"))
    }
  }

  # --- Group 2 observed counts (required unless design = 'uncontrolled') ---
  if (design != 'uncontrolled') {
    for (nm in c("x_c_00", "x_c_01", "x_c_10", "x_c_11")) {
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
    # Warn if the control group data are accidentally supplied for uncontrolled design
    if (!is.null(x_c_00) || !is.null(x_c_01) ||
        !is.null(x_c_10) || !is.null(x_c_11)) {
      warning(paste0("design = 'uncontrolled': 'x_c_00', 'x_c_01', 'x_c_10', ",
                     "'x_c_11' are ignored. Use 'z00', 'z01', 'z10', 'z11' ",
                     "to specify the hypothetical control counts."))
    }
  }

  # --- Dirichlet prior parameters (always required) ---
  for (nm in c("a_t_00", "a_t_01", "a_t_10", "a_t_11",
               "a_c_00", "a_c_01", "a_c_10", "a_c_11")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) || val <= 0) {
      stop(paste0("'", nm, "' must be a single positive numeric value"))
    }
  }

  # --- Future sample sizes (required for predictive) ---
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
    m_t <- as.integer(m_t)
    m_c <- as.integer(m_c)
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
    has_ext1 <- !is.null(xe_t_00) && !is.null(xe_t_01) &&
      !is.null(xe_t_10) && !is.null(xe_t_11) && !is.null(alpha0e_t)
    has_ext2 <- !is.null(xe_c_00) && !is.null(xe_c_01) &&
      !is.null(xe_c_10) && !is.null(xe_c_11) && !is.null(alpha0e_c)

    if (!has_ext1 && !has_ext2) {
      stop(paste0("For design = 'external', at least one complete set of ",
                  "external data (xe1_* + alpha0e_t or xe2_* + alpha0e_c) must be provided"))
    }

    for (nm in c("alpha0e_t", "alpha0e_c")) {
      val <- get(nm)
      if (!is.null(val)) {
        if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
            val <= 0 || val > 1) {
          stop(paste0("'", nm, "' must be a single numeric value in (0, 1]"))
        }
      }
    }

    ext1_nms <- c("xe_t_00", "xe_t_01", "xe_t_10", "xe_t_11")
    ext2_nms <- c("xe_c_00", "xe_c_01", "xe_c_10", "xe_c_11")

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
  xe_t_weight <- if (!is.null(alpha0e_t) && design == 'external') alpha0e_t else 0
  alpha_t <- c(
    a_t_00 + x_t_00 + xe_t_weight * ifelse(!is.null(xe_t_00), xe_t_00, 0),
    a_t_01 + x_t_01 + xe_t_weight * ifelse(!is.null(xe_t_01), xe_t_01, 0),
    a_t_10 + x_t_10 + xe_t_weight * ifelse(!is.null(xe_t_10), xe_t_10, 0),
    a_t_11 + x_t_11 + xe_t_weight * ifelse(!is.null(xe_t_11), xe_t_11, 0)
  )

  # Group 2: parameterisation depends on design
  if (design == 'uncontrolled') {
    # Hypothetical control distribution: prior augmented by pseudo-counts z*
    # (Appendix A1.1.6: P_hyp(p_c) = Dir(beta_c00, ..., beta_c11) where
    #  beta_c** = a2_** + z**)
    alpha_c <- c(
      a_c_00 + z00,
      a_c_01 + z01,
      a_c_10 + z10,
      a_c_11 + z11
    )
  } else {
    # Controlled or external: posterior from observed control data
    xe_c_weight <- if (!is.null(alpha0e_c) && design == 'external') alpha0e_c else 0
    alpha_c <- c(
      a_c_00 + x_c_00 + xe_c_weight * ifelse(!is.null(xe_c_00), xe_c_00, 0),
      a_c_01 + x_c_01 + xe_c_weight * ifelse(!is.null(xe_c_01), xe_c_01, 0),
      a_c_10 + x_c_10 + xe_c_weight * ifelse(!is.null(xe_c_10), xe_c_10, 0),
      a_c_11 + x_c_11 + xe_c_weight * ifelse(!is.null(xe_c_11), xe_c_11, 0)
    )
  }

  # ---------------------------------------------------------------------------
  # Section 3: Monte Carlo estimation
  # ---------------------------------------------------------------------------

  # Draw nMC samples from the posterior (or hypothetical) Dirichlet for each group
  # Columns of p_t/p_c: (p_00, p_01, p_10, p_11)
  p_t <- rdirichlet(nMC, alpha_t)
  p_c <- rdirichlet(nMC, alpha_c)

  if (prob == 'posterior') {

    # Marginal response rates:
    #   pi_k1 = p_k10 + p_k11  (Endpoint 1 success)
    #   pi_k2 = p_k01 + p_k11  (Endpoint 2 success)
    theta1 <- (p_t[, 3L] + p_t[, 4L]) - (p_c[, 3L] + p_c[, 4L])
    theta2 <- (p_t[, 2L] + p_t[, 4L]) - (p_c[, 2L] + p_c[, 4L])

    # Classify each draw into one of 9 regions (row-major: Endpoint 1 slowest)
    # r1: 1 = theta1 > TV1, 2 = MAV1 < theta1 <= TV1, 3 = theta1 <= MAV1
    # r2: 1 = theta2 > TV2, 2 = MAV2 < theta2 <= TV2, 3 = theta2 <= MAV2
    r1 <- 3L -
      as.integer(theta1 > theta_MAV1) -
      as.integer(theta1 > theta_TV1)
    r2 <- 3L -
      as.integer(theta2 > theta_MAV2) -
      as.integer(theta2 > theta_TV2)

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

    x_t_fut <- matrix(0L, nrow = nMC, ncol = 4L)
    x_c_fut <- matrix(0L, nrow = nMC, ncol = 4L)
    rem_t   <- rep(m_t, nMC)
    rem_c   <- rep(m_c, nMC)
    used_t  <- rep(0, nMC)
    used_c  <- rep(0, nMC)

    for (j in seq_len(3L)) {
      # Conditional probability of category j given remaining categories
      denom1 <- pmax(1 - used_t, 0)
      denom2 <- pmax(1 - used_c, 0)
      prob1_j <- pmin(pmax(ifelse(denom1 > 0, p_t[, j] / denom1, 0), 0), 1)
      prob2_j <- pmin(pmax(ifelse(denom2 > 0, p_c[, j] / denom2, 0), 0), 1)

      draw1 <- rbinom(nMC, rem_t, prob1_j)
      draw2 <- rbinom(nMC, rem_c, prob2_j)

      x_t_fut[, j] <- draw1
      x_c_fut[, j] <- draw2
      rem_t         <- rem_t - draw1
      rem_c         <- rem_c - draw2
      used_t        <- used_t + p_t[, j]
      used_c        <- used_c + p_c[, j]
    }
    # Remaining count goes to category 4
    x_t_fut[, 4L] <- rem_t
    x_c_fut[, 4L] <- rem_c

    # Treatment effects based on future sample proportions
    # theta_e1 = (x_et10 + x_et11) / mt - (x_ec10 + x_ec11) / mc
    # theta_e2 = (x_et01 + x_et11) / mt - (x_ec01 + x_ec11) / mc
    theta1 <- (x_t_fut[, 3L] + x_t_fut[, 4L]) / m_t -
      (x_c_fut[, 3L] + x_c_fut[, 4L]) / m_c
    theta2 <- (x_t_fut[, 2L] + x_t_fut[, 4L]) / m_t -
      (x_c_fut[, 2L] + x_c_fut[, 4L]) / m_c

    # Classify each draw into one of 4 regions (Table 2 of reference)
    # theta_TV1 = theta_MAV1 = NULL threshold for Endpoint 1
    # r1: 1 = theta1 > NULL1, 2 = theta1 <= NULL1
    # r2: 1 = theta2 > NULL2, 2 = theta2 <= NULL2
    r1 <- 2L - as.integer(theta1 > theta_TV1)
    r2 <- 2L - as.integer(theta2 > theta_TV2)

    # Region index: row-major (r1 varies slowest)
    # R1 = (1,1), R2 = (1,2), R3 = (2,1), R4 = (2,2)
    region <- (r1 - 1L) * 2L + r2

    Pr_R        <- tabulate(region, nbins = 4L) / nMC
    names(Pr_R) <- paste0("R", 1:4)
  }

  return(Pr_R)
}
