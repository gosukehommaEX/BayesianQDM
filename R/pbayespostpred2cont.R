#' Bayesian Posterior or Posterior Predictive Probability for a Clinical Trial
#' with Two Continuous Endpoints
#'
#' Computes the Bayesian posterior probability or posterior predictive
#' probability for clinical trials with two continuous endpoints under a
#' bivariate Normal model. Sufficient statistics of the observed PoC data
#' are supplied directly. The function returns probabilities for nine decision
#' regions (posterior) or four decision regions (predictive) defined by target
#' values (TV) and minimum acceptable values (MAV) for both endpoints. Three
#' study designs are supported: controlled, uncontrolled (hypothetical
#' control), and external-control (power prior with NIW conjugate
#' representation).
#'
#' @param prob A character string specifying the probability type.
#'        Must be \code{'posterior'} or \code{'predictive'}.
#' @param design A character string specifying the trial design.
#'        Must be \code{'controlled'}, \code{'uncontrolled'}, or
#'        \code{'external'}.
#' @param prior A character string specifying the prior distribution.
#'        Must be \code{'vague'} or \code{'N-Inv-Wishart'}.
#' @param theta.TV1 A numeric scalar giving the target value (TV) for
#'        Endpoint 1. Required when \code{prob = 'posterior'};
#'        set to \code{NULL} for \code{prob = 'predictive'}.
#' @param theta.MAV1 A numeric scalar giving the minimum acceptable value
#'        (MAV) for Endpoint 1. Required when \code{prob = 'posterior'};
#'        must satisfy \code{theta.TV1 > theta.MAV1}.
#'        Set to \code{NULL} for \code{prob = 'predictive'}.
#' @param theta.TV2 A numeric scalar giving the target value (TV) for
#'        Endpoint 2. Required when \code{prob = 'posterior'};
#'        set to \code{NULL} for \code{prob = 'predictive'}.
#' @param theta.MAV2 A numeric scalar giving the minimum acceptable value
#'        (MAV) for Endpoint 2. Required when \code{prob = 'posterior'};
#'        must satisfy \code{theta.TV2 > theta.MAV2}.
#'        Set to \code{NULL} for \code{prob = 'predictive'}.
#' @param theta.NULL1 A numeric scalar giving the null hypothesis value for
#'        Endpoint 1. Required when \code{prob = 'predictive'};
#'        set to \code{NULL} for \code{prob = 'posterior'}.
#' @param theta.NULL2 A numeric scalar giving the null hypothesis value for
#'        Endpoint 2. Required when \code{prob = 'predictive'};
#'        set to \code{NULL} for \code{prob = 'posterior'}.
#' @param n1 A positive integer giving the PoC sample size for group 1
#'        (treatment). Must be >= 3 when \code{prior = 'vague'}.
#' @param n2 A positive integer giving the PoC sample size for group 2
#'        (control). Required for \code{design = 'controlled'} or
#'        \code{'external'}; set to \code{NULL} for \code{'uncontrolled'}.
#'        Must be >= 3 when \code{prior = 'vague'}.
#' @param ybar1 A numeric vector of length 2 giving the sample mean
#'        \eqn{(\bar{y}_{t1}, \bar{y}_{t2})} for group 1.
#' @param S1 A 2x2 positive-definite numeric matrix giving the
#'        sum-of-squares matrix
#'        \eqn{S_t = \sum_j (y_{t,j} - \bar{y}_t)(y_{t,j} - \bar{y}_t)^T}
#'        for group 1.
#' @param ybar2 A numeric vector of length 2 giving the sample mean for
#'        group 2. Required for \code{design = 'controlled'} or
#'        \code{'external'}; set to \code{NULL} for \code{'uncontrolled'}.
#' @param S2 A 2x2 positive-definite numeric matrix giving the
#'        sum-of-squares matrix for group 2. Same availability rules as
#'        \code{ybar2}.
#' @param m1 A positive integer giving the future sample size for group 1.
#'        Required when \code{prob = 'predictive'}; otherwise \code{NULL}.
#' @param m2 A positive integer giving the future sample size for group 2.
#'        Required when \code{prob = 'predictive'}; otherwise \code{NULL}.
#' @param kappa01 A positive numeric scalar. Prior precision parameter for
#'        group 1 when \code{prior = 'N-Inv-Wishart'}; otherwise \code{NULL}.
#' @param nu01 A positive numeric scalar (> 2). Prior degrees of freedom for
#'        group 1 when \code{prior = 'N-Inv-Wishart'}; otherwise \code{NULL}.
#' @param mu01 A numeric vector of length 2. Prior mean for group 1 when
#'        \code{prior = 'N-Inv-Wishart'}; otherwise \code{NULL}.
#' @param Lambda01 A 2x2 positive-definite numeric matrix. Prior scale
#'        matrix for group 1 when \code{prior = 'N-Inv-Wishart'};
#'        otherwise \code{NULL}.
#' @param kappa02 A positive numeric scalar. Prior precision parameter for
#'        group 2 when \code{prior = 'N-Inv-Wishart'}; otherwise \code{NULL}.
#' @param nu02 A positive numeric scalar (> 2). Prior degrees of freedom for
#'        group 2 when \code{prior = 'N-Inv-Wishart'}; otherwise \code{NULL}.
#' @param mu02 A numeric vector of length 2. Prior mean for group 2 when
#'        \code{prior = 'N-Inv-Wishart'}, or the hypothetical control
#'        location when \code{design = 'uncontrolled'}; otherwise \code{NULL}.
#' @param Lambda02 A 2x2 positive-definite numeric matrix. Prior scale
#'        matrix for group 2 when \code{prior = 'N-Inv-Wishart'};
#'        otherwise \code{NULL}.
#' @param r A positive numeric scalar. Variance scaling factor for the
#'        hypothetical control distribution. Required when
#'        \code{design = 'uncontrolled'}; otherwise \code{NULL}.
#'        When \code{r = 1}, the control and treatment have equal variance
#'        scale. The scale matrix of the hypothetical control is set to
#'        \code{r} times the treatment posterior scale matrix.
#' @param ne1 A positive integer giving the external sample size for group 1.
#'        Required when \code{design = 'external'} and external treatment
#'        data are used; otherwise \code{NULL}.
#' @param ne2 A positive integer giving the external sample size for group 2.
#'        Required when \code{design = 'external'} and external control data
#'        are used; otherwise \code{NULL}.
#' @param alpha01e A numeric scalar in \code{(0, 1]}. Power prior weight for
#'        group 1 external data. Required when \code{ne1} is non-NULL;
#'        otherwise \code{NULL}.
#' @param alpha02e A numeric scalar in \code{(0, 1]}. Power prior weight for
#'        group 2 external data. Required when \code{ne2} is non-NULL;
#'        otherwise \code{NULL}.
#' @param ybar_e1 A numeric vector of length 2. Sample mean of external data
#'        for group 1. Required when \code{ne1} is non-NULL; otherwise
#'        \code{NULL}.
#' @param ybar_e2 A numeric vector of length 2. Sample mean of external data
#'        for group 2. Required when \code{ne2} is non-NULL; otherwise
#'        \code{NULL}.
#' @param Se1 A 2x2 positive-definite numeric matrix. Sum-of-squares matrix
#'        of external data for group 1. Required when \code{ne1} is
#'        non-NULL; otherwise \code{NULL}.
#' @param Se2 A 2x2 positive-definite numeric matrix. Sum-of-squares matrix
#'        of external data for group 2. Required when \code{ne2} is
#'        non-NULL; otherwise \code{NULL}.
#' @param nMC A positive integer giving the number of Monte Carlo samples.
#'        Default is \code{10000L}.  Used when \code{method = 'MC'} or when
#'        \code{method = 'MM'} falls back to MC.
#' @param method A character string specifying the computation method.
#'        Must be \code{'MC'} (Monte Carlo; default) or \code{'MM'}
#'        (Moment-Matching approximation based on Theorem 3 of the reference).
#'        When \code{method = 'MM'}, the parameters of the approximating
#'        bivariate t-distribution are determined analytically via moment
#'        matching, and rectangular region probabilities are evaluated with
#'        \code{mvtnorm::pmvt}.  If the degrees-of-freedom condition
#'        \eqn{\nu_k > 4} is not met (e.g., vague prior with small
#'        \eqn{n_k \leq 6}), a warning is issued and the function falls back
#'        to \code{method = 'MC'}.
#'
#' @return A named numeric vector of region probabilities.
#'         For \code{prob = 'posterior'}: length 9, named \code{R1}, ...,
#'         \code{R9}, corresponding to a row-major 3x3 grid where Endpoint 1
#'         varies slowest (see Details for region definitions).
#'         For \code{prob = 'predictive'}: length 4, named \code{R1}, ...,
#'         \code{R4}.
#'         All elements are non-negative and sum to 1.
#'
#' @details
#' \strong{Model.}
#' Both endpoints follow a bivariate Normal distribution
#' \eqn{y_{k,j} \sim N_2(\mu_k, \Sigma_k)} for arm \eqn{k \in \{t, c\}}.
#' The treatment effect is \eqn{\theta = \mu_t - \mu_c}.
#'
#' \strong{Posterior distribution (vague prior).}
#' \deqn{\mu_k | Y_k \sim t_{n_k - 2}\!\left(\bar{y}_k,\;
#'       \frac{S_k}{n_k(n_k - 2)}\right)}
#'
#' \strong{Posterior distribution (NIW prior).}
#' \deqn{\mu_k | Y_k \sim t_{\nu_{nk} - 1}\!\left(\mu_{nk},\;
#'       \frac{\Lambda_{nk}}{\kappa_{nk}(\nu_{nk} - 1)}\right)}
#' with updated hyperparameters
#' \eqn{\kappa_{nk} = \kappa_{0k} + n_k},
#' \eqn{\nu_{nk} = \nu_{0k} + n_k},
#' \eqn{\mu_{nk} = (\kappa_{0k}\mu_{0k} + n_k\bar{y}_k)/\kappa_{nk}}, and
#' \eqn{\Lambda_{nk} = \Lambda_{0k} + S_k +
#' \kappa_{0k}n_k(\bar{y}_k - \mu_{0k})(\bar{y}_k - \mu_{0k})^T / \kappa_{nk}}.
#'
#' \strong{Predictive distribution.}
#' The scale matrix of a single future observation is inflated by
#' \eqn{(1 + n_k)/n_k} (vague) or \eqn{(1 + \kappa_{nk})/\kappa_{nk}} (NIW)
#' relative to the posterior. The mean of \eqn{m_k} future observations has
#' scale divided by \eqn{m_k}.
#'
#' \strong{Posterior probability regions (prob = 'posterior').}
#' Row-major 3x3 grid; Endpoint 1 varies slowest:
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
#' Row-major 2x2 grid; Endpoint 1 varies slowest:
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
#' \strong{Uncontrolled design.}
#' The hypothetical control distribution replaces the control arm posterior.
#' Under vague prior:
#' \eqn{\mu_c \sim t_{n_t - 2}(\mu_{0c},\; rS_t / [n_t(n_t - 2)])}.
#' Under NIW prior:
#' \eqn{\mu_c \sim t_{\nu_{nt} - 1}(\mu_{0c},\;
#' r\Lambda_{nt} / [\kappa_{nt}(\nu_{nt} - 1)])}.
#' The parameter \code{r} allows the variance scale of the hypothetical
#' control to differ from the treatment arm.
#'
#' \strong{External design.}
#' Incorporated via the power prior with NIW conjugate representation.
#' The effective posterior hyperparameters are obtained by constructing the
#' power prior from external data with weight \eqn{a_0}, then updating with
#' current PoC data (see Conjugate_power_prior.pdf, Theorem 5).
#'
#' @examples
#' # Example 1: Posterior probability, controlled design, vague prior
#' S1 <- matrix(c(18.0, 3.6, 3.6, 9.0), 2, 2)
#' S2 <- matrix(c(16.0, 2.8, 2.8, 8.5), 2, 2)
#' pbayespostpred2cont(
#'   prob = 'posterior', design = 'controlled', prior = 'vague',
#'   theta.TV1 = 1.5, theta.MAV1 = 0.5,
#'   theta.TV2 = 1.0, theta.MAV2 = 0.3,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   n1 = 12L, n2 = 12L,
#'   ybar1 = c(3.5, 2.1), S1 = S1,
#'   ybar2 = c(1.8, 1.0), S2 = S2,
#'   m1 = NULL, m2 = NULL,
#'   kappa01 = NULL, nu01 = NULL, mu01 = NULL, Lambda01 = NULL,
#'   kappa02 = NULL, nu02 = NULL, mu02 = NULL, Lambda02 = NULL,
#'   r = NULL,
#'   ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
#'   ybar_e1 = NULL, ybar_e2 = NULL, Se1 = NULL, Se2 = NULL,
#'   nMC = 1000L
#' )
#'
#' # Example 2: Posterior probability, controlled design, NIW prior
#' S1 <- matrix(c(18.0, 3.6, 3.6, 9.0), 2, 2)
#' S2 <- matrix(c(16.0, 2.8, 2.8, 8.5), 2, 2)
#' L0 <- matrix(c(20.0, 0.0, 0.0, 10.0), 2, 2)
#' pbayespostpred2cont(
#'   prob = 'posterior', design = 'controlled', prior = 'N-Inv-Wishart',
#'   theta.TV1 = 1.5, theta.MAV1 = 0.5,
#'   theta.TV2 = 1.0, theta.MAV2 = 0.3,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   n1 = 12L, n2 = 12L,
#'   ybar1 = c(3.5, 2.1), S1 = S1,
#'   ybar2 = c(1.8, 1.0), S2 = S2,
#'   m1 = NULL, m2 = NULL,
#'   kappa01 = 2.0, nu01 = 5.0, mu01 = c(2.0, 1.0), Lambda01 = L0,
#'   kappa02 = 2.0, nu02 = 5.0, mu02 = c(1.0, 0.5), Lambda02 = L0,
#'   r = NULL,
#'   ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
#'   ybar_e1 = NULL, ybar_e2 = NULL, Se1 = NULL, Se2 = NULL,
#'   nMC = 1000L
#' )
#'
#' # Example 3: Posterior probability, uncontrolled design, NIW prior
#' # mu02 specifies the hypothetical control location; r scales the variance
#' S1 <- matrix(c(18.0, 3.6, 3.6, 9.0), 2, 2)
#' L0 <- matrix(c(20.0, 0.0, 0.0, 10.0), 2, 2)
#' pbayespostpred2cont(
#'   prob = 'posterior', design = 'uncontrolled', prior = 'N-Inv-Wishart',
#'   theta.TV1 = 1.5, theta.MAV1 = 0.5,
#'   theta.TV2 = 1.0, theta.MAV2 = 0.3,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   n1 = 12L, n2 = NULL,
#'   ybar1 = c(3.5, 2.1), S1 = S1,
#'   ybar2 = NULL, S2 = NULL,
#'   m1 = NULL, m2 = NULL,
#'   kappa01 = 2.0, nu01 = 5.0, mu01 = c(2.0, 1.0), Lambda01 = L0,
#'   kappa02 = NULL, nu02 = NULL, mu02 = c(1.0, 0.5), Lambda02 = NULL,
#'   r = 1.0,
#'   ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
#'   ybar_e1 = NULL, ybar_e2 = NULL, Se1 = NULL, Se2 = NULL,
#'   nMC = 1000L
#' )
#'
#' # Example 4: Posterior probability, external design (control only), NIW prior
#' S1  <- matrix(c(18.0, 3.6, 3.6, 9.0), 2, 2)
#' S2  <- matrix(c(16.0, 2.8, 2.8, 8.5), 2, 2)
#' L0  <- matrix(c(20.0, 0.0, 0.0, 10.0), 2, 2)
#' Se2 <- matrix(c(15.0, 2.5, 2.5, 7.5), 2, 2)
#' pbayespostpred2cont(
#'   prob = 'posterior', design = 'external', prior = 'N-Inv-Wishart',
#'   theta.TV1 = 1.5, theta.MAV1 = 0.5,
#'   theta.TV2 = 1.0, theta.MAV2 = 0.3,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   n1 = 12L, n2 = 12L,
#'   ybar1 = c(3.5, 2.1), S1 = S1,
#'   ybar2 = c(1.8, 1.0), S2 = S2,
#'   m1 = NULL, m2 = NULL,
#'   kappa01 = 2.0, nu01 = 5.0, mu01 = c(2.0, 1.0), Lambda01 = L0,
#'   kappa02 = 2.0, nu02 = 5.0, mu02 = c(1.0, 0.5), Lambda02 = L0,
#'   r = NULL,
#'   ne1 = NULL, ne2 = 10L, alpha01e = NULL, alpha02e = 0.5,
#'   ybar_e1 = NULL, ybar_e2 = c(1.5, 0.8), Se1 = NULL, Se2 = Se2,
#'   nMC = 1000L
#' )
#'
#' # Example 5: Predictive probability, controlled design, vague prior
#' S1 <- matrix(c(18.0, 3.6, 3.6, 9.0), 2, 2)
#' S2 <- matrix(c(16.0, 2.8, 2.8, 8.5), 2, 2)
#' pbayespostpred2cont(
#'   prob = 'predictive', design = 'controlled', prior = 'vague',
#'   theta.TV1 = NULL, theta.MAV1 = NULL,
#'   theta.TV2 = NULL, theta.MAV2 = NULL,
#'   theta.NULL1 = 0.5, theta.NULL2 = 0.3,
#'   n1 = 12L, n2 = 12L,
#'   ybar1 = c(3.5, 2.1), S1 = S1,
#'   ybar2 = c(1.8, 1.0), S2 = S2,
#'   m1 = 30L, m2 = 30L,
#'   kappa01 = NULL, nu01 = NULL, mu01 = NULL, Lambda01 = NULL,
#'   kappa02 = NULL, nu02 = NULL, mu02 = NULL, Lambda02 = NULL,
#'   r = NULL,
#'   ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
#'   ybar_e1 = NULL, ybar_e2 = NULL, Se1 = NULL, Se2 = NULL,
#'   nMC = 1000L
#' )
#'
#' # Example 6: Posterior probability, controlled design, NIW prior, MM method
#' S1 <- matrix(c(18.0, 3.6, 3.6, 9.0), 2, 2)
#' S2 <- matrix(c(16.0, 2.8, 2.8, 8.5), 2, 2)
#' L0 <- matrix(c(20.0, 0.0, 0.0, 10.0), 2, 2)
#' pbayespostpred2cont(
#'   prob = 'posterior', design = 'controlled', prior = 'N-Inv-Wishart',
#'   theta.TV1 = 1.5, theta.MAV1 = 0.5,
#'   theta.TV2 = 1.0, theta.MAV2 = 0.3,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   n1 = 12L, n2 = 12L,
#'   ybar1 = c(3.5, 2.1), S1 = S1,
#'   ybar2 = c(1.8, 1.0), S2 = S2,
#'   m1 = NULL, m2 = NULL,
#'   kappa01 = 2.0, nu01 = 5.0, mu01 = c(2.0, 1.0), Lambda01 = L0,
#'   kappa02 = 2.0, nu02 = 5.0, mu02 = c(1.0, 0.5), Lambda02 = L0,
#'   r = NULL,
#'   ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
#'   ybar_e1 = NULL, ybar_e2 = NULL, Se1 = NULL, Se2 = NULL,
#'   nMC = 1000L, method = 'MM'
#' )
#'
#' @importFrom stats rnorm rchisq
#' @importFrom mvtnorm pmvt
#' @export
pbayespostpred2cont <- function(prob,
                                design,
                                prior,
                                theta.TV1   = NULL, theta.MAV1  = NULL,
                                theta.TV2   = NULL, theta.MAV2  = NULL,
                                theta.NULL1 = NULL, theta.NULL2 = NULL,
                                n1, n2      = NULL,
                                ybar1, S1,
                                ybar2       = NULL, S2          = NULL,
                                m1          = NULL, m2          = NULL,
                                kappa01     = NULL, nu01        = NULL,
                                mu01        = NULL, Lambda01    = NULL,
                                kappa02     = NULL, nu02        = NULL,
                                mu02        = NULL, Lambda02    = NULL,
                                r           = NULL,
                                ne1         = NULL, ne2         = NULL,
                                alpha01e    = NULL, alpha02e    = NULL,
                                ybar_e1     = NULL, ybar_e2     = NULL,
                                Se1         = NULL, Se2         = NULL,
                                nMC         = 10000L,
                                method      = 'MC') {

  # ---------------------------------------------------------------------------
  # Section 1: Input validation
  # ---------------------------------------------------------------------------

  if (!is.character(prob) || length(prob) != 1L ||
      !prob %in% c('posterior', 'predictive'))
    stop("'prob' must be either 'posterior' or 'predictive'")

  if (!is.character(design) || length(design) != 1L ||
      !design %in% c('controlled', 'uncontrolled', 'external'))
    stop("'design' must be 'controlled', 'uncontrolled', or 'external'")

  if (!is.character(prior) || length(prior) != 1L ||
      !prior %in% c('vague', 'N-Inv-Wishart'))
    stop("'prior' must be 'vague' or 'N-Inv-Wishart'")

  # Threshold validation
  if (prob == 'posterior') {
    for (nm in c('theta.TV1', 'theta.MAV1', 'theta.TV2', 'theta.MAV2')) {
      val <- get(nm)
      if (is.null(val) || !is.numeric(val) || length(val) != 1L || is.na(val))
        stop(paste0("'", nm, "' must be a single numeric value when prob = 'posterior'"))
    }
    if (theta.TV1 <= theta.MAV1)
      stop("'theta.TV1' must be strictly greater than 'theta.MAV1'")
    if (theta.TV2 <= theta.MAV2)
      stop("'theta.TV2' must be strictly greater than 'theta.MAV2'")
  } else {
    for (nm in c('theta.NULL1', 'theta.NULL2')) {
      val <- get(nm)
      if (is.null(val) || !is.numeric(val) || length(val) != 1L || is.na(val))
        stop(paste0("'", nm, "' must be a single numeric value when prob = 'predictive'"))
    }
    for (nm in c('m1', 'm2')) {
      val <- get(nm)
      if (is.null(val) || !is.numeric(val) || length(val) != 1L || is.na(val) ||
          val != floor(val) || val < 1L)
        stop(paste0("'", nm, "' must be a single positive integer when prob = 'predictive'"))
    }
    m1 <- as.integer(m1); m2 <- as.integer(m2)
  }

  # n1 validation
  if (!is.numeric(n1) || length(n1) != 1L || is.na(n1) ||
      n1 != floor(n1) || n1 < 1L)
    stop("'n1' must be a single positive integer")
  n1 <- as.integer(n1)
  if (prior == 'vague' && n1 < 3L)
    stop("'n1' must be >= 3 when prior = 'vague' (requires df = n1 - 2 > 0)")

  # ybar1 and S1 validation
  if (!is.numeric(ybar1) || length(ybar1) != 2L || any(is.na(ybar1)))
    stop("'ybar1' must be a numeric vector of length 2")
  if (!is.matrix(S1) || !is.numeric(S1) || nrow(S1) != 2L || ncol(S1) != 2L ||
      any(is.na(S1)))
    stop("'S1' must be a 2x2 numeric matrix")

  # n2, ybar2, S2 for controlled/external
  if (design %in% c('controlled', 'external')) {
    if (is.null(n2) || !is.numeric(n2) || length(n2) != 1L || is.na(n2) ||
        n2 != floor(n2) || n2 < 1L)
      stop("'n2' must be a single positive integer for controlled/external design")
    n2 <- as.integer(n2)
    if (prior == 'vague' && n2 < 3L)
      stop("'n2' must be >= 3 when prior = 'vague'")
    if (is.null(ybar2) || !is.numeric(ybar2) || length(ybar2) != 2L ||
        any(is.na(ybar2)))
      stop("'ybar2' must be a numeric vector of length 2 for controlled/external design")
    if (is.null(S2) || !is.matrix(S2) || !is.numeric(S2) ||
        nrow(S2) != 2L || ncol(S2) != 2L || any(is.na(S2)))
      stop("'S2' must be a 2x2 numeric matrix for controlled/external design")
  }

  # NIW prior parameter validation
  if (prior == 'N-Inv-Wishart') {
    for (nm in c('kappa01', 'nu01')) {
      val <- get(nm)
      if (is.null(val) || !is.numeric(val) || length(val) != 1L || is.na(val) ||
          val <= 0)
        stop(paste0("'", nm, "' must be a single positive numeric when prior = 'N-Inv-Wishart'"))
    }
    if (nu01 <= 2)
      stop("'nu01' must be > 2 for a proper NIW prior")
    if (is.null(mu01) || !is.numeric(mu01) || length(mu01) != 2L ||
        any(is.na(mu01)))
      stop("'mu01' must be a numeric vector of length 2 when prior = 'N-Inv-Wishart'")
    if (is.null(Lambda01) || !is.matrix(Lambda01) || !is.numeric(Lambda01) ||
        nrow(Lambda01) != 2L || ncol(Lambda01) != 2L || any(is.na(Lambda01)))
      stop("'Lambda01' must be a 2x2 numeric matrix when prior = 'N-Inv-Wishart'")

    if (design %in% c('controlled', 'external')) {
      for (nm in c('kappa02', 'nu02')) {
        val <- get(nm)
        if (is.null(val) || !is.numeric(val) || length(val) != 1L || is.na(val) ||
            val <= 0)
          stop(paste0("'", nm, "' must be a single positive numeric when prior = 'N-Inv-Wishart'"))
      }
      if (nu02 <= 2)
        stop("'nu02' must be > 2 for a proper NIW prior")
      if (is.null(mu02) || !is.numeric(mu02) || length(mu02) != 2L ||
          any(is.na(mu02)))
        stop("'mu02' must be a numeric vector of length 2 when prior = 'N-Inv-Wishart'")
      if (is.null(Lambda02) || !is.matrix(Lambda02) || !is.numeric(Lambda02) ||
          nrow(Lambda02) != 2L || ncol(Lambda02) != 2L || any(is.na(Lambda02)))
        stop("'Lambda02' must be a 2x2 numeric matrix when prior = 'N-Inv-Wishart'")
    }
  }

  # Uncontrolled design: mu02 and r
  if (design == 'uncontrolled') {
    if (is.null(mu02) || !is.numeric(mu02) || length(mu02) != 2L ||
        any(is.na(mu02)))
      stop("'mu02' must be a numeric vector of length 2 when design = 'uncontrolled'")
    if (is.null(r) || !is.numeric(r) || length(r) != 1L || is.na(r) || r <= 0)
      stop("'r' must be a single positive numeric when design = 'uncontrolled'")
  }

  # External design: power prior parameters
  if (design == 'external') {
    has_ext1 <- !is.null(ne1)
    has_ext2 <- !is.null(ne2)
    if (!has_ext1 && !has_ext2)
      stop("For design = 'external', at least one of 'ne1' or 'ne2' must be non-NULL")

    if (has_ext1) {
      if (is.null(alpha01e) || is.null(ybar_e1) || is.null(Se1))
        stop("'alpha01e', 'ybar_e1', and 'Se1' must be non-NULL when 'ne1' is specified")
      if (!is.numeric(ne1) || length(ne1) != 1L || ne1 != floor(ne1) || ne1 < 1L)
        stop("'ne1' must be a single positive integer")
      if (!is.numeric(alpha01e) || length(alpha01e) != 1L ||
          alpha01e <= 0 || alpha01e > 1)
        stop("'alpha01e' must be a single numeric in (0, 1]")
      if (!is.numeric(ybar_e1) || length(ybar_e1) != 2L || any(is.na(ybar_e1)))
        stop("'ybar_e1' must be a numeric vector of length 2")
      if (!is.matrix(Se1) || !is.numeric(Se1) || nrow(Se1) != 2L ||
          ncol(Se1) != 2L || any(is.na(Se1)))
        stop("'Se1' must be a 2x2 numeric matrix")
      ne1 <- as.integer(ne1)
    }

    if (has_ext2) {
      if (is.null(alpha02e) || is.null(ybar_e2) || is.null(Se2))
        stop("'alpha02e', 'ybar_e2', and 'Se2' must be non-NULL when 'ne2' is specified")
      if (!is.numeric(ne2) || length(ne2) != 1L || ne2 != floor(ne2) || ne2 < 1L)
        stop("'ne2' must be a single positive integer")
      if (!is.numeric(alpha02e) || length(alpha02e) != 1L ||
          alpha02e <= 0 || alpha02e > 1)
        stop("'alpha02e' must be a single numeric in (0, 1]")
      if (!is.numeric(ybar_e2) || length(ybar_e2) != 2L || any(is.na(ybar_e2)))
        stop("'ybar_e2' must be a numeric vector of length 2")
      if (!is.matrix(Se2) || !is.numeric(Se2) || nrow(Se2) != 2L ||
          ncol(Se2) != 2L || any(is.na(Se2)))
        stop("'Se2' must be a 2x2 numeric matrix")
      ne2 <- as.integer(ne2)
    }
  }

  if (!is.numeric(nMC) || length(nMC) != 1L || is.na(nMC) ||
      nMC != floor(nMC) || nMC < 1L)
    stop("'nMC' must be a single positive integer")
  nMC <- as.integer(nMC)

  if (!is.character(method) || length(method) != 1L ||
      !method %in% c('MC', 'MM'))
    stop("'method' must be either 'MC' or 'MM'")

  # ---------------------------------------------------------------------------
  # Section 2: Compute posterior hyperparameters for each arm
  # ---------------------------------------------------------------------------

  use_pred <- (prob == 'predictive')

  # NIW posterior update: returns list(df, mu_n, V)
  # V = posterior scale matrix            when for_predictive = FALSE
  # V = single-observation predictive scale matrix  when for_predictive = TRUE
  .niw_post <- function(n_k, ybar_k, S_k, kappa0, nu0, mu0, Lambda0,
                        for_predictive = FALSE) {
    kappa_n  <- kappa0 + n_k
    nu_n     <- nu0 + n_k
    mu_n     <- (kappa0 * mu0 + n_k * ybar_k) / kappa_n
    Lambda_n <- Lambda0 + S_k +
      (kappa0 * n_k / kappa_n) * tcrossprod(ybar_k - mu0)
    df <- nu_n - 1
    V  <- if (for_predictive) {
      ((1 + kappa_n) / kappa_n) * Lambda_n / (nu_n - 1)
    } else {
      Lambda_n / (kappa_n * (nu_n - 1))
    }
    list(df = df, mu_n = mu_n, V = V)
  }

  # Vague posterior update
  .vague_post <- function(n_k, ybar_k, S_k, for_predictive = FALSE) {
    df <- n_k - 2L
    V  <- if (for_predictive) {
      ((1 + n_k) / n_k) * S_k / (n_k - 2)
    } else {
      S_k / (n_k * (n_k - 2))
    }
    list(df = df, mu_n = ybar_k, V = V)
  }

  # NIW power prior then PoC update (Conjugate_power_prior.pdf, Theorem 5)
  .niw_ext_post <- function(n_k, ybar_k, S_k, kappa0, nu0, mu0, Lambda0,
                            n_e, alpha0e, ybar_e, S_e,
                            for_predictive = FALSE) {
    # Power prior: effective NIW hyperparameters from external data
    kappa_e  <- alpha0e * n_e + kappa0
    nu_e     <- alpha0e * n_e + nu0
    mu_e     <- (alpha0e * n_e * ybar_e + kappa0 * mu0) / kappa_e
    Lambda_e <- alpha0e * S_e + Lambda0 +
      (kappa0 * alpha0e * n_e / kappa_e) * tcrossprod(mu0 - ybar_e)
    # Update with current PoC data
    .niw_post(n_k, ybar_k, S_k,
              kappa0 = kappa_e, nu0 = nu_e,
              mu0 = mu_e, Lambda0 = Lambda_e,
              for_predictive = for_predictive)
  }

  # Arm 1 (treatment)
  if (design == 'external' && !is.null(ne1)) {
    post1 <- .niw_ext_post(n1, ybar1, S1,
                           kappa0 = kappa01, nu0 = nu01,
                           mu0 = mu01, Lambda0 = Lambda01,
                           n_e = ne1, alpha0e = alpha01e,
                           ybar_e = ybar_e1, S_e = Se1,
                           for_predictive = use_pred)
  } else if (prior == 'N-Inv-Wishart') {
    post1 <- .niw_post(n1, ybar1, S1,
                       kappa0 = kappa01, nu0 = nu01,
                       mu0 = mu01, Lambda0 = Lambda01,
                       for_predictive = use_pred)
  } else {
    post1 <- .vague_post(n1, ybar1, S1, for_predictive = use_pred)
  }

  # Arm 2 (control or hypothetical)
  if (design == 'uncontrolled') {
    # Hypothetical control: same df as treatment arm, scale = r * treatment scale
    # (V already includes predictive inflation if use_pred = TRUE)
    post2 <- list(df = post1$df, mu_n = mu02, V = r * post1$V)
  } else if (design == 'external' && !is.null(ne2)) {
    post2 <- .niw_ext_post(n2, ybar2, S2,
                           kappa0 = kappa02, nu0 = nu02,
                           mu0 = mu02, Lambda0 = Lambda02,
                           n_e = ne2, alpha0e = alpha02e,
                           ybar_e = ybar_e2, S_e = Se2,
                           for_predictive = use_pred)
  } else if (prior == 'N-Inv-Wishart') {
    post2 <- .niw_post(n2, ybar2, S2,
                       kappa0 = kappa02, nu0 = nu02,
                       mu0 = mu02, Lambda0 = Lambda02,
                       for_predictive = use_pred)
  } else {
    post2 <- .vague_post(n2, ybar2, S2, for_predictive = use_pred)
  }

  # ---------------------------------------------------------------------------
  # Section 3: Region probability computation (MC or MM)
  # ---------------------------------------------------------------------------

  # --- MM method: Theorem 3 moment-matching helper ---
  # Computes the parameters of the approximating bivariate t-distribution for
  # the difference Z = mu_t - mu_c, where mu_k ~ t_p(mu_nk, V_k, df_k).
  # Returns list(mu_diff, Sigma_star, nu_star).
  .mm_general <- function(mu1, V1, nu1, mu2, V2, nu2) {
    beta1  <- nu1 / (nu1 - 2)
    beta2  <- nu2 / (nu2 - 2)
    alpha1 <- nu1 ^ 2 / ((nu1 - 2) * (nu1 - 4))
    alpha2 <- nu2 ^ 2 / ((nu2 - 2) * (nu2 - 4))
    gamma  <- 2 * nu1 * nu2 / ((nu1 - 2) * (nu2 - 2))

    # Variance-covariance sum of the two arms
    Vsum <- V1 * beta1 + V2 * beta2

    # A = Vsum^{-1} via Cholesky for numerical stability
    R_chol <- tryCatch(chol(Vsum), error = function(e) NULL)
    A <- if (!is.null(R_chol)) chol2inv(R_chol) else solve(Vsum)

    # Trace quantities required for Qm (Theorem 3)
    AV1      <- A %*% V1
    AV2      <- A %*% V2
    tr_AV1   <- sum(diag(AV1))
    tr_AV2   <- sum(diag(AV2))
    tr_AV1V1 <- sum(diag(AV1 %*% AV1))
    tr_AV2V2 <- sum(diag(AV2 %*% AV2))
    tr_AV1V2 <- sum(diag(AV1 %*% AV2))

    p  <- length(mu1)
    Qm <- (alpha1 * (tr_AV1 ^ 2 + 2 * tr_AV1V1) +
             alpha2 * (tr_AV2 ^ 2 + 2 * tr_AV2V2) +
             gamma  * (tr_AV1 * tr_AV2 + 2 * tr_AV1V2)) / (p * (p + 2))

    nu_star    <- (2 - 4 * Qm) / (1 - Qm)
    Sigma_star <- Vsum * (nu_star - 2) / nu_star

    list(mu_diff = mu1 - mu2, Sigma_star = Sigma_star, nu_star = nu_star)
  }

  # --- Helper: rectangular region probability via inclusion-exclusion ---
  # Computes P(lo1 < Z1 <= hi1, lo2 < Z2 <= hi2) for Z ~ t(mu_diff, Sigma*, nu*)
  # using mvtnorm::pmvt on the centred variable Z - mu_diff.
  # -Inf / +Inf boundaries are handled naturally by pmvt.
  .rect_prob <- function(lo1, hi1, lo2, hi2, mm) {
    nu_r   <- round(mm$nu_star)
    sig    <- mm$Sigma_star
    mu_d   <- mm$mu_diff
    # pmvt integrates P(lower < Z - mu_diff < upper)
    corners <- rbind(
      c(lo1 - mu_d[1L], lo2 - mu_d[2L]),   # lower-lower  (+)
      c(hi1 - mu_d[1L], lo2 - mu_d[2L]),   # upper-lower  (-)
      c(lo1 - mu_d[1L], hi2 - mu_d[2L]),   # lower-upper  (-)
      c(hi1 - mu_d[1L], hi2 - mu_d[2L])    # upper-upper  (+)
    )
    signs <- c(1, -1, -1, 1)
    val   <- 0
    for (k in seq_len(4L)) {
      lo_k <- c(-Inf, -Inf)
      hi_k <- corners[k, ]
      val  <- val + signs[k] *
        as.numeric(mvtnorm::pmvt(lower = lo_k, upper = hi_k,
                                 delta = c(0, 0),
                                 sigma = sig, df = nu_r))
    }
    pmax(val, 0)   # guard against small negative numerical noise
  }

  # --- Determine whether MM is feasible (both df > 4) ---
  use_mm <- (method == 'MM')
  if (use_mm && (post1$df <= 4 || post2$df <= 4)) {
    warning(
      "MM method requires df > 4 for both arms (df1 = ", post1$df,
      ", df2 = ", post2$df, "). Falling back to method = 'MC'."
    )
    use_mm <- FALSE
  }

  if (!use_mm) {

    # -------------------------------------------------------------------------
    # MC path (original implementation)
    # -------------------------------------------------------------------------

    if (prob == 'posterior') {

      # Draw nMC samples from each arm's marginal posterior distribution
      mu_t <- rnsbt(nMC, df = post1$df, mu = post1$mu_n, V = post1$V)
      mu_c <- rnsbt(nMC, df = post2$df, mu = post2$mu_n, V = post2$V)

      theta1 <- mu_t[, 1L] - mu_c[, 1L]
      theta2 <- mu_t[, 2L] - mu_c[, 2L]

      # 3x3 region grid, row-major, Endpoint 1 slow index
      r1 <- 3L - as.integer(theta1 > theta.MAV1) - as.integer(theta1 > theta.TV1)
      r2 <- 3L - as.integer(theta2 > theta.MAV2) - as.integer(theta2 > theta.TV2)
      region <- (r1 - 1L) * 3L + r2

      Pr_R        <- tabulate(region, nbins = 9L) / nMC
      names(Pr_R) <- paste0("R", 1:9)

    } else {

      # V already inflated for one future observation; divide by m_k for the
      # distribution of the mean of m_k future observations
      V1_pred <- post1$V / m1
      V2_pred <- post2$V / m2

      ytilde_t <- rnsbt(nMC, df = post1$df, mu = post1$mu_n, V = V1_pred)
      ytilde_c <- rnsbt(nMC, df = post2$df, mu = post2$mu_n, V = V2_pred)

      theta1 <- ytilde_t[, 1L] - ytilde_c[, 1L]
      theta2 <- ytilde_t[, 2L] - ytilde_c[, 2L]

      # 2x2 region grid, row-major, Endpoint 1 slow index
      r1 <- 2L - as.integer(theta1 > theta.NULL1)
      r2 <- 2L - as.integer(theta2 > theta.NULL2)
      region <- (r1 - 1L) * 2L + r2

      Pr_R        <- tabulate(region, nbins = 4L) / nMC
      names(Pr_R) <- paste0("R", 1:4)
    }

  } else {

    # -------------------------------------------------------------------------
    # MM path: Theorem 3 moment-matching + mvtnorm::pmvt
    # -------------------------------------------------------------------------

    if (prob == 'posterior') {

      mm <- .mm_general(post1$mu_n, post1$V, post1$df,
                        post2$mu_n, post2$V, post2$df)

      # Endpoint 1 cut-points: -Inf < MAV1 < TV1 < +Inf  (3 bands)
      # Endpoint 2 cut-points: -Inf < MAV2 < TV2 < +Inf  (3 bands)
      # 3x3 = 9 regions, row-major, Endpoint 1 slow index
      cuts1 <- c(-Inf, theta.MAV1, theta.TV1,  Inf)
      cuts2 <- c(-Inf, theta.MAV2, theta.TV2,  Inf)

      Pr_R <- numeric(9L)
      idx  <- 1L
      # Endpoint 1 band: from high to low to match row-major order R1..R9
      # Row 1: theta1 > TV1  (cuts1[3] to cuts1[4])
      # Row 2: MAV1 < theta1 <= TV1  (cuts1[2] to cuts1[3])
      # Row 3: theta1 <= MAV1  (cuts1[1] to cuts1[2])
      for (i in 3L:1L) {
        lo1 <- cuts1[i];  hi1 <- cuts1[i + 1L]
        for (j in 3L:1L) {
          lo2 <- cuts2[j];  hi2 <- cuts2[j + 1L]
          Pr_R[idx] <- .rect_prob(lo1, hi1, lo2, hi2, mm)
          idx <- idx + 1L
        }
      }
      # Normalise to ensure exact sum-to-1 (corrects for pmvt rounding)
      Pr_R        <- Pr_R / sum(Pr_R)
      names(Pr_R) <- paste0("R", 1:9)

    } else {

      # Predictive: scale matrices already inflated for single obs; divide by m_k
      V1_pred <- post1$V / m1
      V2_pred <- post2$V / m2

      mm <- .mm_general(post1$mu_n, V1_pred, post1$df,
                        post2$mu_n, V2_pred, post2$df)

      # 2x2 regions (row-major, Endpoint 1 slow index):
      # R1: theta1 > NULL1  AND  theta2 > NULL2
      # R2: theta1 > NULL1  AND  theta2 <= NULL2
      # R3: theta1 <= NULL1 AND  theta2 > NULL2
      # R4: theta1 <= NULL1 AND  theta2 <= NULL2
      cuts1 <- c(-Inf, theta.NULL1, Inf)
      cuts2 <- c(-Inf, theta.NULL2, Inf)

      Pr_R <- numeric(4L)
      idx  <- 1L
      for (i in 2L:1L) {
        lo1 <- cuts1[i];  hi1 <- cuts1[i + 1L]
        for (j in 2L:1L) {
          lo2 <- cuts2[j];  hi2 <- cuts2[j + 1L]
          Pr_R[idx] <- .rect_prob(lo1, hi1, lo2, hi2, mm)
          idx <- idx + 1L
        }
      }
      Pr_R        <- Pr_R / sum(Pr_R)
      names(Pr_R) <- paste0("R", 1:4)
    }
  }

  Pr_R
}
