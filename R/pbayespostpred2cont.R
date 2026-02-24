#' Region Probabilities for Two Continuous Endpoints
#'
#' Computes the posterior or posterior predictive probability that the
#' bivariate treatment effect \eqn{\theta = \mu_t - \mu_c} falls in each of
#' the 9 (posterior) or 4 (predictive) rectangular regions defined by the
#' decision thresholds.
#'
#' @param prob Character scalar: \code{'posterior'} or \code{'predictive'}.
#' @param design Character scalar: \code{'controlled'}, \code{'uncontrolled'},
#'        or \code{'external'}.
#' @param prior Character scalar: \code{'vague'} or \code{'N-Inv-Wishart'}.
#' @param theta.TV1 Numeric scalar. Target value for Endpoint 1
#'        (required when \code{prob = 'posterior'}).
#' @param theta.MAV1 Numeric scalar. Minimum acceptable value for Endpoint 1
#'        (required when \code{prob = 'posterior'}).
#' @param theta.TV2 Numeric scalar. Target value for Endpoint 2
#'        (required when \code{prob = 'posterior'}).
#' @param theta.MAV2 Numeric scalar. Minimum acceptable value for Endpoint 2
#'        (required when \code{prob = 'posterior'}).
#' @param theta.NULL1 Numeric scalar. Null threshold for Endpoint 1
#'        (required when \code{prob = 'predictive'}).
#' @param theta.NULL2 Numeric scalar. Null threshold for Endpoint 2
#'        (required when \code{prob = 'predictive'}).
#' @param n1 Positive integer. Sample size for the treatment arm.
#' @param n2 Positive integer or \code{NULL}. Sample size for the control arm
#'        (not used when \code{design = 'uncontrolled'}).
#' @param ybar1 Numeric vector of length 2 \strong{or} a numeric matrix with
#'        2 columns. Sample mean(s) for the treatment arm. When a matrix with
#'        \eqn{N} rows is supplied, the function computes probabilities for
#'        all \eqn{N} observations simultaneously and returns an
#'        \eqn{N \times n_{\rm regions}} matrix.
#' @param S1 A 2x2 numeric matrix (single observation) \strong{or} a list of
#'        \eqn{N} such matrices. Sum-of-squares matrix/matrices for the
#'        treatment arm. Must be consistent with \code{ybar1}: if
#'        \code{ybar1} is a matrix, \code{S1} must be a list of the same
#'        length.
#' @param ybar2 Numeric vector of length 2, a numeric matrix with 2 columns,
#'        or \code{NULL}. Sample mean(s) for the control arm.
#' @param S2 A 2x2 numeric matrix, a list of 2x2 matrices, or \code{NULL}.
#'        Sum-of-squares matrix/matrices for the control arm.
#' @param m1 Positive integer or \code{NULL}. Pivotal study sample size for
#'        the treatment arm (required when \code{prob = 'predictive'}).
#' @param m2 Positive integer or \code{NULL}. Pivotal study sample size for
#'        the control arm (required when \code{prob = 'predictive'} and
#'        \code{design != 'uncontrolled'}).
#' @param kappa01 Positive numeric scalar. NIW prior concentration for the
#'        treatment arm (required when \code{prior = 'N-Inv-Wishart'}).
#' @param nu01 Numeric scalar \eqn{> 3}. NIW prior degrees of freedom for
#'        the treatment arm (required when \code{prior = 'N-Inv-Wishart'}).
#' @param mu01 Length-2 numeric vector. NIW prior mean for the treatment arm
#'        (required when \code{prior = 'N-Inv-Wishart'}).
#' @param Lambda01 A 2x2 positive-definite numeric matrix. NIW prior scale
#'        matrix for the treatment arm (required when
#'        \code{prior = 'N-Inv-Wishart'}).
#' @param kappa02 Positive numeric scalar. NIW prior concentration for the
#'        control arm.
#' @param nu02 Numeric scalar \eqn{> 3}. NIW prior degrees of freedom for
#'        the control arm.
#' @param mu02 Length-2 numeric vector. NIW prior mean for the control arm or
#'        hypothetical control location (\code{design = 'uncontrolled'}).
#' @param Lambda02 A 2x2 positive-definite numeric matrix. NIW prior scale
#'        matrix for the control arm.
#' @param r Positive numeric scalar. Variance scaling factor for the
#'        hypothetical control (required when
#'        \code{design = 'uncontrolled'}).
#' @param ne1 Positive integer or \code{NULL}. External treatment sample size
#'        (used when \code{design = 'external'}).
#' @param ne2 Positive integer or \code{NULL}. External control sample size
#'        (used when \code{design = 'external'}).
#' @param alpha01e Numeric scalar in \code{(0, 1]}. Power prior weight for
#'        external treatment data.
#' @param alpha02e Numeric scalar in \code{(0, 1]}. Power prior weight for
#'        external control data.
#' @param ybar_e1 Length-2 numeric vector. External treatment sample mean.
#' @param ybar_e2 Length-2 numeric vector. External control sample mean.
#' @param Se1 A 2x2 numeric matrix. External treatment sum-of-squares matrix.
#' @param Se2 A 2x2 numeric matrix. External control sum-of-squares matrix.
#' @param nMC Positive integer or \code{NULL}. Number of Monte Carlo draws.
#'        Default \code{10000L}. Required when \code{method = 'MC'}.  May be
#'        set to \code{NULL} when \code{method = 'MM'} and \eqn{\nu_k > 4}
#'        (the MM method uses \code{mvtnorm::pmvt} analytically); if
#'        \code{method = 'MM'} but \eqn{\nu_k \le 4} causes a fallback to MC,
#'        \code{nMC} must be a positive integer.
#' @param method Character scalar: \code{'MC'} (default) or \code{'MM'}
#'        (Moment-Matching via \code{mvtnorm::pmvt}). When
#'        \code{method = 'MM'} and \eqn{\nu_k \le 4}, a warning is issued
#'        and the function falls back to \code{method = 'MC'}.
#'
#' @return When \code{ybar1} is a length-2 vector (single observation): a
#'         named numeric vector of length 9 (\code{R1}--\code{R9}) for
#'         \code{prob = 'posterior'} or length 4 (\code{R1}--\code{R4}) for
#'         \code{prob = 'predictive'}. All elements are non-negative and sum
#'         to 1.
#'
#'         When \code{ybar1} is an \eqn{N \times 2} matrix (vectorised call):
#'         a numeric matrix with \eqn{N} rows and 9 (or 4) columns, with
#'         column names \code{R1}--\code{R9} (or \code{R1}--\code{R4}). Each
#'         row sums to 1.
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
#'   \item R1: \eqn{\tilde\theta_1 > \theta_{\rm NULL1}} AND
#'             \eqn{\tilde\theta_2 > \theta_{\rm NULL2}}
#'   \item R2: \eqn{\tilde\theta_1 > \theta_{\rm NULL1}} AND
#'             \eqn{\tilde\theta_2 \le \theta_{\rm NULL2}}
#'   \item R3: \eqn{\tilde\theta_1 \le \theta_{\rm NULL1}} AND
#'             \eqn{\tilde\theta_2 > \theta_{\rm NULL2}}
#'   \item R4: \eqn{\tilde\theta_1 \le \theta_{\rm NULL1}} AND
#'             \eqn{\tilde\theta_2 \le \theta_{\rm NULL2}}
#' }
#'
#' \strong{Uncontrolled design.}
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
#' \strong{Vectorised usage.}
#' When \code{ybar1} is supplied as an \eqn{N \times 2} matrix and
#' \code{S1} as a list of \eqn{N} scatter matrices, the function computes
#' region probabilities for all \eqn{N} observations in a single call,
#' returning an \eqn{N \times n_{\rm regions}} matrix.  For
#' \code{method = 'MC'}, standard normal and chi-squared variates are
#' pre-generated once (size \code{nMC}) and reused across all observations,
#' with only the Cholesky factor of the replicate-specific scale matrix
#' recomputed per observation.  For \code{method = 'MM'}, the
#' moment-matching parameters are computed per observation (since they depend
#' on the replicate-specific \eqn{V_k}) and \code{mvtnorm::pmvt} is called
#' once per region per observation.
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
#' # Example 7: Vectorised call -- N = 3 observations simultaneously
#' S1_list <- list(
#'   matrix(c(18.0, 3.6, 3.6, 9.0), 2, 2),
#'   matrix(c(20.0, 4.0, 4.0, 10.0), 2, 2),
#'   matrix(c(15.0, 2.5, 2.5, 7.5), 2, 2)
#' )
#' S2_list <- list(
#'   matrix(c(16.0, 2.8, 2.8, 8.5), 2, 2),
#'   matrix(c(17.0, 3.0, 3.0, 9.0), 2, 2),
#'   matrix(c(14.0, 2.2, 2.2, 7.0), 2, 2)
#' )
#' ybar1_mat <- rbind(c(3.5, 2.1), c(4.0, 2.5), c(2.5, 1.5))
#' ybar2_mat <- rbind(c(1.8, 1.0), c(1.9, 1.1), c(1.7, 0.9))
#' pbayespostpred2cont(
#'   prob = 'posterior', design = 'controlled', prior = 'vague',
#'   theta.TV1 = 1.5, theta.MAV1 = 0.5,
#'   theta.TV2 = 1.0, theta.MAV2 = 0.3,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   n1 = 12L, n2 = 12L,
#'   ybar1 = ybar1_mat, S1 = S1_list,
#'   ybar2 = ybar2_mat, S2 = S2_list,
#'   m1 = NULL, m2 = NULL,
#'   kappa01 = NULL, nu01 = NULL, mu01 = NULL, Lambda01 = NULL,
#'   kappa02 = NULL, nu02 = NULL, mu02 = NULL, Lambda02 = NULL,
#'   r = NULL,
#'   ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
#'   ybar_e1 = NULL, ybar_e2 = NULL, Se1 = NULL, Se2 = NULL,
#'   nMC = 1000L
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

  if (prob == 'posterior') {
    for (nm in c('theta.TV1', 'theta.MAV1', 'theta.TV2', 'theta.MAV2')) {
      val <- get(nm)
      if (is.null(val) || !is.numeric(val) || length(val) != 1L || is.na(val))
        stop(paste0("'", nm, "' must be a single numeric when prob = 'posterior'"))
    }
    if (theta.TV1 <= theta.MAV1)
      stop("'theta.TV1' must be strictly greater than 'theta.MAV1'")
    if (theta.TV2 <= theta.MAV2)
      stop("'theta.TV2' must be strictly greater than 'theta.MAV2'")
  } else {
    for (nm in c('theta.NULL1', 'theta.NULL2')) {
      val <- get(nm)
      if (is.null(val) || !is.numeric(val) || length(val) != 1L || is.na(val))
        stop(paste0("'", nm, "' must be a single numeric when prob = 'predictive'"))
    }
    for (nm in c('m1', 'm2')) {
      val <- get(nm)
      if (is.null(val) || !is.numeric(val) || length(val) != 1L || is.na(val) ||
          val != floor(val) || val < 1L)
        stop(paste0("'", nm, "' must be a positive integer when prob = 'predictive'"))
    }
    m1 <- as.integer(m1); m2 <- as.integer(m2)
  }

  if (!is.numeric(n1) || length(n1) != 1L || is.na(n1) ||
      n1 != floor(n1) || n1 < 2L)
    stop("'n1' must be a single integer >= 2")
  n1 <- as.integer(n1)

  if (design %in% c('controlled', 'external')) {
    if (is.null(n2) || !is.numeric(n2) || length(n2) != 1L || is.na(n2) ||
        n2 != floor(n2) || n2 < 2L)
      stop("'n2' must be a single integer >= 2 for controlled/external design")
    n2 <- as.integer(n2)
  }

  # Determine whether vectorised (N > 1) or scalar call
  vec_call <- is.matrix(ybar1) && nrow(ybar1) > 1L
  if (vec_call) {
    if (ncol(ybar1) != 2L)
      stop("When 'ybar1' is a matrix it must have exactly 2 columns")
    N <- nrow(ybar1)
    if (!is.list(S1) || length(S1) != N)
      stop("When 'ybar1' is a matrix, 'S1' must be a list of the same length")
    if (design %in% c('controlled', 'external')) {
      if (is.null(ybar2) || !is.matrix(ybar2) || nrow(ybar2) != N ||
          ncol(ybar2) != 2L)
        stop("When 'ybar1' is a matrix, 'ybar2' must be a matrix with the same number of rows")
      if (!is.list(S2) || length(S2) != N)
        stop("When 'ybar1' is a matrix, 'S2' must be a list of the same length")
    }
  } else {
    N <- 1L
    # Coerce scalar ybar1 to length-2 vector
    if (!is.numeric(ybar1) || length(ybar1) != 2L || any(is.na(ybar1)))
      stop("'ybar1' must be a numeric vector of length 2 or a matrix with 2 columns")
    if (!is.matrix(S1) || !is.numeric(S1) || nrow(S1) != 2L || ncol(S1) != 2L)
      stop("'S1' must be a 2x2 numeric matrix (or a list for vectorised call)")
    # Wrap as list for uniform handling below
    ybar1 <- matrix(ybar1, nrow = 1L)
    S1    <- list(S1)
    if (design %in% c('controlled', 'external')) {
      if (!is.numeric(ybar2) || length(ybar2) != 2L || any(is.na(ybar2)))
        stop("'ybar2' must be a numeric vector of length 2")
      if (!is.matrix(S2) || !is.numeric(S2) || nrow(S2) != 2L || ncol(S2) != 2L)
        stop("'S2' must be a 2x2 numeric matrix")
      ybar2 <- matrix(ybar2, nrow = 1L)
      S2    <- list(S2)
    }
  }

  if (prior == 'N-Inv-Wishart') {
    for (nm in c('kappa01', 'nu01')) {
      val <- get(nm)
      if (is.null(val) || !is.numeric(val) || length(val) != 1L ||
          is.na(val) || val <= 0)
        stop(paste0("'", nm, "' must be a single positive numeric for NIW prior"))
    }
    if (is.null(mu01) || !is.numeric(mu01) || length(mu01) != 2L)
      stop("'mu01' must be a length-2 numeric vector for NIW prior")
    if (is.null(Lambda01) || !is.matrix(Lambda01) || nrow(Lambda01) != 2L)
      stop("'Lambda01' must be a 2x2 numeric matrix for NIW prior")
    if (design %in% c('controlled', 'external')) {
      for (nm in c('kappa02', 'nu02')) {
        val <- get(nm)
        if (is.null(val) || !is.numeric(val) || length(val) != 1L ||
            is.na(val) || val <= 0)
          stop(paste0("'", nm, "' must be a single positive numeric for NIW prior"))
      }
      if (is.null(mu02) || !is.numeric(mu02) || length(mu02) != 2L)
        stop("'mu02' must be a length-2 numeric vector for NIW prior")
      if (is.null(Lambda02) || !is.matrix(Lambda02) || nrow(Lambda02) != 2L)
        stop("'Lambda02' must be a 2x2 numeric matrix for NIW prior")
    }
  }

  if (design == 'uncontrolled') {
    if (is.null(r) || !is.numeric(r) || length(r) != 1L || is.na(r) || r <= 0)
      stop("'r' must be a single positive numeric when design = 'uncontrolled'")
    if (is.null(mu02) || !is.numeric(mu02) || length(mu02) != 2L)
      stop("'mu02' must be a length-2 numeric vector when design = 'uncontrolled'")
  }

  if (design == 'external') {
    has_ext1 <- !is.null(ne1) && !is.null(alpha01e) &&
      !is.null(ybar_e1) && !is.null(Se1)
    has_ext2 <- !is.null(ne2) && !is.null(alpha02e) &&
      !is.null(ybar_e2) && !is.null(Se2)
    if (!has_ext1 && !has_ext2)
      stop(paste0("For design = 'external', at least one complete set of ",
                  "external data must be provided"))
    if (!is.null(ne1)) ne1 <- as.integer(ne1)
    if (!is.null(ne2)) ne2 <- as.integer(ne2)
  }

  if (!is.character(method) || length(method) != 1L ||
      !method %in% c('MC', 'MM'))
    stop("'method' must be either 'MC' or 'MM'")

  # nMC validation: required for method = 'MC', optional for method = 'MM'
  if (method == 'MC') {
    if (is.null(nMC))
      stop("'nMC' must be non-NULL when method = 'MC'")
    if (!is.numeric(nMC) || length(nMC) != 1L || is.na(nMC) ||
        nMC != floor(nMC) || nMC < 1L)
      stop("'nMC' must be a single positive integer")
    nMC <- as.integer(nMC)
  } else {
    # method == 'MM': nMC may be NULL or a positive integer
    if (!is.null(nMC)) {
      if (!is.numeric(nMC) || length(nMC) != 1L || is.na(nMC) ||
          nMC != floor(nMC) || nMC < 1L)
        stop("'nMC' must be a single positive integer or NULL")
      nMC <- as.integer(nMC)
    }
  }

  # ---------------------------------------------------------------------------
  # Section 2: Posterior hyperparameter helpers
  #
  # Each helper accepts a single (ybar_k, S_k) pair and returns list(df, mu_n, V).
  # V is the posterior scale matrix (prob = 'posterior') or the single-
  # observation predictive scale matrix (prob = 'predictive'); for the
  # predictive mean of m_k observations it is further divided by m_k in
  # Section 3.
  # ---------------------------------------------------------------------------

  use_pred <- (prob == 'predictive')

  # NIW posterior update
  .niw_post <- function(n_k, ybar_k, S_k, kappa0, nu0, mu0, Lambda0,
                        for_predictive = FALSE) {
    kappa_n  <- kappa0 + n_k
    nu_n     <- nu0 + n_k
    mu_n     <- (kappa0 * mu0 + n_k * ybar_k) / kappa_n
    Lambda_n <- Lambda0 + S_k +
      (kappa0 * n_k / kappa_n) * tcrossprod(ybar_k - mu0)
    df <- nu_n - 1L
    V  <- if (for_predictive) {
      ((1 + kappa_n) / kappa_n) * Lambda_n / (nu_n - 1L)
    } else {
      Lambda_n / (kappa_n * (nu_n - 1L))
    }
    list(df = df, mu_n = mu_n, V = V)
  }

  # Vague posterior update
  .vague_post <- function(n_k, ybar_k, S_k, for_predictive = FALSE) {
    df <- n_k - 2L
    V  <- if (for_predictive) {
      ((1 + n_k) / n_k) * S_k / (n_k - 2L)
    } else {
      S_k / (n_k * (n_k - 2L))
    }
    list(df = df, mu_n = ybar_k, V = V)
  }

  # NIW power prior then PoC update (Theorem 5)
  .niw_ext_post <- function(n_k, ybar_k, S_k, kappa0, nu0, mu0, Lambda0,
                            n_e, alpha0e, ybar_e, S_e,
                            for_predictive = FALSE) {
    kappa_e  <- alpha0e * n_e + kappa0
    nu_e     <- alpha0e * n_e + nu0
    mu_e     <- (alpha0e * n_e * ybar_e + kappa0 * mu0) / kappa_e
    Lambda_e <- alpha0e * S_e + Lambda0 +
      (kappa0 * alpha0e * n_e / kappa_e) * tcrossprod(mu0 - ybar_e)
    .niw_post(n_k, ybar_k, S_k,
              kappa0 = kappa_e, nu0 = nu_e,
              mu0 = mu_e, Lambda0 = Lambda_e,
              for_predictive = for_predictive)
  }

  # Pre-compute effective external hyperparameters (scalar, done once)
  if (design == 'external') {
    if (!is.null(ne1)) {
      kappa_e1  <- alpha01e * ne1 + kappa01
      nu_e1     <- alpha01e * ne1 + nu01
      mu_e1     <- (alpha01e * ne1 * ybar_e1 + kappa01 * mu01) / kappa_e1
      Lambda_e1 <- alpha01e * Se1 + Lambda01 +
        (kappa01 * alpha01e * ne1 / kappa_e1) * tcrossprod(mu01 - ybar_e1)
    }
    if (!is.null(ne2)) {
      kappa_e2  <- alpha02e * ne2 + kappa02
      nu_e2     <- alpha02e * ne2 + nu02
      mu_e2     <- (alpha02e * ne2 * ybar_e2 + kappa02 * mu02) / kappa_e2
      Lambda_e2 <- alpha02e * Se2 + Lambda02 +
        (kappa02 * alpha02e * ne2 / kappa_e2) * tcrossprod(mu02 - ybar_e2)
    }
  }

  # ---------------------------------------------------------------------------
  # Section 3: MM helpers (constructed only when method = 'MM')
  # ---------------------------------------------------------------------------

  # MM: moment-matching parameters for the difference of two bivariate t's
  .mm_general <- function(mu1_loc, V1, nu1, mu2_loc, V2, nu2) {
    beta1  <- nu1 / (nu1 - 2)
    beta2  <- nu2 / (nu2 - 2)
    alpha1 <- nu1 ^ 2 / ((nu1 - 2) * (nu1 - 4))
    alpha2 <- nu2 ^ 2 / ((nu2 - 2) * (nu2 - 4))
    gamma  <- 2 * nu1 * nu2 / ((nu1 - 2) * (nu2 - 2))

    Vsum   <- V1 * beta1 + V2 * beta2
    R_chol <- tryCatch(chol(Vsum), error = function(e) NULL)
    A      <- if (!is.null(R_chol)) chol2inv(R_chol) else solve(Vsum)

    AV1      <- A %*% V1
    AV2      <- A %*% V2
    tr_AV1   <- sum(diag(AV1))
    tr_AV2   <- sum(diag(AV2))
    tr_AV1V1 <- sum(diag(AV1 %*% AV1))
    tr_AV2V2 <- sum(diag(AV2 %*% AV2))
    tr_AV1V2 <- sum(diag(AV1 %*% AV2))

    p  <- length(mu1_loc)
    Qm <- (alpha1 * (tr_AV1 ^ 2 + 2 * tr_AV1V1) +
             alpha2 * (tr_AV2 ^ 2 + 2 * tr_AV2V2) +
             gamma  * (tr_AV1 * tr_AV2 + 2 * tr_AV1V2)) / (p * (p + 2))

    nu_star    <- (2 - 4 * Qm) / (1 - Qm)
    Sigma_star <- Vsum * (nu_star - 2) / nu_star

    list(mu_diff = mu1_loc - mu2_loc, Sigma_star = Sigma_star,
         nu_star = nu_star)
  }

  # MM: rectangular region probability via inclusion-exclusion
  .rect_prob <- function(lo1, hi1, lo2, hi2, mm) {
    nu_r <- round(mm$nu_star)
    sig  <- mm$Sigma_star
    mu_d <- mm$mu_diff
    corners <- rbind(
      c(lo1 - mu_d[1L], lo2 - mu_d[2L]),
      c(hi1 - mu_d[1L], lo2 - mu_d[2L]),
      c(lo1 - mu_d[1L], hi2 - mu_d[2L]),
      c(hi1 - mu_d[1L], hi2 - mu_d[2L])
    )
    signs <- c(1, -1, -1, 1)
    val   <- 0
    for (k in seq_len(4L)) {
      val <- val + signs[k] *
        as.numeric(mvtnorm::pmvt(lower = c(-Inf, -Inf), upper = corners[k, ],
                                 delta = c(0, 0), sigma = sig, df = nu_r))
    }
    pmax(val, 0)
  }

  # ---------------------------------------------------------------------------
  # Section 4: Pre-generate raw MC variates (shared across all N observations)
  #
  # df is observation-invariant (depends only on n_k and prior hyperparameters).
  # Z and W are drawn once and reused for every observation in the vectorised
  # call, with only the Cholesky factor of V_k changing per observation.
  # ---------------------------------------------------------------------------

  # Compute df (scalar, same for all observations in a vectorised call)
  if (prior == 'N-Inv-Wishart') {
    df1 <- as.integer(nu01 + n1 - 1L)
    if (design %in% c('controlled', 'external')) df2 <- as.integer(nu02 + n2 - 1L)
    if (design == 'uncontrolled')                df2 <- df1
  } else {
    df1 <- n1 - 2L
    if (design %in% c('controlled', 'external')) df2 <- n2 - 2L
    if (design == 'uncontrolled')                df2 <- df1
  }
  if (design == 'external') {
    if (!is.null(ne1)) df1 <- as.integer(nu_e1 + n1 - 1L)
    if (!is.null(ne2)) df2 <- as.integer(nu_e2 + n2 - 1L)
  }

  # Check MM feasibility
  use_mm <- (method == 'MM')
  if (use_mm && (df1 <= 4L || df2 <= 4L)) {
    warning(
      "MM method requires df > 4 for both arms (df1 = ", df1,
      ", df2 = ", df2, "). Falling back to method = 'MC'."
    )
    use_mm <- FALSE
  }

  # Pre-generate MC raw variates (only for MC path; also used as fallback for MM)
  if (!use_mm) {
    if (is.null(nMC))
      stop(paste0("'nMC' must be a positive integer when method = 'MM' falls ",
                  "back to MC (nu_k <= 4)"))
    Z_t <- matrix(rnorm(nMC * 2L), nrow = nMC, ncol = 2L)
    W_t <- rchisq(nMC, df = df1)
    sc_t <- sqrt(W_t / df1)

    Z_c <- matrix(rnorm(nMC * 2L), nrow = nMC, ncol = 2L)
    W_c <- rchisq(nMC, df = df2)
    sc_c <- sqrt(W_c / df2)
  }

  # ---------------------------------------------------------------------------
  # Section 5: Region probability computation for N observations
  # ---------------------------------------------------------------------------

  n_regions <- if (prob == 'posterior') 9L else 4L
  Pr_R_mat  <- matrix(0, nrow = N, ncol = n_regions)

  # Cut-points (fixed across all observations)
  if (prob == 'posterior') {
    cuts1 <- c(-Inf, theta.MAV1, theta.TV1, Inf)
    cuts2 <- c(-Inf, theta.MAV2, theta.TV2, Inf)
  } else {
    cuts1 <- c(-Inf, theta.NULL1, Inf)
    cuts2 <- c(-Inf, theta.NULL2, Inf)
  }

  for (i in seq_len(N)) {

    yb1 <- ybar1[i, ]
    S1i <- S1[[i]]

    # Posterior parameters for arm 1
    if (design == 'external' && !is.null(ne1)) {
      post1 <- .niw_post(n1, yb1, S1i,
                         kappa0 = kappa_e1, nu0 = nu_e1,
                         mu0 = mu_e1, Lambda0 = Lambda_e1,
                         for_predictive = use_pred)
    } else if (prior == 'N-Inv-Wishart') {
      post1 <- .niw_post(n1, yb1, S1i,
                         kappa0 = kappa01, nu0 = nu01,
                         mu0 = mu01, Lambda0 = Lambda01,
                         for_predictive = use_pred)
    } else {
      post1 <- .vague_post(n1, yb1, S1i, for_predictive = use_pred)
    }

    # Posterior parameters for arm 2
    if (design == 'uncontrolled') {
      post2 <- list(df = post1$df, mu_n = mu02, V = r * post1$V)
    } else if (design == 'external' && !is.null(ne2)) {
      yb2 <- ybar2[i, ]
      S2i <- S2[[i]]
      post2 <- .niw_post(n2, yb2, S2i,
                         kappa0 = kappa_e2, nu0 = nu_e2,
                         mu0 = mu_e2, Lambda0 = Lambda_e2,
                         for_predictive = use_pred)
    } else if (prior == 'N-Inv-Wishart') {
      yb2 <- ybar2[i, ]
      S2i <- S2[[i]]
      post2 <- .niw_post(n2, yb2, S2i,
                         kappa0 = kappa02, nu0 = nu02,
                         mu0 = mu02, Lambda0 = Lambda02,
                         for_predictive = use_pred)
    } else {
      yb2 <- ybar2[i, ]
      S2i <- S2[[i]]
      post2 <- .vague_post(n2, yb2, S2i, for_predictive = use_pred)
    }

    if (!use_mm) {

      # -----------------------------------------------------------------------
      # MC path: re-use pre-generated Z and W; only chol(V) changes per obs
      # -----------------------------------------------------------------------

      V1_use <- if (use_pred) post1$V / m1 else post1$V
      V2_use <- if (use_pred) post2$V / m2 else post2$V

      R1 <- chol(V1_use + diag(1e-10, 2L))
      R2 <- chol(V2_use + diag(1e-10, 2L))

      mu_t <- sweep(Z_t %*% R1, 1L, sc_t, '/') +
        matrix(post1$mu_n, nrow = nMC, ncol = 2L, byrow = TRUE)
      mu_c <- sweep(Z_c %*% R2, 1L, sc_c, '/') +
        matrix(post2$mu_n, nrow = nMC, ncol = 2L, byrow = TRUE)

      theta1 <- mu_t[, 1L] - mu_c[, 1L]
      theta2 <- mu_t[, 2L] - mu_c[, 2L]

      if (prob == 'posterior') {
        r1 <- 3L - as.integer(theta1 > theta.MAV1) -
          as.integer(theta1 > theta.TV1)
        r2 <- 3L - as.integer(theta2 > theta.MAV2) -
          as.integer(theta2 > theta.TV2)
        region <- (r1 - 1L) * 3L + r2
        Pr_R_mat[i, ] <- tabulate(region, nbins = 9L) / nMC
      } else {
        r1 <- 2L - as.integer(theta1 > theta.NULL1)
        r2 <- 2L - as.integer(theta2 > theta.NULL2)
        region <- (r1 - 1L) * 2L + r2
        Pr_R_mat[i, ] <- tabulate(region, nbins = 4L) / nMC
      }

    } else {

      # -----------------------------------------------------------------------
      # MM path: moment-matching + mvtnorm::pmvt per observation
      # -----------------------------------------------------------------------

      V1_use <- if (use_pred) post1$V / m1 else post1$V
      V2_use <- if (use_pred) post2$V / m2 else post2$V

      mm <- .mm_general(post1$mu_n, V1_use, post1$df,
                        post2$mu_n, V2_use, post2$df)

      n_bands <- length(cuts1) - 1L
      Pr_R <- numeric(n_regions)
      idx  <- 1L
      for (ii in n_bands:1L) {
        for (jj in n_bands:1L) {
          Pr_R[idx] <- .rect_prob(cuts1[ii], cuts1[ii + 1L],
                                  cuts2[jj], cuts2[jj + 1L], mm)
          idx <- idx + 1L
        }
      }
      # Normalise to correct for pmvt rounding
      Pr_R_mat[i, ] <- Pr_R / sum(Pr_R)
    }
  }

  # ---------------------------------------------------------------------------
  # Section 6: Return
  # ---------------------------------------------------------------------------

  region_names <- paste0("R", seq_len(n_regions))
  colnames(Pr_R_mat) <- region_names

  if (vec_call) {
    # Vectorised call: return N x n_regions matrix
    Pr_R_mat
  } else {
    # Scalar call: return named vector (original behaviour)
    out <- Pr_R_mat[1L, ]
    names(out) <- region_names
    out
  }
}
