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
#' @param theta_TV1 Numeric scalar. Target value for Endpoint 1
#'        (required when \code{prob = 'posterior'}).
#' @param theta_MAV1 Numeric scalar. Minimum acceptable value for Endpoint 1
#'        (required when \code{prob = 'posterior'}).
#' @param theta_TV2 Numeric scalar. Target value for Endpoint 2
#'        (required when \code{prob = 'posterior'}).
#' @param theta_MAV2 Numeric scalar. Minimum acceptable value for Endpoint 2
#'        (required when \code{prob = 'posterior'}).
#' @param theta_NULL1 Numeric scalar. Null threshold for Endpoint 1
#'        (required when \code{prob = 'predictive'}).
#' @param theta_NULL2 Numeric scalar. Null threshold for Endpoint 2
#'        (required when \code{prob = 'predictive'}).
#' @param n_t Positive integer. Sample size for the treatment group.
#' @param n_c Positive integer or \code{NULL}. Sample size for the control group
#'        (not used when \code{design = 'uncontrolled'}).
#' @param ybar_t Numeric vector of length 2 \strong{or} a numeric matrix with
#'        2 columns. Sample mean(s) for the treatment group. When a matrix with
#'        \eqn{N} rows is supplied, the function computes probabilities for
#'        all \eqn{N} observations simultaneously and returns an
#'        \eqn{N \times n_{\rm regions}} matrix.
#' @param S_t A 2x2 numeric matrix (single observation) \strong{or} a list of
#'        \eqn{N} such matrices. Sum-of-squares matrix/matrices for the
#'        treatment group. Must be consistent with \code{ybar_t}: if
#'        \code{ybar_t} is a matrix, \code{S_t} must be a list of the same
#'        length.
#' @param ybar_c Numeric vector of length 2, a numeric matrix with 2 columns,
#'        or \code{NULL}. Sample mean(s) for the control group.
#' @param S_c A 2x2 numeric matrix, a list of 2x2 matrices, or \code{NULL}.
#'        Sum-of-squares matrix/matrices for the control group.
#' @param m_t Positive integer or \code{NULL}. Pivotal study sample size for
#'        the treatment group (required when \code{prob = 'predictive'}).
#' @param m_c Positive integer or \code{NULL}. Pivotal study sample size for
#'        the control group (required when \code{prob = 'predictive'} and
#'        \code{design != 'uncontrolled'}).
#' @param kappa0_t Positive numeric scalar. NIW prior concentration for the
#'        treatment group (required when \code{prior = 'N-Inv-Wishart'}).
#' @param nu0_t Numeric scalar \eqn{> 3}. NIW prior degrees of freedom for
#'        the treatment group (required when \code{prior = 'N-Inv-Wishart'}).
#' @param mu0_t Length-2 numeric vector. NIW prior mean for the treatment group
#'        (required when \code{prior = 'N-Inv-Wishart'}).
#' @param Lambda0_t A 2x2 positive-definite numeric matrix. NIW prior scale
#'        matrix for the treatment group (required when
#'        \code{prior = 'N-Inv-Wishart'}).
#' @param kappa0_c Positive numeric scalar. NIW prior concentration for the
#'        control group.
#' @param nu0_c Numeric scalar \eqn{> 3}. NIW prior degrees of freedom for
#'        the control group.
#' @param mu0_c Length-2 numeric vector. NIW prior mean for the control group or
#'        hypothetical control location (\code{design = 'uncontrolled'}).
#' @param Lambda0_c A 2x2 positive-definite numeric matrix. NIW prior scale
#'        matrix for the control group.
#' @param r Positive numeric scalar. Variance scaling factor for the
#'        hypothetical control (required when
#'        \code{design = 'uncontrolled'}).
#' @param ne_t Positive integer or \code{NULL}. External treatment sample size
#'        (used when \code{design = 'external'}).
#' @param ne_c Positive integer or \code{NULL}. External control sample size
#'        (used when \code{design = 'external'}).
#' @param alpha0e_t Numeric scalar in \code{(0, 1]}. Power prior weight for
#'        external treatment data.
#' @param alpha0e_c Numeric scalar in \code{(0, 1]}. Power prior weight for
#'        external control data.
#' @param bar_ye_t Length-2 numeric vector. External treatment sample mean.
#' @param bar_ye_c Length-2 numeric vector. External control sample mean.
#' @param se_t A 2x2 numeric matrix. External treatment sum-of-squares matrix.
#' @param se_c A 2x2 numeric matrix. External control sum-of-squares matrix.
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
#' @return When \code{ybar_t} is a length-2 vector (single observation): a
#'         named numeric vector of length 9 (\code{R_chol_t}--\code{R9}) for
#'         \code{prob = 'posterior'} or length 4 (\code{R_chol_t}--\code{R4}) for
#'         \code{prob = 'predictive'}. All elements are non-negative and sum
#'         to 1.
#'
#'         When \code{ybar_t} is an \eqn{N \times 2} matrix (vectorised call):
#'         a numeric matrix with \eqn{N} rows and 9 (or 4) columns, with
#'         column names \code{R_chol_t}--\code{R9} (or \code{R_chol_t}--\code{R4}). Each
#'         row sums to 1.
#'
#' @details
#' \strong{Model.}
#' Both endpoints follow a bivariate Normal distribution
#' \eqn{y_{k,j} \sim N_2(\mu_k, \Sigma_k)} for group \eqn{k \in \{t, c\}}.
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
#'   \item R_chol_t: \eqn{\theta_1 > TV_1} AND \eqn{\theta_2 > TV_2}
#'   \item R_chol_c: \eqn{\theta_1 > TV_1} AND \eqn{TV_2 \ge \theta_2 > MAV_2}
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
#'   \item R_chol_t: \eqn{\tilde\theta_1 > \theta_{\rm NULL1}} AND
#'             \eqn{\tilde\theta_2 > \theta_{\rm NULL2}}
#'   \item R_chol_c: \eqn{\tilde\theta_1 > \theta_{\rm NULL1}} AND
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
#' control to differ from the treatment group.
#'
#' \strong{External design.}
#' Incorporated via the power prior with NIW conjugate representation.
#' The effective posterior hyperparameters are obtained by constructing the
#' power prior from external data with weight \eqn{a_0}, then updating with
#' current PoC data (see Conjugate_power_prior.pdf, Theorem 5).
#'
#' \strong{Vectorised usage.}
#' When \code{ybar_t} is supplied as an \eqn{N \times 2} matrix and
#' \code{S_t} as a list of \eqn{N} scatter matrices, the function computes
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
#' S_t <- matrix(c(18.0, 3.6, 3.6, 9.0), 2, 2)
#' S_c <- matrix(c(16.0, 2.8, 2.8, 8.5), 2, 2)
#' pbayespostpred2cont(
#'   prob = 'posterior', design = 'controlled', prior = 'vague',
#'   theta_TV1 = 1.5, theta_MAV1 = 0.5,
#'   theta_TV2 = 1.0, theta_MAV2 = 0.3,
#'   theta_NULL1 = NULL, theta_NULL2 = NULL,
#'   n_t = 12L, n_c = 12L,
#'   ybar_t = c(3.5, 2.1), S_t = S_t,
#'   ybar_c = c(1.8, 1.0), S_c = S_c,
#'   m_t = NULL, m_c = NULL,
#'   kappa0_t = NULL, nu0_t = NULL, mu0_t = NULL, Lambda0_t = NULL,
#'   kappa0_c = NULL, nu0_c = NULL, mu0_c = NULL, Lambda0_c = NULL,
#'   r = NULL,
#'   ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
#'   bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
#'   nMC = 1000L
#' )
#'
#' # Example 2: Posterior probability, controlled design, NIW prior
#' S_t <- matrix(c(18.0, 3.6, 3.6, 9.0), 2, 2)
#' S_c <- matrix(c(16.0, 2.8, 2.8, 8.5), 2, 2)
#' L0 <- matrix(c(20.0, 0.0, 0.0, 10.0), 2, 2)
#' pbayespostpred2cont(
#'   prob = 'posterior', design = 'controlled', prior = 'N-Inv-Wishart',
#'   theta_TV1 = 1.5, theta_MAV1 = 0.5,
#'   theta_TV2 = 1.0, theta_MAV2 = 0.3,
#'   theta_NULL1 = NULL, theta_NULL2 = NULL,
#'   n_t = 12L, n_c = 12L,
#'   ybar_t = c(3.5, 2.1), S_t = S_t,
#'   ybar_c = c(1.8, 1.0), S_c = S_c,
#'   m_t = NULL, m_c = NULL,
#'   kappa0_t = 2.0, nu0_t = 5.0, mu0_t = c(2.0, 1.0), Lambda0_t = L0,
#'   kappa0_c = 2.0, nu0_c = 5.0, mu0_c = c(1.0, 0.5), Lambda0_c = L0,
#'   r = NULL,
#'   ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
#'   bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
#'   nMC = 1000L
#' )
#'
#' # Example 3: Posterior probability, uncontrolled design, NIW prior
#' # mu0_c specifies the hypothetical control location; r scales the variance
#' S_t <- matrix(c(18.0, 3.6, 3.6, 9.0), 2, 2)
#' L0 <- matrix(c(20.0, 0.0, 0.0, 10.0), 2, 2)
#' pbayespostpred2cont(
#'   prob = 'posterior', design = 'uncontrolled', prior = 'N-Inv-Wishart',
#'   theta_TV1 = 1.5, theta_MAV1 = 0.5,
#'   theta_TV2 = 1.0, theta_MAV2 = 0.3,
#'   theta_NULL1 = NULL, theta_NULL2 = NULL,
#'   n_t = 12L, n_c = NULL,
#'   ybar_t = c(3.5, 2.1), S_t = S_t,
#'   ybar_c = NULL, S_c = NULL,
#'   m_t = NULL, m_c = NULL,
#'   kappa0_t = 2.0, nu0_t = 5.0, mu0_t = c(2.0, 1.0), Lambda0_t = L0,
#'   kappa0_c = NULL, nu0_c = NULL, mu0_c = c(1.0, 0.5), Lambda0_c = NULL,
#'   r = 1.0,
#'   ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
#'   bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
#'   nMC = 1000L
#' )
#'
#' # Example 4: Posterior probability, external design (control only), NIW prior
#' S_t  <- matrix(c(18.0, 3.6, 3.6, 9.0), 2, 2)
#' S_c  <- matrix(c(16.0, 2.8, 2.8, 8.5), 2, 2)
#' L0  <- matrix(c(20.0, 0.0, 0.0, 10.0), 2, 2)
#' se_c <- matrix(c(15.0, 2.5, 2.5, 7.5), 2, 2)
#' pbayespostpred2cont(
#'   prob = 'posterior', design = 'external', prior = 'N-Inv-Wishart',
#'   theta_TV1 = 1.5, theta_MAV1 = 0.5,
#'   theta_TV2 = 1.0, theta_MAV2 = 0.3,
#'   theta_NULL1 = NULL, theta_NULL2 = NULL,
#'   n_t = 12L, n_c = 12L,
#'   ybar_t = c(3.5, 2.1), S_t = S_t,
#'   ybar_c = c(1.8, 1.0), S_c = S_c,
#'   m_t = NULL, m_c = NULL,
#'   kappa0_t = 2.0, nu0_t = 5.0, mu0_t = c(2.0, 1.0), Lambda0_t = L0,
#'   kappa0_c = 2.0, nu0_c = 5.0, mu0_c = c(1.0, 0.5), Lambda0_c = L0,
#'   r = NULL,
#'   ne_t = NULL, ne_c = 10L, alpha0e_t = NULL, alpha0e_c = 0.5,
#'   bar_ye_t = NULL, bar_ye_c = c(1.5, 0.8), se_t = NULL, se_c = se_c,
#'   nMC = 1000L
#' )
#'
#' # Example 5: Predictive probability, controlled design, vague prior
#' S_t <- matrix(c(18.0, 3.6, 3.6, 9.0), 2, 2)
#' S_c <- matrix(c(16.0, 2.8, 2.8, 8.5), 2, 2)
#' pbayespostpred2cont(
#'   prob = 'predictive', design = 'controlled', prior = 'vague',
#'   theta_TV1 = NULL, theta_MAV1 = NULL,
#'   theta_TV2 = NULL, theta_MAV2 = NULL,
#'   theta_NULL1 = 0.5, theta_NULL2 = 0.3,
#'   n_t = 12L, n_c = 12L,
#'   ybar_t = c(3.5, 2.1), S_t = S_t,
#'   ybar_c = c(1.8, 1.0), S_c = S_c,
#'   m_t = 30L, m_c = 30L,
#'   kappa0_t = NULL, nu0_t = NULL, mu0_t = NULL, Lambda0_t = NULL,
#'   kappa0_c = NULL, nu0_c = NULL, mu0_c = NULL, Lambda0_c = NULL,
#'   r = NULL,
#'   ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
#'   bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
#'   nMC = 1000L
#' )
#'
#' # Example 6: Posterior probability, controlled design, NIW prior, MM method
#' S_t <- matrix(c(18.0, 3.6, 3.6, 9.0), 2, 2)
#' S_c <- matrix(c(16.0, 2.8, 2.8, 8.5), 2, 2)
#' L0 <- matrix(c(20.0, 0.0, 0.0, 10.0), 2, 2)
#' pbayespostpred2cont(
#'   prob = 'posterior', design = 'controlled', prior = 'N-Inv-Wishart',
#'   theta_TV1 = 1.5, theta_MAV1 = 0.5,
#'   theta_TV2 = 1.0, theta_MAV2 = 0.3,
#'   theta_NULL1 = NULL, theta_NULL2 = NULL,
#'   n_t = 12L, n_c = 12L,
#'   ybar_t = c(3.5, 2.1), S_t = S_t,
#'   ybar_c = c(1.8, 1.0), S_c = S_c,
#'   m_t = NULL, m_c = NULL,
#'   kappa0_t = 2.0, nu0_t = 5.0, mu0_t = c(2.0, 1.0), Lambda0_t = L0,
#'   kappa0_c = 2.0, nu0_c = 5.0, mu0_c = c(1.0, 0.5), Lambda0_c = L0,
#'   r = NULL,
#'   ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
#'   bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
#'   nMC = 1000L, method = 'MM'
#' )
#'
#' # Example 7: Vectorised call -- N = 3 observations simultaneously
#' S_t_list <- list(
#'   matrix(c(18.0, 3.6, 3.6, 9.0), 2, 2),
#'   matrix(c(20.0, 4.0, 4.0, 10.0), 2, 2),
#'   matrix(c(15.0, 2.5, 2.5, 7.5), 2, 2)
#' )
#' S_c_list <- list(
#'   matrix(c(16.0, 2.8, 2.8, 8.5), 2, 2),
#'   matrix(c(17.0, 3.0, 3.0, 9.0), 2, 2),
#'   matrix(c(14.0, 2.2, 2.2, 7.0), 2, 2)
#' )
#' ybar_t_mat <- rbind(c(3.5, 2.1), c(4.0, 2.5), c(2.5, 1.5))
#' ybar_c_mat <- rbind(c(1.8, 1.0), c(1.9, 1.1), c(1.7, 0.9))
#' pbayespostpred2cont(
#'   prob = 'posterior', design = 'controlled', prior = 'vague',
#'   theta_TV1 = 1.5, theta_MAV1 = 0.5,
#'   theta_TV2 = 1.0, theta_MAV2 = 0.3,
#'   theta_NULL1 = NULL, theta_NULL2 = NULL,
#'   n_t = 12L, n_c = 12L,
#'   ybar_t = ybar_t_mat, S_t = S_t_list,
#'   ybar_c = ybar_c_mat, S_c = S_c_list,
#'   m_t = NULL, m_c = NULL,
#'   kappa0_t = NULL, nu0_t = NULL, mu0_t = NULL, Lambda0_t = NULL,
#'   kappa0_c = NULL, nu0_c = NULL, mu0_c = NULL, Lambda0_c = NULL,
#'   r = NULL,
#'   ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
#'   bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
#'   nMC = 1000L
#' )
#'
#' @importFrom stats rnorm rchisq
#' @importFrom mvtnorm pmvt
#' @export
pbayespostpred2cont <- function(prob,
                                design,
                                prior,
                                theta_TV1   = NULL, theta_MAV1  = NULL,
                                theta_TV2   = NULL, theta_MAV2  = NULL,
                                theta_NULL1 = NULL, theta_NULL2 = NULL,
                                n_t, n_c      = NULL,
                                ybar_t, S_t,
                                ybar_c       = NULL, S_c          = NULL,
                                m_t          = NULL, m_c          = NULL,
                                kappa0_t     = NULL, nu0_t        = NULL,
                                mu0_t        = NULL, Lambda0_t    = NULL,
                                kappa0_c     = NULL, nu0_c        = NULL,
                                mu0_c        = NULL, Lambda0_c    = NULL,
                                r           = NULL,
                                ne_t         = NULL, ne_c         = NULL,
                                alpha0e_t    = NULL, alpha0e_c    = NULL,
                                bar_ye_t     = NULL, bar_ye_c     = NULL,
                                se_t         = NULL, se_c         = NULL,
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
    for (nm in c('theta_TV1', 'theta_MAV1', 'theta_TV2', 'theta_MAV2')) {
      val <- get(nm)
      if (is.null(val) || !is.numeric(val) || length(val) != 1L || is.na(val))
        stop(paste0("'", nm, "' must be a single numeric when prob = 'posterior'"))
    }
    if (theta_TV1 <= theta_MAV1)
      stop("'theta_TV1' must be strictly greater than 'theta_MAV1'")
    if (theta_TV2 <= theta_MAV2)
      stop("'theta_TV2' must be strictly greater than 'theta_MAV2'")
  } else {
    for (nm in c('theta_NULL1', 'theta_NULL2')) {
      val <- get(nm)
      if (is.null(val) || !is.numeric(val) || length(val) != 1L || is.na(val))
        stop(paste0("'", nm, "' must be a single numeric when prob = 'predictive'"))
    }
    for (nm in c('m_t', 'm_c')) {
      val <- get(nm)
      if (is.null(val) || !is.numeric(val) || length(val) != 1L || is.na(val) ||
          val != floor(val) || val < 1L)
        stop(paste0("'", nm, "' must be a positive integer when prob = 'predictive'"))
    }
    m_t <- as.integer(m_t); m_c <- as.integer(m_c)
  }

  if (!is.numeric(n_t) || length(n_t) != 1L || is.na(n_t) ||
      n_t != floor(n_t) || n_t < 2L)
    stop("'n_t' must be a single integer >= 2")
  n_t <- as.integer(n_t)

  if (design %in% c('controlled', 'external')) {
    if (is.null(n_c) || !is.numeric(n_c) || length(n_c) != 1L || is.na(n_c) ||
        n_c != floor(n_c) || n_c < 2L)
      stop("'n_c' must be a single integer >= 2 for controlled/external design")
    n_c <- as.integer(n_c)
  }

  # Determine whether vectorised (N > 1) or scalar call
  vec_call <- is.matrix(ybar_t) && nrow(ybar_t) > 1L
  if (vec_call) {
    if (ncol(ybar_t) != 2L)
      stop("When 'ybar_t' is a matrix it must have exactly 2 columns")
    N <- nrow(ybar_t)
    if (!is.list(S_t) || length(S_t) != N)
      stop("When 'ybar_t' is a matrix, 'S_t' must be a list of the same length")
    if (design %in% c('controlled', 'external')) {
      if (is.null(ybar_c) || !is.matrix(ybar_c) || nrow(ybar_c) != N ||
          ncol(ybar_c) != 2L)
        stop("When 'ybar_t' is a matrix, 'ybar_c' must be a matrix with the same number of rows")
      if (!is.list(S_c) || length(S_c) != N)
        stop("When 'ybar_t' is a matrix, 'S_c' must be a list of the same length")
    }
  } else {
    N <- 1L
    # Coerce scalar ybar_t to length-2 vector
    if (!is.numeric(ybar_t) || length(ybar_t) != 2L || any(is.na(ybar_t)))
      stop("'ybar_t' must be a numeric vector of length 2 or a matrix with 2 columns")
    if (!is.matrix(S_t) || !is.numeric(S_t) || nrow(S_t) != 2L || ncol(S_t) != 2L)
      stop("'S_t' must be a 2x2 numeric matrix (or a list for vectorised call)")
    # Wrap as list for uniform handling below
    ybar_t <- matrix(ybar_t, nrow = 1L)
    S_t    <- list(S_t)
    if (design %in% c('controlled', 'external')) {
      if (!is.numeric(ybar_c) || length(ybar_c) != 2L || any(is.na(ybar_c)))
        stop("'ybar_c' must be a numeric vector of length 2")
      if (!is.matrix(S_c) || !is.numeric(S_c) || nrow(S_c) != 2L || ncol(S_c) != 2L)
        stop("'S_c' must be a 2x2 numeric matrix")
      ybar_c <- matrix(ybar_c, nrow = 1L)
      S_c    <- list(S_c)
    }
  }

  if (prior == 'N-Inv-Wishart') {
    for (nm in c('kappa0_t', 'nu0_t')) {
      val <- get(nm)
      if (is.null(val) || !is.numeric(val) || length(val) != 1L ||
          is.na(val) || val <= 0)
        stop(paste0("'", nm, "' must be a single positive numeric for NIW prior"))
    }
    if (is.null(mu0_t) || !is.numeric(mu0_t) || length(mu0_t) != 2L)
      stop("'mu0_t' must be a length-2 numeric vector for NIW prior")
    if (is.null(Lambda0_t) || !is.matrix(Lambda0_t) || nrow(Lambda0_t) != 2L)
      stop("'Lambda0_t' must be a 2x2 numeric matrix for NIW prior")
    if (design %in% c('controlled', 'external')) {
      for (nm in c('kappa0_c', 'nu0_c')) {
        val <- get(nm)
        if (is.null(val) || !is.numeric(val) || length(val) != 1L ||
            is.na(val) || val <= 0)
          stop(paste0("'", nm, "' must be a single positive numeric for NIW prior"))
      }
      if (is.null(mu0_c) || !is.numeric(mu0_c) || length(mu0_c) != 2L)
        stop("'mu0_c' must be a length-2 numeric vector for NIW prior")
      if (is.null(Lambda0_c) || !is.matrix(Lambda0_c) || nrow(Lambda0_c) != 2L)
        stop("'Lambda0_c' must be a 2x2 numeric matrix for NIW prior")
    }
  }

  if (design == 'uncontrolled') {
    if (is.null(r) || !is.numeric(r) || length(r) != 1L || is.na(r) || r <= 0)
      stop("'r' must be a single positive numeric when design = 'uncontrolled'")
    if (is.null(mu0_c) || !is.numeric(mu0_c) || length(mu0_c) != 2L)
      stop("'mu0_c' must be a length-2 numeric vector when design = 'uncontrolled'")
  }

  if (design == 'external') {
    has_ext_t <- !is.null(ne_t) && !is.null(alpha0e_t) &&
      !is.null(bar_ye_t) && !is.null(se_t)
    has_ext_c <- !is.null(ne_c) && !is.null(alpha0e_c) &&
      !is.null(bar_ye_c) && !is.null(se_c)
    if (!has_ext_t && !has_ext_c)
      stop(paste0("For design = 'external', at least one complete set of ",
                  "external data must be provided"))
    if (!is.null(ne_t)) ne_t <- as.integer(ne_t)
    if (!is.null(ne_c)) ne_c <- as.integer(ne_c)
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
    if (!is.null(ne_t)) {
      kappa_e_t  <- alpha0e_t * ne_t + kappa0_t
      nu_e_t     <- alpha0e_t * ne_t + nu0_t
      mu_e_t     <- (alpha0e_t * ne_t * bar_ye_t + kappa0_t * mu0_t) / kappa_e_t
      Lambda_e_t <- alpha0e_t * se_t + Lambda0_t +
        (kappa0_t * alpha0e_t * ne_t / kappa_e_t) * tcrossprod(mu0_t - bar_ye_t)
    }
    if (!is.null(ne_c)) {
      kappa_e_c  <- alpha0e_c * ne_c + kappa0_c
      nu_e_c     <- alpha0e_c * ne_c + nu0_c
      mu_e_c     <- (alpha0e_c * ne_c * bar_ye_c + kappa0_c * mu0_c) / kappa_e_c
      Lambda_e_c <- alpha0e_c * se_c + Lambda0_c +
        (kappa0_c * alpha0e_c * ne_c / kappa_e_c) * tcrossprod(mu0_c - bar_ye_c)
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
    df_t <- as.integer(nu0_t + n_t - 1L)
    if (design %in% c('controlled', 'external')) df_c <- as.integer(nu0_c + n_c - 1L)
    if (design == 'uncontrolled')                df_c <- df_t
  } else {
    df_t <- n_t - 2L
    if (design %in% c('controlled', 'external')) df_c <- n_c - 2L
    if (design == 'uncontrolled')                df_c <- df_t
  }
  if (design == 'external') {
    if (!is.null(ne_t)) df_t <- as.integer(nu_e_t + n_t - 1L)
    if (!is.null(ne_c)) df_c <- as.integer(nu_e_c + n_c - 1L)
  }

  # Check MM feasibility
  use_mm <- (method == 'MM')
  if (use_mm && (df_t <= 4L || df_c <= 4L)) {
    warning(
      "MM method requires df > 4 for both arms (df_t = ", df_t,
      ", df_c = ", df_c, "). Falling back to method = 'MC'."
    )
    use_mm <- FALSE
  }

  # Pre-generate MC raw variates (only for MC path; also used as fallback for MM)
  if (!use_mm) {
    if (is.null(nMC))
      stop(paste0("'nMC' must be a positive integer when method = 'MM' falls ",
                  "back to MC (nu_k <= 4)"))
    Z_t <- matrix(rnorm(nMC * 2L), nrow = nMC, ncol = 2L)
    W_t <- rchisq(nMC, df = df_t)
    sc_t <- sqrt(W_t / df_t)

    Z_c <- matrix(rnorm(nMC * 2L), nrow = nMC, ncol = 2L)
    W_c <- rchisq(nMC, df = df_c)
    sc_c <- sqrt(W_c / df_c)
  }

  # ---------------------------------------------------------------------------
  # Section 5: Region probability computation for N observations
  # ---------------------------------------------------------------------------

  n_regions <- if (prob == 'posterior') 9L else 4L
  Pr_R_mat  <- matrix(0, nrow = N, ncol = n_regions)

  # Cut-points (fixed across all observations)
  if (prob == 'posterior') {
    cuts_t <- c(-Inf, theta_MAV1, theta_TV1, Inf)
    cuts_c <- c(-Inf, theta_MAV2, theta_TV2, Inf)
  } else {
    cuts_t <- c(-Inf, theta_NULL1, Inf)
    cuts_c <- c(-Inf, theta_NULL2, Inf)
  }

  for (i in seq_len(N)) {

    yb_t <- ybar_t[i, ]
    S_t_i <- S_t[[i]]

    # Posterior parameters for group 1
    if (design == 'external' && !is.null(ne_t)) {
      post_t <- .niw_post(n_t, yb_t, S_t_i,
                          kappa0 = kappa_e_t, nu0 = nu_e_t,
                          mu0 = mu_e_t, Lambda0 = Lambda_e_t,
                          for_predictive = use_pred)
    } else if (prior == 'N-Inv-Wishart') {
      post_t <- .niw_post(n_t, yb_t, S_t_i,
                          kappa0 = kappa0_t, nu0 = nu0_t,
                          mu0 = mu0_t, Lambda0 = Lambda0_t,
                          for_predictive = use_pred)
    } else {
      post_t <- .vague_post(n_t, yb_t, S_t_i, for_predictive = use_pred)
    }

    # Posterior parameters for group 2
    if (design == 'uncontrolled') {
      post_c <- list(df = post_t$df, mu_n = mu0_c, V = r * post_t$V)
    } else if (design == 'external' && !is.null(ne_c)) {
      yb_c <- ybar_c[i, ]
      S_c_i <- S_c[[i]]
      post_c <- .niw_post(n_c, yb_c, S_c_i,
                          kappa0 = kappa_e_c, nu0 = nu_e_c,
                          mu0 = mu_e_c, Lambda0 = Lambda_e_c,
                          for_predictive = use_pred)
    } else if (prior == 'N-Inv-Wishart') {
      yb_c <- ybar_c[i, ]
      S_c_i <- S_c[[i]]
      post_c <- .niw_post(n_c, yb_c, S_c_i,
                          kappa0 = kappa0_c, nu0 = nu0_c,
                          mu0 = mu0_c, Lambda0 = Lambda0_c,
                          for_predictive = use_pred)
    } else {
      yb_c <- ybar_c[i, ]
      S_c_i <- S_c[[i]]
      post_c <- .vague_post(n_c, yb_c, S_c_i, for_predictive = use_pred)
    }

    if (!use_mm) {

      # -----------------------------------------------------------------------
      # MC path: re-use pre-generated Z and W; only chol(V) changes per obs
      # -----------------------------------------------------------------------

      V_t_use <- if (use_pred) post_t$V / m_t else post_t$V
      V_c_use <- if (use_pred) post_c$V / m_c else post_c$V

      R_chol_t <- chol(V_t_use + diag(1e-10, 2L))
      R_chol_c <- chol(V_c_use + diag(1e-10, 2L))

      mu_t <- sweep(Z_t %*% R_chol_t, 1L, sc_t, '/') +
        matrix(post_t$mu_n, nrow = nMC, ncol = 2L, byrow = TRUE)
      mu_c <- sweep(Z_c %*% R_chol_c, 1L, sc_c, '/') +
        matrix(post_c$mu_n, nrow = nMC, ncol = 2L, byrow = TRUE)

      theta1 <- mu_t[, 1L] - mu_c[, 1L]
      theta2 <- mu_t[, 2L] - mu_c[, 2L]

      if (prob == 'posterior') {
        r1 <- 3L - as.integer(theta1 > theta_MAV1) -
          as.integer(theta1 > theta_TV1)
        r2 <- 3L - as.integer(theta2 > theta_MAV2) -
          as.integer(theta2 > theta_TV2)
        region <- (r1 - 1L) * 3L + r2
        Pr_R_mat[i, ] <- tabulate(region, nbins = 9L) / nMC
      } else {
        r1 <- 2L - as.integer(theta1 > theta_NULL1)
        r2 <- 2L - as.integer(theta2 > theta_NULL2)
        region <- (r1 - 1L) * 2L + r2
        Pr_R_mat[i, ] <- tabulate(region, nbins = 4L) / nMC
      }

    } else {

      # -----------------------------------------------------------------------
      # MM path: moment-matching + mvtnorm::pmvt per observation
      # -----------------------------------------------------------------------

      V_t_use <- if (use_pred) post_t$V / m_t else post_t$V
      V_c_use <- if (use_pred) post_c$V / m_c else post_c$V

      mm <- .mm_general(post_t$mu_n, V_t_use, post_t$df,
                        post_c$mu_n, V_c_use, post_c$df)

      n_bands <- length(cuts_t) - 1L
      Pr_R <- numeric(n_regions)
      idx  <- 1L
      for (ii in n_bands:1L) {
        for (jj in n_bands:1L) {
          Pr_R[idx] <- .rect_prob(cuts_t[ii], cuts_t[ii + 1L],
                                  cuts_c[jj], cuts_c[jj + 1L], mm)
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
