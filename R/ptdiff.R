#' Distribution function of the difference of two t variables
#'
#' Computes \eqn{\Pr(T_t - T_c \le q)} (or the upper tail
#' \eqn{\Pr(T_t - T_c > q)}) for independent location-scale t variables
#' \eqn{T_k = \mu_k + s_k X_k}, where \eqn{X_k} follows a t distribution with
#' \eqn{\nu_k} degrees of freedom. This is the posterior (or posterior
#' predictive) probability used for a single continuous endpoint.
#'
#' Three methods are available.
#' \describe{
#'   \item{\code{"NI"}}{Numerical integration over the whole real line by
#'     adaptive Gauss-Kronrod quadrature (the QUADPACK routine used by
#'     \code{\link[stats]{integrate}}). The integral is taken over the
#'     variable of the group with the smaller scale, for example
#'     \eqn{\int f_c(u) \Pr(T_t > u + q) \, du} when
#'     \eqn{s_c \le s_t}. The range is split around the peak of the density
#'     and around the change of the other factor, and the parts far from the
#'     origin are integrated on a logarithmic scale, so that heavy tails and
#'     very different scales are handled. The absolute error is below
#'     \eqn{10^{-9}} in the tests, which include scale ratios from
#'     \eqn{10^{-9}} to \eqn{10^3} and degrees of freedom from 0.5 to
#'     \code{Inf}.}
#'   \item{\code{"MM"}}{The moment-matching approximation of Yamaguchi et al.
#'     (submitted, Theorem 1), which replaces the difference by a single t
#'     distribution with matched second and fourth moments. It requires
#'     \eqn{\nu_t > 4} and \eqn{\nu_c > 4}, is evaluated in closed form, and is
#'     the fastest method.}
#'   \item{\code{"MC"}}{Monte Carlo simulation with \code{nMC} draws for each
#'     element, using the random number generator of R.}
#' }
#' Degrees of freedom may be \code{Inf} (normal distributions). For the lower
#' tail, the upper tail of \eqn{T_c - T_t} at \eqn{-q} is computed, so both
#' tails are accurate near 0.
#'
#' @param q Numeric vector of quantiles.
#' @param mu_t,scale_t,df_t Location, positive scale and positive degrees of
#'   freedom of the treatment group.
#' @param mu_c,scale_c,df_c Location, positive scale and positive degrees of
#'   freedom of the control group.
#' @param method One of \code{"NI"} (default), \code{"MM"} or \code{"MC"}.
#' @param nMC Number of Monte Carlo draws per element for \code{method = "MC"}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are
#'   \eqn{\Pr(T_t - T_c \le q)}, otherwise \eqn{\Pr(T_t - T_c > q)}.
#'
#' @return A numeric vector of probabilities. All arguments from \code{q} to
#'   \code{df_c} are recycled to a common length; each must have length 1 or
#'   the common length.
#'
#' @family distributions
#' @examples
#' ptdiff(1, mu_t = 3, scale_t = 1, df_t = 15, mu_c = 1, scale_c = 1,
#'        df_c = 15, lower.tail = FALSE)
#' ptdiff(1, mu_t = 3, scale_t = 1, df_t = 15, mu_c = 1, scale_c = 1,
#'        df_c = 15, method = "MM", lower.tail = FALSE)
#' @export
ptdiff <- function(q, mu_t, scale_t, df_t, mu_c, scale_c, df_c,
                   method = c("NI", "MM", "MC"), nMC = 10000L,
                   lower.tail = TRUE) {
  method <- match.arg(method)
  check_numeric(q, "q", len = NULL)
  check_numeric(mu_t, "mu_t", len = NULL)
  check_numeric(scale_t, "scale_t", len = NULL, lower = 0, lower_open = TRUE)
  check_numeric(df_t, "df_t", len = NULL, lower = 0, lower_open = TRUE,
                allow_inf = TRUE)
  check_numeric(mu_c, "mu_c", len = NULL)
  check_numeric(scale_c, "scale_c", len = NULL, lower = 0, lower_open = TRUE)
  check_numeric(df_c, "df_c", len = NULL, lower = 0, lower_open = TRUE,
                allow_inf = TRUE)
  check_flag(lower.tail, "lower.tail")
  args <- recycle_common(list(q = q, mu_t = mu_t, scale_t = scale_t,
                              df_t = df_t, mu_c = mu_c, scale_c = scale_c,
                              df_c = df_c))
  if (method == "MM" && any(args$df_t <= 4 | args$df_c <= 4)) {
    stop("method = \"MM\" requires 'df_t' and 'df_c' greater than 4.",
         call. = FALSE)
  }
  if (method == "MC") {
    check_numeric(nMC, "nMC", lower = 1, upper = .Machine$integer.max,
                  integer = TRUE)
  }
  if (lower.tail) {
    # P(T_t - T_c <= q) = P(T_c - T_t >= -q), the upper tail with the groups
    # exchanged (the distribution is continuous).
    args <- list(q = -args$q, mu_t = args$mu_c, scale_t = args$scale_c,
                 df_t = args$df_c, mu_c = args$mu_t, scale_c = args$scale_t,
                 df_c = args$df_t)
  }
  switch(
    method,
    NI = cpp_ptdiff_ni_upper(args$q, args$mu_t, args$scale_t, args$df_t,
                             args$mu_c, args$scale_c, args$df_c),
    MM = ptdiff_mm_upper(args$q, args$mu_t, args$scale_t, args$df_t,
                         args$mu_c, args$scale_c, args$df_c),
    MC = cpp_ptdiff_mc_upper(args$q, args$mu_t, args$scale_t, args$df_t,
                             args$mu_c, args$scale_c, args$df_c,
                             as.integer(round(nMC)))
  )
}
