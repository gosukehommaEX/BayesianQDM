#' Hypothetical control for an uncontrolled study
#'
#' Specifies the distribution of a hypothetical control group, used when the
#' proof-of-concept study has no concurrent control group (Kang et al., 2026).
#'
#' For binary and time-to-event endpoints, give \code{data}: the
#' hypothetical control data are combined with the prior of the control
#' group, for example \eqn{\pi_c \sim \textrm{Beta}(a_c + y, b_c + n - y)}
#' for one binary endpoint, or a Gamma distribution of the hazard with shape
#' \eqn{a_c + d} and rate \eqn{b_c + E} for \eqn{d} events in the total
#' observed time \eqn{E}. For continuous endpoints, give \code{mean} and \code{r}: the
#' control mean has the same t distribution as the posterior distribution of
#' the treatment mean, with location \code{mean} and squared scale multiplied
#' by \code{r}, that is
#' \eqn{\mu_c \sim t_{\nu_{nt}}(\mu_{0c}, r \sigma_{nt}^2 / \kappa_{nt})}.
#'
#' @param data Hypothetical control data created by
#'   \code{\link{data_binary}}, \code{\link{data_binary2}} (binary
#'   endpoints) or \code{\link{data_tte}} (time-to-event endpoint).
#' @param mean Location of the hypothetical control mean, a number for one
#'   continuous endpoint or a numeric vector of length 2 for two continuous
#'   endpoints.
#' @param r Positive ratio of the squared scale of the hypothetical control to
#'   that of the treatment group (continuous endpoints).
#'
#' @return An object of class \code{"qdm_hypothetical"} with elements
#'   \code{type} (\code{"data"} or \code{"location"}), \code{data},
#'   \code{mean} and \code{r}.
#'
#' @references
#' Kang, R., Yamaguchi, Y. and Han, C. (2026). Bayesian quantitative
#' decision-making framework in cell and gene therapy development for rare
#' diseases. \emph{Journal of Biopharmaceutical Statistics}.
#' \doi{10.1080/10543406.2026.2655410}
#'
#' @examples
#' # Binary endpoint: 5 responders out of 15 hypothetical control subjects
#' hypothetical(data = data_binary(n = 15, y = 5))
#'
#' # Continuous endpoint: control mean located at 1 with the same scale
#' hypothetical(mean = 1, r = 1)
#'
#' # Time-to-event endpoint: 8 events in a total observed time of 60
#' hypothetical(data = data_tte(events = 8, exposure = 60))
#' @export
hypothetical <- function(data = NULL, mean = NULL, r = NULL) {
  has_data <- !is.null(data)
  has_location <- !is.null(mean) || !is.null(r)
  if (has_data == has_location) {
    stop(paste("Specify either 'data' (binary and time-to-event endpoints)",
               "or both 'mean' and 'r' (continuous endpoints)."),
         call. = FALSE)
  }
  if (has_data) {
    if (!inherits(data, c("qdm_data_binary", "qdm_data_binary2",
                          "qdm_data_tte"))) {
      stop(paste("'data' must be created by data_binary(), data_binary2()",
                 "or data_tte()."), call. = FALSE)
    }
    out <- list(type = "data", data = data, mean = NULL, r = NULL)
  } else {
    if (is.null(mean) || is.null(r)) {
      stop("Both 'mean' and 'r' must be specified for continuous endpoints.",
           call. = FALSE)
    }
    if (!(length(mean) %in% c(1L, 2L))) {
      stop("'mean' must have length 1 or 2.", call. = FALSE)
    }
    check_numeric(mean, "mean", len = length(mean))
    check_numeric(r, "r", lower = 0, lower_open = TRUE)
    out <- list(type = "location", data = NULL, mean = as.numeric(mean),
                r = r)
  }
  structure(out, class = "qdm_hypothetical")
}
