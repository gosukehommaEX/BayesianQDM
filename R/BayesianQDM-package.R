#' BayesianQDM: Bayesian Quantitative Decision-Making for Clinical Trials
#'
#' Tools for Go/NoGo/Gray decision-making in proof-of-concept clinical trials
#' under a Bayesian framework. A trial is described by small specification
#' objects (priors, data, external data with power priors, hypothetical
#' controls and decision criteria), and the decision rule is evaluated for
#' observed data or for its operating characteristics.
#'
#' @keywords internal
#' @useDynLib BayesianQDM, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom ggplot2 .data
"_PACKAGE"
