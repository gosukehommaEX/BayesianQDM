#' Write a report of a design and its results
#'
#' Writes a report of a design, or of the result of an analysis of a design,
#' as a draft of the corresponding sections of a study protocol: the study
#' design, the statistical model and the prior distributions, the decision
#' criteria with their formulas, the computation of the probabilities and the
#' result. The report is meant to be copied into a protocol and edited there.
#' It can end a pipe of the design steps (see \code{\link{qdm_endpoint}}) and
#' of the functions that use a design.
#'
#' The content of the report depends on \code{x}:
#' \describe{
#'   \item{A design or a design specification}{The design.}
#'   \item{A result of \code{\link{qdm_oc}}}{The design and the table of the
#'     operating characteristics.}
#'   \item{A result of \code{\link{qdm_calibrate}}}{The design with the
#'     calibrated thresholds, how the thresholds were chosen, and the decision
#'     probabilities under the two scenarios of the calibration.}
#'   \item{A result of \code{\link{qdm_boundary}}}{The design and the decision
#'     boundaries.}
#'   \item{A result of \code{\link{qdm_decide}}}{The design, the observed
#'     data, the probabilities of the criteria, the decision and the posterior
#'     distributions.}
#' }
#'
#' The formulas are written in LaTeX between dollar signs. In the HTML
#' report they are displayed by MathJax, which the web browser loads from a
#' content delivery network, so an internet connection is needed to display
#' them (without one, the LaTeX code is shown). The text report is plain text
#' in Markdown format, with the formulas in LaTeX. A PDF file can be made by
#' printing the HTML report from the web browser.
#'
#' @param x A design (created by a \code{design_*()} function such as
#'   \code{\link{design_binary}}), a design specification (started by
#'   \code{\link{qdm_endpoint}}), or a result of \code{\link{qdm_oc}},
#'   \code{\link{qdm_calibrate}}, \code{\link{qdm_boundary}} or
#'   \code{\link{qdm_decide}}.
#' @param file_type \code{"html"} (the default) for an HTML file or
#'   \code{"text"} for a plain text file.
#' @param file Name of the file to write, or \code{NULL} (the default) for a
#'   temporary file with the extension \code{.html} or \code{.txt}.
#' @param title Title of the report, or \code{NULL} for a title that depends
#'   on \code{x}.
#' @param digits Number of decimal places of the probabilities (a whole
#'   number from 1 to 10). The parameters are shown with up to four
#'   significant digits.
#' @param open Whether to show the report: an HTML report is opened in the
#'   viewer of RStudio, if the file is in the temporary folder of the R
#'   session, or else in the web browser, and a text report is printed in
#'   the console. By default the report is shown when R is used
#'   interactively.
#'
#' @return The name of the file, invisibly.
#'
#' @examples
#' spec <- qdm_endpoint("binary") |>
#'   qdm_sample_size(n_t = 12, n_c = 12) |>
#'   qdm_prior(prior_t = prior_beta(0.5, 0.5),
#'             prior_c = prior_beta(0.5, 0.5)) |>
#'   qdm_criterion(criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                     gamma_go = 0.5, gamma_nogo = 0.8))
#'
#' # Text report of the operating characteristics, printed in the console
#' spec |>
#'   qdm_oc(data.frame(pi_t = c(0.3, 0.5), pi_c = 0.3)) |>
#'   show_report(file_type = "text", open = TRUE)
#'
#' # HTML report of a decision in a temporary file, which is opened when R
#' # is used interactively
#' path <- spec |>
#'   qdm_decide(data_t = data_binary(12, 9), data_c = data_binary(12, 4)) |>
#'   show_report()
#' @export
show_report <- function(x, file_type = c("html", "text"), file = NULL,
                        title = NULL, digits = 3L, open = interactive()) {
  file_type <- match.arg(file_type)
  if (!is.null(file) && (!is.character(file) || length(file) != 1L ||
                         is.na(file) || !nzchar(file))) {
    stop("'file' must be NULL or the name of a file.", call. = FALSE)
  }
  if (!is.null(title) && (!is.character(title) || length(title) != 1L ||
                          is.na(title))) {
    stop("'title' must be NULL or a character string.", call. = FALSE)
  }
  check_numeric(digits, "digits", lower = 1, upper = 10, integer = TRUE)
  check_flag(open, "open")
  blocks <- report_blocks(x, as.integer(round(digits)))
  if (is.null(title)) {
    title <- report_title(x)
  }
  note <- sprintf("Report written by BayesianQDM %s.",
                  utils::packageVersion("BayesianQDM"))
  lines <- if (file_type == "html") {
    report_html(title, blocks, note)
  } else {
    report_text(title, blocks, note)
  }
  if (is.null(file)) {
    file <- tempfile("qdm_report_",
                     fileext = if (file_type == "html") ".html" else ".txt")
  }
  writeLines(enc2utf8(lines), file, useBytes = TRUE)
  if (open) {
    if (file_type == "text") {
      cat(lines, sep = "\n")
    } else {
      open_html(file)
    }
  }
  invisible(file)
}
