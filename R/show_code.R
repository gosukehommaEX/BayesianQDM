#' R code of a design
#'
#' Writes the R code that creates a design: a pipe of the design steps that
#' starts with \code{\link{qdm_endpoint}} and ends with
#' \code{\link{qdm_build}}, or a call of the \code{design_*()} function of
#' the endpoint. Evaluating the code gives a design identical to the
#' original. The code can be kept with a study protocol or a report (see the
#' argument \code{code} of \code{\link{show_report}}), or used as the start
#' of a new analysis.
#'
#' The arguments are written with the values that the design keeps (the
#' numbers with the fewest digits that give them exactly), and the arguments
#' equal to their defaults are left out. For a design specification, the
#' code of the pipe gives the specification (without \code{qdm_build()}).
#' For a result of an analysis, the code creates the design of the result;
#' the scenarios, the number of simulated studies and the seed of the
#' analysis are not part of the code.
#'
#' @param x A design (created by a \code{design_*()} function such as
#'   \code{\link{design_binary}}), a design specification (started by
#'   \code{\link{qdm_endpoint}}), or a result of \code{\link{qdm_oc}},
#'   \code{\link{qdm_calibrate}}, \code{\link{qdm_find_n}},
#'   \code{\link{qdm_boundary}} or \code{\link{qdm_decide}}.
#' @param style \code{"pipe"} (the default) for the pipe of the design
#'   steps, or \code{"design"} for the call of the \code{design_*()}
#'   function.
#' @param name Name to which the code assigns the design, or \code{NULL} for
#'   no assignment.
#'
#' @return An object of class \code{"qdm_code"}: a character vector with one
#'   line of code per element, which is printed as the code.
#'
#' @examples
#' des <- design_binary(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_beta(0.5, 0.5), prior_c = prior_beta(0.5, 0.5),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                   gamma_go = 0.5, gamma_nogo = 0.8),
#'   external_c = borrow(data_binary(n = 40, y = 12), a0 = 0.5)
#' )
#' show_code(des)
#' show_code(des, style = "design")
#'
#' # The code gives the same design
#' code <- show_code(des, name = NULL)
#' identical(eval(parse(text = code)), des)
#' @export
show_code <- function(x, style = c("pipe", "design"), name = "design") {
  style <- match.arg(style)
  if (!is.null(name) && (!is.character(name) || length(name) != 1L ||
                         is.na(name) || make.names(name) != name)) {
    stop("'name' must be NULL or a syntactic name, such as \"design\".",
         call. = FALSE)
  }
  target <- code_target(x)
  spec <- if (inherits(target, "qdm_spec")) target else design_spec(target)
  args <- code_design_args(spec)
  vals <- lapply(names(args), function(nm) {
    code_value(args[[nm]], keep_type = !(nm %in% c("n_t", "n_c", "nMC")))
  })
  names(vals) <- names(args)
  prefix <- if (is.null(name)) "" else paste(name, "<- ")
  if (style == "design") {
    lines <- code_format(code_node(design_function(spec$endpoint), vals),
                         nchar(prefix), hanging = TRUE)
    lines[1L] <- paste0(prefix, lines[1L])
    return(structure(lines, class = "qdm_code"))
  }
  # Arguments of the design_*() functions set by each design step, in the
  # order of the pipe
  step_args <- list(
    qdm_sample_size = c("n_t", "n_c"),
    qdm_prior = c("prior_t", "prior_c"),
    qdm_criterion = "criterion",
    qdm_external = c("external_t", "external_c"),
    qdm_control = c("hypothetical_c", "fixed_c"),
    qdm_tte = c("milestone", "effect", "conduct", "conduct_future"),
    qdm_computation = c("method", "nMC", "enumerate")
  )
  steps <- list(code_node("qdm_endpoint", list(code_value(spec$endpoint))))
  for (step in names(step_args)) {
    present <- intersect(step_args[[step]], names(vals))
    if (length(present) == 0L) {
      next
    }
    a <- vals[present]
    if (step == "qdm_criterion") {
      names(a) <- ""
    }
    steps <- c(steps, list(code_node(step, a)))
  }
  if (inherits(target, "qdm_design")) {
    steps <- c(steps, list(code_node("qdm_build")))
  }
  lines <- character(0)
  n <- length(steps)
  for (k in seq_len(n)) {
    lead <- if (k == 1L) prefix else "  "
    # Room for " |>" at the end of all steps but the last
    part <- code_format(steps[[k]], nchar(lead),
                        width = if (k < n) 77L else 80L)
    part[1L] <- paste0(lead, part[1L])
    if (k < n) {
      part[length(part)] <- paste0(part[length(part)], " |>")
    }
    lines <- c(lines, part)
  }
  structure(lines, class = "qdm_code")
}

