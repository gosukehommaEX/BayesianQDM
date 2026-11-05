#' Operating characteristics of several designs
#'
#' Computes the operating characteristics of several designs under the same
#' scenarios, for example designs with different sample sizes, prior
#' distributions, external data or thresholds, and collects them in one
#' table with a column that names the design.
#'
#' Each design is computed by \code{\link{qdm_oc}} with the same
#' \code{scenarios}, \code{nsim} and \code{seed}, so the rows of a design
#' are those of \code{qdm_oc()} for that design alone. With a seed, the
#' simulated studies of all designs start from the same random numbers. The
#' warnings of \code{qdm_oc()} (Miss decisions, scenarios that are not
#' evaluable) start with the name of the design.
#'
#' Designs that differ in one part can be made from one design with the
#' design steps (see \code{\link{qdm_endpoint}}), which accept a design; for
#' example, for several sample sizes,
#' \code{lapply(c(12, 18, 24), function(n) qdm_sample_size(des, n, n))}.
#'
#' @param designs A list of designs (created by a \code{design_*()} function
#'   such as \code{\link{design_binary}}) or design specifications (started
#'   by \code{\link{qdm_endpoint}}), all for the same endpoint. The names of
#'   the list name the designs; designs without a name are called
#'   \code{"Design 1"}, \code{"Design 2"} and so on, by their position.
#' @param scenarios A data frame of scenarios with the columns needed by all
#'   the designs (see \code{\link{qdm_oc}} and \code{\link{qdm_scenarios}}).
#'   It must not have a column named \code{design}.
#' @param nsim,seed As in \code{\link{qdm_oc}}.
#'
#' @return An object of class \code{c("qdm_comparison", "data.frame")}: the
#'   column \code{design} (a factor with the names of the designs as levels,
#'   in the order of \code{designs}), followed by the columns of the results
#'   of \code{\link{qdm_oc}} (the scenarios and the probabilities of the
#'   decisions). The attribute \code{designs} is the list of the designs
#'   (built from the specifications) and the attribute \code{method} the
#'   methods of computation, named by the designs.
#'
#' @family operations
#' @examples
#' des <- design_binary(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_beta(0.5, 0.5), prior_c = prior_beta(0.5, 0.5),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                   gamma_go = 0.5, gamma_nogo = 0.8)
#' )
#' n <- c(12, 18, 24)
#' designs <- lapply(n, function(k) qdm_sample_size(des, n_t = k, n_c = k))
#' names(designs) <- paste(n, "per group")
#' cmp <- qdm_compare(designs, qdm_scenarios(des, theta = seq(0, 0.4, by = 0.1),
#'                                           pi_c = 0.3))
#' cmp
#' plot(cmp)
#'
#' # Borrowing of external control data
#' ext <- qdm_external(des, external_c = borrow(data_binary(40, 12), 0.5))
#' qdm_compare(list(Concurrent = des, Borrowing = ext),
#'             data.frame(pi_t = c(0.3, 0.6), pi_c = 0.3))
#' @export
qdm_compare <- function(designs, scenarios, nsim = 10000L, seed = NULL) {
  if (!is.list(designs) || inherits(designs, c("qdm_design", "qdm_spec")) ||
      length(designs) == 0L) {
    stop("'designs' must be a list of designs or design specifications.",
         call. = FALSE)
  }
  ok <- vapply(designs, inherits, logical(1),
               what = c("qdm_design", "qdm_spec"))
  if (!all(ok)) {
    stop(paste("Each element of 'designs' must be a design (created by a",
               "design_*() function) or a design specification (started by",
               "qdm_endpoint())."), call. = FALSE)
  }
  labels <- names(designs)
  if (is.null(labels)) {
    labels <- rep("", length(designs))
  }
  unnamed <- is.na(labels) | labels == ""
  labels[unnamed] <- paste("Design", which(unnamed))
  if (anyDuplicated(labels) > 0L) {
    stop("The names of the designs must be unique.", call. = FALSE)
  }
  built <- stats::setNames(lapply(designs, check_design), labels)
  endpoints <- unique(vapply(built, function(d) d$endpoint, character(1)))
  if (length(endpoints) > 1L) {
    stop("The designs must have the same endpoint.", call. = FALSE)
  }
  if (is.data.frame(scenarios) && "design" %in% names(scenarios)) {
    stop("'scenarios' must not have a column named design.", call. = FALSE)
  }
  results <- lapply(seq_along(built), function(k) {
    withCallingHandlers(
      qdm_oc(built[[k]], scenarios, nsim = nsim, seed = seed),
      warning = function(w) {
        warning(sprintf("%s: %s", labels[k], conditionMessage(w)),
                call. = FALSE)
        invokeRestart("muffleWarning")
      }
    )
  })
  tables <- lapply(seq_along(results), function(k) {
    r <- results[[k]]
    attr(r, "design") <- NULL
    attr(r, "method") <- NULL
    class(r) <- "data.frame"
    cbind(data.frame(design = rep(labels[k], nrow(r))), r)
  })
  out <- do.call(rbind, tables)
  out$design <- factor(out$design, levels = labels)
  rownames(out) <- NULL
  structure(out, class = c("qdm_comparison", "data.frame"), designs = built,
            method = stats::setNames(vapply(results, attr, character(1),
                                            which = "method"), labels))
}
