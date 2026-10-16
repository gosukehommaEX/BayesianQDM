#' Format a decision boundary
#'
#' Formats the decision intervals of \code{\link{qdm_boundary}} as tables
#' with one row per condition and one column per decision: the range of the
#' observed quantity (numbers of responders or events, or observed means),
#' and the corresponding observed response rates (binary endpoint) or
#' milestone survival probabilities (time-to-event endpoint) and observed
#' treatment effects.
#'
#' @param x An object of class \code{"qdm_boundary"}.
#' @param digits Number of significant digits for observed means, milestone
#'   survival probabilities, effects and conditions.
#' @param n Maximum number of conditions shown in each table.
#' @param ... Not used.
#'
#' @return A character vector with one element per line of the tables
#'   (\code{\link{print.qdm_boundary}} prints them). For rows subset from
#'   the object, which no longer have its attributes, the formatted data
#'   frame.
#'
#' @examples
#' des <- design_binary(
#'   n_t = 20, n_c = 0, prior_t = prior_beta(1, 1),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.05,
#'                                   gamma_go = 0.8, gamma_nogo = 0.8),
#'   fixed_c = 0.2
#' )
#' format(qdm_boundary(des, n_t = c(10, 20)))
#' @export
format.qdm_boundary <- function(x, digits = 4L, n = 20L, ...) {
  axes <- attr(x, "axes")
  check_numeric(digits, "digits", lower = 1, integer = TRUE)
  check_numeric(n, "n", lower = 1, integer = TRUE)
  if (is.null(axes)) {
    return(format(as.data.frame(x), digits = digits))
  }
  titles <- attr(x, "titles")
  design <- attr(x, "design")
  discrete <- attr(x, "discrete")
  num <- function(v) vapply(v, format, character(1), digits = digits)
  cond_col <- axes$x
  conds <- unique(x[[cond_col]])
  shown <- conds[seq_len(min(n, length(conds)))]
  decisions <- c("NoGo", "Gray", "Miss", "Go")
  if (!attr(x, "increasing")) {
    decisions <- rev(decisions)
  }
  decisions <- decisions[decisions %in% as.character(x$decision)]
  # Text of an interval: a range of whole numbers, of percentages or of
  # values (with <= and >= for infinite ends)
  cell <- function(from, to, kind) {
    lo <- pmin(from, to)
    hi <- pmax(from, to)
    if (kind == "count") {
      return(ifelse(lo == hi, num(lo), paste0(num(lo), "-", num(hi))))
    }
    f <- if (kind == "percent") {
      function(v) {
        paste0(vapply(round(100 * v, 1), format, character(1), nsmall = 1),
               "%")
      }
    } else {
      num
    }
    ifelse(is.infinite(lo) & is.infinite(hi), "all",
           ifelse(is.infinite(lo), paste("<=", f(hi)),
                  ifelse(is.infinite(hi), paste(">=", f(lo)),
                         ifelse(lo == hi, f(lo),
                                paste(f(lo), "to", f(hi))))))
  }
  table <- function(from, to, kind) {
    txt <- cell(from, to, kind)
    rows <- lapply(shown, function(cv) {
      sel <- x[[cond_col]] == cv
      vapply(decisions, function(d) {
        k <- sel & x$decision == d
        if (any(k)) paste(txt[k], collapse = ", ") else "-"
      }, character(1))
    })
    m <- rbind(c(cond_col, decisions),
               cbind(num(shown), do.call(rbind, rows)))
    widths <- apply(nchar(m), 2L, max)
    lines <- apply(m, 1L, function(r) {
      paste0("  ", paste(sprintf("%-*s", widths, r), collapse = "   "))
    })
    sub("\\s+$", "", lines)
  }
  endpoint <- switch(design$endpoint,
                     binary = "one binary endpoint",
                     continuous = "one continuous endpoint",
                     continuous2 = "two continuous endpoints",
                     tte = if (identical(design$effect, "hazard_ratio")) {
                       "one time-to-event endpoint (hazard ratio)"
                     } else {
                       "one time-to-event endpoint (milestone survival)"
                     })
  out <- sprintf("Decision boundary of a design with %s", endpoint)
  fixed_cond <- setdiff(intersect(c("n_t", "n_c", "sd_t"), names(x)),
                        cond_col)
  if (length(fixed_cond) > 0L) {
    out <- c(out, paste0("  ", paste(fixed_cond, "=",
                                     vapply(fixed_cond, function(nm) {
                                       num(x[[nm]][1L])
                                     }, character(1)), collapse = ", ")))
  }
  if (!is.null(axes$caption)) {
    out <- c(out, paste0("  ", axes$caption))
  }
  kind_cut <- if (discrete) "count" else "value"
  out <- c(out, "", titles[["cut"]], table(x$from, x$to, kind_cut))
  if ("rate_from" %in% names(x)) {
    out <- c(out, "", titles[["rate"]],
             table(x$rate_from, x$rate_to, "percent"))
  }
  if ("surv_from" %in% names(x)) {
    out <- c(out, "", titles[["surv"]],
             table(x$surv_from, x$surv_to, "value"))
  }
  out <- c(out, "", titles[["effect"]],
           table(x$effect_from, x$effect_to, "value"))
  k <- length(conds) - length(shown)
  if (k > 0L) {
    out <- c(out, sprintf("... %d more value%s of %s (see the data frame)", k,
                          if (k == 1L) "" else "s", cond_col))
  }
  out
}
