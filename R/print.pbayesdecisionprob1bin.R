#' Print Method for pbayesdecisionprob1bin Objects
#'
#' Displays a formatted summary of Go/NoGo/Gray decision probabilities
#' for binary endpoint results returned by \code{\link{pbayesdecisionprob1bin}}.
#'
#' @param x An object of class \code{pbayesdecisionprob1bin}.
#' @param digits A positive integer specifying the number of decimal places
#'        for probability values. Default is 4.
#' @param ... Further arguments passed to or from other methods (ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.pbayesdecisionprob1bin <- function(x, digits = 4, ...) {
  # Helper: format value as string (NULL -> "NULL")
  fmt <- function(v) if (is.null(v)) "NULL" else as.character(v)

  # Extract metadata from attributes
  prob          <- attr(x, "prob")
  design        <- attr(x, "design")
  gamma_go      <- attr(x, "gamma_go")
  gamma_nogo    <- attr(x, "gamma_nogo")
  n_t           <- attr(x, "n_t")
  n_c           <- attr(x, "n_c")
  a_t           <- attr(x, "a_t")
  a_c           <- attr(x, "a_c")
  b_t           <- attr(x, "b_t")
  b_c           <- attr(x, "b_c")
  z             <- attr(x, "z")
  m_t           <- attr(x, "m_t")
  m_c           <- attr(x, "m_c")
  ne_t          <- attr(x, "ne_t")
  ne_c          <- attr(x, "ne_c")
  ye_t          <- attr(x, "ye_t")
  ye_c          <- attr(x, "ye_c")
  alpha0e_t          <- attr(x, "alpha0e_t")
  alpha0e_c          <- attr(x, "alpha0e_c")
  error_if_Miss <- attr(x, "error_if_Miss")
  Gray_inc_Miss <- attr(x, "Gray_inc_Miss")

  # Build threshold string based on probability type
  if (prob == "posterior") {
    theta_str <- sprintf("TV = %s, MAV = %s",
                         fmt(attr(x, "theta_TV")), fmt(attr(x, "theta_MAV")))
  } else {
    theta_str <- sprintf("NULL = %s", fmt(attr(x, "theta_NULL")))
  }

  # Prior label depends on probability type
  prior_label <- if (prob == "posterior") "Prior (Beta)    " else "Prior (Beta-bin)"

  # Build info lines with fixed label width (lw) for consistent alignment
  lw  <- 17L   # label field width
  pad <- "  "  # left margin

  lines <- character(0)
  lines <- c(lines, sprintf("%s%-*s: %s", pad, lw, "Probability type", prob))
  lines <- c(lines, sprintf("%s%-*s: %s", pad, lw, "Design",           design))
  lines <- c(lines, sprintf("%s%-*s: %s", pad, lw, "Threshold(s)",     theta_str))
  lines <- c(lines, sprintf("%s%-*s: gamma_go = %s",
                            pad, lw, "Go  threshold",  fmt(gamma_go)))
  lines <- c(lines, sprintf("%s%-*s: gamma_nogo = %s",
                            pad, lw, "NoGo threshold", fmt(gamma_nogo)))
  lines <- c(lines, sprintf("%s%-*s: n_t = %s, n_c = %s",
                            pad, lw, "Sample size", fmt(n_t), fmt(n_c)))
  lines <- c(lines, sprintf("%s%-*s: a_t = %s, a_c = %s, b_t = %s, b_c = %s",
                            pad, lw, prior_label,
                            fmt(a_t), fmt(a_c), fmt(b_t), fmt(b_c)))
  if (design == "uncontrolled") {
    lines <- c(lines, sprintf("%s%-*s: z = %s",
                              pad, lw, "Uncontrolled", fmt(z)))
  }
  if (prob == "predictive") {
    lines <- c(lines, sprintf("%s%-*s: m_t = %s, m_c = %s",
                              pad, lw, "Future trial", fmt(m_t), fmt(m_c)))
  }
  if (design == "external") {
    # Split external parameters across two lines to avoid overflow
    lines <- c(lines, sprintf("%s%-*s: ne_t = %s, ne_c = %s, ye_t = %s, ye_c = %s",
                              pad, lw, "External data",
                              fmt(ne_t), fmt(ne_c), fmt(ye_t), fmt(ye_c)))
    lines <- c(lines, sprintf("%s%-*s  alpha0e_t = %s, alpha0e_c = %s",
                              pad, lw, "", fmt(alpha0e_t), fmt(alpha0e_c)))
  }
  lines <- c(lines, sprintf("%s%-*s: error_if_Miss = %s, Gray_inc_Miss = %s",
                            pad, lw, "Miss handling",
                            fmt(error_if_Miss), fmt(Gray_inc_Miss)))

  # Determine separator width dynamically from the longest line
  title     <- "Go/NoGo/Gray Decision Probabilities (Single Binary Endpoint)"
  sep_width <- max(nchar(title), max(nchar(lines)))
  sep       <- strrep("-", sep_width)

  # Print header block
  cat(title, "\n")
  cat(sep, "\n")
  for (ln in lines) cat(ln, "\n")
  cat(sep, "\n")

  # Format probability columns only (not pi_t / pi_c)
  prob_cols <- names(x)[!names(x) %in% c("pi_t", "pi_c")]
  x_print   <- x
  x_print[prob_cols] <- lapply(x[prob_cols], function(col) {
    formatC(col, digits = digits, format = "f")
  })

  # Print table without row names (explicit call to avoid recursion)
  print.data.frame(x_print, row.names = FALSE, quote = FALSE)
  cat(sep, "\n")

  invisible(x)
}
