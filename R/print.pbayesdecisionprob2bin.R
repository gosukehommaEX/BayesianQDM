#' Print Method for pbayesdecisionprob2bin Objects
#'
#' Displays a formatted summary of Go/NoGo/Gray decision probabilities for
#' two-binary-endpoint results returned by \code{\link{pbayesdecisionprob2bin}}.
#'
#' @param x An object of class \code{pbayesdecisionprob2bin}.
#' @param digits A positive integer specifying the number of decimal places
#'        for probability values.  Default is 4.
#' @param ... Further arguments passed to or from other methods (ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.pbayesdecisionprob2bin <- function(x, digits = 4, ...) {

  # Helper: format a value as string (NULL -> "NULL")
  fmt <- function(v) if (is.null(v)) "NULL" else as.character(v)

  # Extract metadata
  prob          <- attr(x, "prob")
  design        <- attr(x, "design")
  GoRegions     <- attr(x, "GoRegions")
  NoGoRegions   <- attr(x, "NoGoRegions")
  gamma_go      <- attr(x, "gamma_go")
  gamma_nogo    <- attr(x, "gamma_nogo")
  n_t           <- attr(x, "n_t")
  n_c           <- attr(x, "n_c")
  a_t_00        <- attr(x, "a_t_00")
  a_t_01        <- attr(x, "a_t_01")
  a_t_10        <- attr(x, "a_t_10")
  a_t_11        <- attr(x, "a_t_11")
  a_c_00        <- attr(x, "a_c_00")
  a_c_01        <- attr(x, "a_c_01")
  a_c_10        <- attr(x, "a_c_10")
  a_c_11        <- attr(x, "a_c_11")
  m_t           <- attr(x, "m_t")
  m_c           <- attr(x, "m_c")
  theta_TV1     <- attr(x, "theta_TV1")
  theta_MAV1    <- attr(x, "theta_MAV1")
  theta_TV2     <- attr(x, "theta_TV2")
  theta_MAV2    <- attr(x, "theta_MAV2")
  theta_NULL1   <- attr(x, "theta_NULL1")
  theta_NULL2   <- attr(x, "theta_NULL2")
  z00           <- attr(x, "z00")
  z01           <- attr(x, "z01")
  z10           <- attr(x, "z10")
  z11           <- attr(x, "z11")
  xe_t_00       <- attr(x, "xe_t_00")
  xe_t_01       <- attr(x, "xe_t_01")
  xe_t_10       <- attr(x, "xe_t_10")
  xe_t_11       <- attr(x, "xe_t_11")
  xe_c_00       <- attr(x, "xe_c_00")
  xe_c_01       <- attr(x, "xe_c_01")
  xe_c_10       <- attr(x, "xe_c_10")
  xe_c_11       <- attr(x, "xe_c_11")
  alpha0e_t     <- attr(x, "alpha0e_t")
  alpha0e_c     <- attr(x, "alpha0e_c")
  nMC           <- attr(x, "nMC")
  nsim          <- attr(x, "nsim")
  CalcMethod    <- attr(x, "CalcMethod")
  error_if_Miss <- attr(x, "error_if_Miss")
  Gray_inc_Miss <- attr(x, "Gray_inc_Miss")

  # Build info lines with fixed label width (lw) for consistent alignment
  lw  <- 17L   # label field width
  pad <- "  "  # left margin

  lines <- character(0)
  lines <- c(lines, sprintf("%s%-*s: %s",           pad, lw, "Probability type", prob))
  lines <- c(lines, sprintf("%s%-*s: %s",           pad, lw, "Design",           design))

  # Threshold(s): split posterior across two lines to avoid overflow
  if (prob == "posterior") {
    lines <- c(lines, sprintf("%s%-*s: TV1 = %s, MAV1 = %s",
                              pad, lw, "Threshold(s)",
                              fmt(theta_TV1), fmt(theta_MAV1)))
    lines <- c(lines, sprintf("%s%-*s  TV2 = %s, MAV2 = %s",
                              pad, lw, "",
                              fmt(theta_TV2), fmt(theta_MAV2)))
  } else {
    lines <- c(lines, sprintf("%s%-*s: NULL1 = %s, NULL2 = %s",
                              pad, lw, "Threshold(s)",
                              fmt(theta_NULL1), fmt(theta_NULL2)))
  }

  lines <- c(lines, sprintf("%s%-*s: gamma_go = %s",
                            pad, lw, "Go  threshold",  fmt(gamma_go)))
  lines <- c(lines, sprintf("%s%-*s: gamma_nogo = %s",
                            pad, lw, "NoGo threshold", fmt(gamma_nogo)))
  lines <- c(lines, sprintf("%s%-*s: {%s}",
                            pad, lw, "Go  regions",
                            paste(GoRegions,   collapse = ", ")))
  lines <- c(lines, sprintf("%s%-*s: {%s}",
                            pad, lw, "NoGo regions",
                            paste(NoGoRegions, collapse = ", ")))
  lines <- c(lines, sprintf("%s%-*s: n_t = %s, n_c = %s",
                            pad, lw, "Sample size", fmt(n_t), fmt(n_c)))

  # Dirichlet prior: treatment and control on separate lines
  lines <- c(lines, sprintf(
    "%s%-*s: a_t = (%s, %s, %s, %s)  [a_00, a_01, a_10, a_11]",
    pad, lw, "Prior (treatment)",
    fmt(a_t_00), fmt(a_t_01), fmt(a_t_10), fmt(a_t_11)))
  lines <- c(lines, sprintf(
    "%s%-*s: a_c = (%s, %s, %s, %s)  [a_00, a_01, a_10, a_11]",
    pad, lw, "Prior (control)  ",
    fmt(a_c_00), fmt(a_c_01), fmt(a_c_10), fmt(a_c_11)))

  if (design == "uncontrolled") {
    lines <- c(lines, sprintf(
      "%s%-*s: z = (%s, %s, %s, %s)  [z00, z01, z10, z11]",
      pad, lw, "Hyp. control",
      fmt(z00), fmt(z01), fmt(z10), fmt(z11)))
  }
  if (prob == "predictive") {
    lines <- c(lines, sprintf("%s%-*s: m_t = %s, m_c = %s",
                              pad, lw, "Future trial", fmt(m_t), fmt(m_c)))
  }
  if (design == "external") {
    # External data: treatment and control counts on separate lines,
    # followed by power prior weights
    lines <- c(lines, sprintf(
      "%s%-*s: xe_t = (%s, %s, %s, %s)  [xe_00, xe_01, xe_10, xe_11]",
      pad, lw, "External (treat.)",
      fmt(xe_t_00), fmt(xe_t_01), fmt(xe_t_10), fmt(xe_t_11)))
    lines <- c(lines, sprintf(
      "%s%-*s: xe_c = (%s, %s, %s, %s)  [xe_00, xe_01, xe_10, xe_11]",
      pad, lw, "External (cont.) ",
      fmt(xe_c_00), fmt(xe_c_01), fmt(xe_c_10), fmt(xe_c_11)))
    lines <- c(lines, sprintf("%s%-*s: alpha0e_t = %s, alpha0e_c = %s",
                              pad, lw, "Power prior",
                              fmt(alpha0e_t), fmt(alpha0e_c)))
  }
  lines <- c(lines, sprintf("%s%-*s: %s",           pad, lw, "Method",           fmt(CalcMethod)))
  lines <- c(lines, sprintf("%s%-*s: nMC = %s",     pad, lw, "MC draws",         fmt(nMC)))
  if (CalcMethod == "MC") {
    lines <- c(lines, sprintf("%s%-*s: nsim = %s",  pad, lw, "Sim size",         fmt(nsim)))
  }
  lines <- c(lines, sprintf("%s%-*s: error_if_Miss = %s, Gray_inc_Miss = %s",
                            pad, lw, "Miss handling",
                            fmt(error_if_Miss), fmt(Gray_inc_Miss)))

  # Determine separator width dynamically from the longest line
  title     <- "Go/NoGo/Gray Decision Probabilities (Two Binary Endpoints)"
  sep_width <- max(nchar(title), max(nchar(lines)))
  sep       <- strrep("-", sep_width)

  # Print header block
  cat(title, "\n")
  cat(sep, "\n")
  for (ln in lines) cat(ln, "\n")
  cat(sep, "\n")

  # Format probability columns only (not scenario columns)
  scenario_cols <- c("pi_t1", "pi_t2", "rho_t", "pi_c1", "pi_c2", "rho_c")
  prob_cols     <- names(x)[!names(x) %in% scenario_cols]

  x_print <- x
  x_print[prob_cols] <- lapply(x[prob_cols], function(col) {
    formatC(col, digits = digits, format = "f")
  })

  # Print table without row names (explicit call to avoid recursion)
  print.data.frame(x_print, row.names = FALSE, quote = FALSE)
  cat(sep, "\n")

  invisible(x)
}
