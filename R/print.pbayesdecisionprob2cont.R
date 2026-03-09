#' Print Method for pbayesdecisionprob2cont Objects
#'
#' Displays a formatted summary of Go/NoGo/Gray decision probabilities for
#' two-continuous-endpoint results returned by
#' \code{\link{pbayesdecisionprob2cont}}.
#'
#' @param x An object of class \code{pbayesdecisionprob2cont}.
#' @param digits A positive integer specifying the number of decimal places
#'        for probability values. Default is 4.
#' @param ... Further arguments passed to or from other methods (ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.pbayesdecisionprob2cont <- function(x, digits = 4, ...) {

  # Helper: format a scalar/vector/NULL as string
  fmt <- function(v) {
    if (is.null(v))    return("NULL")
    if (length(v) > 1) return(paste0("(", paste(v, collapse = ", "), ")"))
    as.character(v)
  }

  # Helper: format a 2x2 matrix as "[r1c1, r1c2; r2c1, r2c2]"
  fmt_mat <- function(m) {
    if (is.null(m)) return("NULL")
    sprintf("[%s, %s; %s, %s]", m[1,1], m[1,2], m[2,1], m[2,2])
  }

  # Extract metadata from attributes
  prob          <- attr(x, "prob")
  design        <- attr(x, "design")
  prior         <- attr(x, "prior")
  nsim          <- attr(x, "nsim")
  nMC           <- attr(x, "nMC")
  method        <- attr(x, "method")
  GoRegions     <- attr(x, "GoRegions")
  NoGoRegions   <- attr(x, "NoGoRegions")
  gamma_go      <- attr(x, "gamma_go")
  gamma_nogo    <- attr(x, "gamma_nogo")
  theta_TV1     <- attr(x, "theta_TV1")
  theta_MAV1    <- attr(x, "theta_MAV1")
  theta_TV2     <- attr(x, "theta_TV2")
  theta_MAV2    <- attr(x, "theta_MAV2")
  theta_NULL1   <- attr(x, "theta_NULL1")
  theta_NULL2   <- attr(x, "theta_NULL2")
  n_t           <- attr(x, "n_t")
  n_c           <- attr(x, "n_c")
  m_t           <- attr(x, "m_t")
  m_c           <- attr(x, "m_c")
  Sigma_t       <- attr(x, "Sigma_t")
  Sigma_c       <- attr(x, "Sigma_c")
  kappa0_t      <- attr(x, "kappa0_t")
  nu0_t         <- attr(x, "nu0_t")
  mu0_t         <- attr(x, "mu0_t")
  Lambda0_t     <- attr(x, "Lambda0_t")
  kappa0_c      <- attr(x, "kappa0_c")
  nu0_c         <- attr(x, "nu0_c")
  mu0_c         <- attr(x, "mu0_c")
  Lambda0_c     <- attr(x, "Lambda0_c")
  r             <- attr(x, "r")
  ne_t          <- attr(x, "ne_t")
  ne_c          <- attr(x, "ne_c")
  alpha0e_t     <- attr(x, "alpha0e_t")
  alpha0e_c     <- attr(x, "alpha0e_c")
  bar_ye_t      <- attr(x, "bar_ye_t")
  bar_ye_c      <- attr(x, "bar_ye_c")
  se_t          <- attr(x, "se_t")
  se_c          <- attr(x, "se_c")
  error_if_Miss <- attr(x, "error_if_Miss")
  Gray_inc_Miss <- attr(x, "Gray_inc_Miss")
  seed          <- attr(x, "seed")

  # Build info lines with fixed label width (lw) for consistent alignment
  lw  <- 17L   # label field width
  pad <- "  "  # left margin

  lines <- character(0)
  lines <- c(lines, sprintf("%s%-*s: %s",        pad, lw, "Probability type", prob))
  lines <- c(lines, sprintf("%s%-*s: %s",        pad, lw, "Design",           design))
  lines <- c(lines, sprintf("%s%-*s: %s",        pad, lw, "Prior",            prior))
  lines <- c(lines, sprintf("%s%-*s: %s",        pad, lw, "Method",           fmt(method)))
  lines <- c(lines, sprintf("%s%-*s: nsim = %s", pad, lw, "Simulations",      fmt(nsim)))
  lines <- c(lines, sprintf("%s%-*s: nMC = %s",  pad, lw, "MC draws",         fmt(nMC)))
  lines <- c(lines, sprintf("%s%-*s: %s",        pad, lw, "Seed",             fmt(seed)))

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
  lines <- c(lines, sprintf("%s%-*s: Sigma_t = %s",
                            pad, lw, "True cov (treat.)", fmt_mat(Sigma_t)))
  if (design != "uncontrolled") {
    lines <- c(lines, sprintf("%s%-*s: Sigma_c = %s",
                              pad, lw, "True cov (cont.) ", fmt_mat(Sigma_c)))
  }

  # N-Inv-Wishart prior: treatment and control each split across two lines
  if (prior == "N-Inv-Wishart") {
    lines <- c(lines, sprintf(
      "%s%-*s: kappa0_t = %s, nu0_t = %s, mu0_t = %s",
      pad, lw, "Prior (treatment)", fmt(kappa0_t), fmt(nu0_t), fmt(mu0_t)))
    lines <- c(lines, sprintf("%s%-*s  Lambda0_t = %s",
                              pad, lw, "", fmt_mat(Lambda0_t)))
    if (design %in% c("controlled", "external")) {
      lines <- c(lines, sprintf(
        "%s%-*s: kappa0_c = %s, nu0_c = %s, mu0_c = %s",
        pad, lw, "Prior (control)  ", fmt(kappa0_c), fmt(nu0_c), fmt(mu0_c)))
      lines <- c(lines, sprintf("%s%-*s  Lambda0_c = %s",
                                pad, lw, "", fmt_mat(Lambda0_c)))
    }
  }

  if (design == "uncontrolled") {
    lines <- c(lines, sprintf("%s%-*s: mu0_c = %s, r = %s",
                              pad, lw, "Hyp. control", fmt(mu0_c), fmt(r)))
  }
  if (prob == "predictive") {
    lines <- c(lines, sprintf("%s%-*s: m_t = %s, m_c = %s",
                              pad, lw, "Future trial", fmt(m_t), fmt(m_c)))
  }
  if (design == "external") {
    # External data: treatment and control on separate lines
    lines <- c(lines, sprintf(
      "%s%-*s: ne_t = %s, alpha0e_t = %s",
      pad, lw, "External (treat.)", fmt(ne_t), fmt(alpha0e_t)))
    if (!is.null(bar_ye_t)) {
      lines <- c(lines, sprintf("%s%-*s  bar_ye_t = %s, se_t = %s",
                                pad, lw, "", fmt(bar_ye_t), fmt(se_t)))
    }
    lines <- c(lines, sprintf(
      "%s%-*s: ne_c = %s, alpha0e_c = %s",
      pad, lw, "External (cont.) ", fmt(ne_c), fmt(alpha0e_c)))
    if (!is.null(bar_ye_c)) {
      lines <- c(lines, sprintf("%s%-*s  bar_ye_c = %s, se_c = %s",
                                pad, lw, "", fmt(bar_ye_c), fmt(se_c)))
    }
  }
  lines <- c(lines, sprintf("%s%-*s: error_if_Miss = %s, Gray_inc_Miss = %s",
                            pad, lw, "Miss handling",
                            fmt(error_if_Miss), fmt(Gray_inc_Miss)))

  # Determine separator width dynamically from the longest line
  title     <- "Go/NoGo/Gray Decision Probabilities (Two Continuous Endpoints)"
  sep_width <- max(nchar(title), max(nchar(lines)))
  sep       <- strrep("-", sep_width)

  # Print header block
  cat(title, "\n")
  cat(sep, "\n")
  for (ln in lines) cat(ln, "\n")
  cat(sep, "\n")

  # Format probability columns only (not scenario columns)
  scenario_cols <- c("mu_t1", "mu_t2", "mu_c1", "mu_c2")
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
