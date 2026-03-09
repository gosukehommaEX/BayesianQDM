#' Print Method for pbayesdecisionprob1cont Objects
#'
#' Displays a formatted summary of Go/NoGo/Gray decision probabilities
#' for continuous endpoint results returned by \code{\link{pbayesdecisionprob1cont}}.
#'
#' @param x An object of class \code{pbayesdecisionprob1cont}.
#' @param digits A positive integer specifying the number of decimal places
#'        for probability values. Default is 4.
#' @param ... Further arguments passed to or from other methods (ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.pbayesdecisionprob1cont <- function(x, digits = 4, ...) {
  # Helper to format a value as string (NULL -> "NULL")
  fmt <- function(v) if (is.null(v)) "NULL" else as.character(v)

  # Extract metadata from attributes
  prob          <- attr(x, "prob")
  design        <- attr(x, "design")
  prior         <- attr(x, "prior")
  CalcMethod    <- attr(x, "CalcMethod")
  nsim          <- attr(x, "nsim")
  nMC           <- attr(x, "nMC")
  gamma_go      <- attr(x, "gamma_go")
  gamma_nogo    <- attr(x, "gamma_nogo")
  n_t           <- attr(x, "n_t")
  n_c           <- attr(x, "n_c")
  sigma_t       <- attr(x, "sigma_t")
  sigma_c       <- attr(x, "sigma_c")
  r             <- attr(x, "r")
  m_t           <- attr(x, "m_t")
  m_c           <- attr(x, "m_c")
  kappa0_t      <- attr(x, "kappa0_t")
  kappa0_c      <- attr(x, "kappa0_c")
  nu0_t         <- attr(x, "nu0_t")
  nu0_c         <- attr(x, "nu0_c")
  mu0_t         <- attr(x, "mu0_t")
  mu0_c         <- attr(x, "mu0_c")
  sigma0_t      <- attr(x, "sigma0_t")
  sigma0_c      <- attr(x, "sigma0_c")
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

  # Build threshold string based on probability type
  if (prob == "posterior") {
    theta_str <- sprintf("TV = %s, MAV = %s",
                         fmt(attr(x, "theta_TV")), fmt(attr(x, "theta_MAV")))
  } else {
    theta_str <- sprintf("NULL = %s", fmt(attr(x, "theta_NULL")))
  }

  # Build info lines with fixed label width (lw) for consistent alignment
  lw  <- 17L   # label field width
  pad <- "  "  # left margin

  lines <- character(0)
  lines <- c(lines, sprintf("%s%-*s: %s",           pad, lw, "Probability type", prob))
  lines <- c(lines, sprintf("%s%-*s: %s",           pad, lw, "Design",           design))
  lines <- c(lines, sprintf("%s%-*s: %s",           pad, lw, "Prior",            prior))
  lines <- c(lines, sprintf("%s%-*s: %s",           pad, lw, "Calc method",      CalcMethod))
  lines <- c(lines, sprintf("%s%-*s: nsim = %s",    pad, lw, "Simulations",      fmt(nsim)))
  if (!is.null(nMC)) {
    lines <- c(lines, sprintf("%s%-*s: nMC = %s",   pad, lw, "MC draws",         fmt(nMC)))
  }
  lines <- c(lines, sprintf("%s%-*s: %s",           pad, lw, "Threshold(s)",     theta_str))
  lines <- c(lines, sprintf("%s%-*s: gamma_go = %s",   pad, lw, "Go  threshold",  fmt(gamma_go)))
  lines <- c(lines, sprintf("%s%-*s: gamma_nogo = %s", pad, lw, "NoGo threshold", fmt(gamma_nogo)))
  lines <- c(lines, sprintf("%s%-*s: n_t = %s, n_c = %s",
                            pad, lw, "Sample size", fmt(n_t), fmt(n_c)))
  lines <- c(lines, sprintf("%s%-*s: sigma_t = %s, sigma_c = %s",
                            pad, lw, "True SD", fmt(sigma_t), fmt(sigma_c)))
  if (design == "uncontrolled") {
    lines <- c(lines, sprintf("%s%-*s: r = %s",     pad, lw, "Variance ratio",   fmt(r)))
  }
  if (!is.null(m_t) || !is.null(m_c)) {
    lines <- c(lines, sprintf("%s%-*s: m_t = %s, m_c = %s",
                              pad, lw, "Future size", fmt(m_t), fmt(m_c)))
  }
  if (prior == "N-Inv-Chisq") {
    # Split treatment and control prior parameters across two lines
    lines <- c(lines, sprintf("%s%-*s: kappa0_t = %s, nu0_t = %s, mu0_t = %s, sigma0_t = %s",
                              pad, lw, "Prior (treatment)",
                              fmt(kappa0_t), fmt(nu0_t), fmt(mu0_t), fmt(sigma0_t)))
    lines <- c(lines, sprintf("%s%-*s: kappa0_c = %s, nu0_c = %s, mu0_c = %s, sigma0_c = %s",
                              pad, lw, "Prior (control)  ",
                              fmt(kappa0_c), fmt(nu0_c), fmt(mu0_c), fmt(sigma0_c)))
  }
  if (design == "external") {
    # Split treatment and control external parameters across two lines each
    lines <- c(lines, sprintf("%s%-*s: ne_t = %s, alpha0e_t = %s, bar_ye_t = %s, se_t = %s",
                              pad, lw, "External (treat.)",
                              fmt(ne_t), fmt(alpha0e_t), fmt(bar_ye_t), fmt(se_t)))
    lines <- c(lines, sprintf("%s%-*s: ne_c = %s, alpha0e_c = %s, bar_ye_c = %s, se_c = %s",
                              pad, lw, "External (cont.) ",
                              fmt(ne_c), fmt(alpha0e_c), fmt(bar_ye_c), fmt(se_c)))
  }
  lines <- c(lines, sprintf("%s%-*s: error_if_Miss = %s, Gray_inc_Miss = %s",
                            pad, lw, "Miss handling",
                            fmt(error_if_Miss), fmt(Gray_inc_Miss)))
  lines <- c(lines, sprintf("%s%-*s: %s",           pad, lw, "Seed",             fmt(seed)))

  # Determine separator width dynamically from the longest line
  title     <- "Go/NoGo/Gray Decision Probabilities (Single Continuous Endpoint)"
  sep_width <- max(nchar(title), max(nchar(lines)))
  sep       <- strrep("-", sep_width)

  # Print header block
  cat(title, "\n")
  cat(sep, "\n")
  for (ln in lines) cat(ln, "\n")
  cat(sep, "\n")

  # Format probability columns only (not mu_t / mu_c)
  df_print  <- as.data.frame(x)
  prob_cols <- names(df_print)[!names(df_print) %in% c("mu_t", "mu_c")]
  df_print[prob_cols] <- lapply(df_print[prob_cols],
                                function(col) round(col, digits))

  # Print table without row names (explicit call to avoid recursion)
  print.data.frame(df_print, row.names = FALSE)
  cat(sep, "\n")

  invisible(x)
}
