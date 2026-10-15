#' Decision boundary data for one binary endpoint
#'
#' @param design An object of class \code{"qdm_design_binary"}.
#' @param n_t Numbers of subjects in the treatment group (a single number
#'   with a concurrent control group).
#'
#' @return A list with elements \code{grid} (data frame), \code{p_go},
#'   \code{p_nogo} and \code{axes} (see \code{qdm_boundary()}).
#' @noRd
boundary_binary <- function(design, n_t) {
  check_numeric(n_t, "n_t", len = NULL, lower = 1, integer = TRUE)
  n_t <- sort(unique(as.integer(round(n_t))))
  if (design$control == "concurrent") {
    if (length(n_t) != 1L) {
      stop(paste("'n_t' must be a single number for a design with a",
                 "concurrent control group."), call. = FALSE)
    }
    n_c <- design$n_c
    grid <- expand.grid(y_t = 0:n_t, y_c = 0:n_c)
    grid <- data.frame(n_t = n_t, y_t = grid$y_t, n_c = n_c, y_c = grid$y_c)
    post_c <- posterior_binary(design$prior_c, n_c, grid$y_c,
                               external = design$external_c)
    axes <- list(x = "y_c", y = "y_t",
                 xlab = sprintf("Responders in the control group (n = %d)",
                                n_c),
                 ylab = sprintf("Responders in the treatment group (n = %d)",
                                n_t),
                 caption = NULL)
  } else {
    grid <- data.frame(n_t = rep(n_t, n_t + 1L),
                       y_t = unlist(lapply(n_t, function(n) 0:n)))
    post_c <- if (design$control == "fixed") {
      NULL
    } else {
      posterior_binary(design$prior_c, 0, 0, external = design$external_c,
                       hypothetical = design$hypothetical_c)
    }
    axes <- list(x = "n_t", y = "y_t", xlab = "Treatment group size",
                 ylab = "Responders in the treatment group", caption = NULL)
  }
  post_t <- posterior_binary(design$prior_t, grid$n_t, grid$y_t,
                             external = design$external_t)
  p <- region_prob_binary(design$criterion, post_t$a, post_t$b, post_c$a,
                          post_c$b, fixed = design$fixed_c)
  list(grid = grid, p_go = p$go, p_nogo = p$nogo, axes = axes)
}
