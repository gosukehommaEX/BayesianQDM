#' Decision boundary data for one continuous endpoint
#'
#' @param design An object of class \code{"qdm_design_continuous"}.
#' @param n_t,mean_t,sd_t Values of the grid (at most one of \code{n_t} and
#'   \code{sd_t} with more than one value).
#' @param data_c Observed control data (concurrent control only).
#'
#' @return A list with elements \code{grid} (data frame), \code{p_go},
#'   \code{p_nogo} and \code{axes} (see \code{qdm_boundary()}).
#' @noRd
boundary_continuous <- function(design, n_t, mean_t, sd_t, data_c) {
  if (is.null(mean_t) || is.null(sd_t)) {
    stop("'mean_t' and 'sd_t' must be given for a continuous endpoint.",
         call. = FALSE)
  }
  check_numeric(n_t, "n_t", len = NULL, lower = 2, integer = TRUE)
  check_numeric(mean_t, "mean_t", len = NULL)
  check_numeric(sd_t, "sd_t", len = NULL, lower = 0, lower_open = TRUE)
  n_t <- sort(unique(as.integer(round(n_t))))
  mean_t <- sort(unique(mean_t))
  sd_t <- sort(unique(sd_t))
  if (length(n_t) > 1L && length(sd_t) > 1L) {
    stop("Only one of 'n_t' and 'sd_t' may have more than one value.",
         call. = FALSE)
  }
  concurrent <- design$control == "concurrent"
  if (concurrent && !inherits(data_c, "qdm_data_continuous")) {
    stop("'data_c' must be created by data_continuous().", call. = FALSE)
  }
  grid <- expand.grid(mean_t = mean_t, sd_t = sd_t, n_t = n_t,
                      KEEP.OUT.ATTRS = FALSE)
  grid <- grid[, c("n_t", "mean_t", "sd_t")]
  p_go <- numeric(nrow(grid))
  p_nogo <- numeric(nrow(grid))
  for (n in n_t) {
    rows <- which(grid$n_t == n)
    post_t <- posterior_continuous(design$prior_t, n, grid$mean_t[rows],
                                   grid$sd_t[rows],
                                   external = design$external_t)
    post_c <- if (concurrent) {
      control_posterior_continuous(design, post_t, data_c$n, data_c$mean,
                                   data_c$sd)
    } else {
      control_posterior_continuous(design, post_t)
    }
    if (!is.null(post_c)) {
      check_df_mm(post_t$nu[1L], post_c$nu[1L], design$method)
    }
    p <- region_prob_continuous(design$criterion, post_t, post_c,
                                design$method, design$nMC,
                                fixed = design$fixed_c)
    p_go[rows] <- p$go
    p_nogo[rows] <- p$nogo
  }
  x <- if (length(sd_t) > 1L) "sd_t" else "n_t"
  xlab <- if (x == "sd_t") {
    "Observed treatment standard deviation"
  } else {
    "Treatment group size"
  }
  if (concurrent) {
    grid$diff <- grid$mean_t - data_c$mean
    axes <- list(x = x, y = "diff", xlab = xlab,
                 ylab = "Observed difference in means",
                 caption = sprintf("Control data: %s",
                                   format(data_c, digits = 4L)))
  } else {
    axes <- list(x = x, y = "mean_t", xlab = xlab,
                 ylab = "Observed treatment mean", caption = NULL)
  }
  list(grid = grid, p_go = p_go, p_nogo = p_nogo, axes = axes)
}
