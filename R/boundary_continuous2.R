#' Decision boundary specification for two continuous endpoints
#'
#' The conditions are the observed treatment means of endpoint 1
#' (\code{mean_t[[1]]}), and the cut variable is the observed treatment mean
#' of endpoint 2, searched on the grid \code{mean_t[[2]]} (with two values
#' far beyond each end added).
#'
#' @param design An object of class \code{"qdm_design_continuous2"}.
#' @param n_t Number of subjects in the treatment group.
#' @param mean_t List of the values of endpoint 1 and the grid of endpoint 2.
#' @param cov_t Observed treatment covariance matrix.
#' @param data_c Observed control data (concurrent control only).
#'
#' @return A list with the elements used by \code{qdm_boundary()} (see
#'   \code{boundary_binary()}) and \code{range}, the range of the grid.
#' @noRd
boundary_continuous2 <- function(design, n_t, mean_t, cov_t, data_c) {
  if (!is.list(mean_t) || length(mean_t) != 2L || is.null(cov_t)) {
    stop(paste("For two continuous endpoints, 'mean_t' must be a list of two",
               "numeric vectors and 'cov_t' must be given."), call. = FALSE)
  }
  check_numeric(n_t, "n_t", lower = 3, integer = TRUE)
  check_numeric(mean_t[[1L]], "mean_t[[1]]", len = NULL)
  check_numeric(mean_t[[2L]], "mean_t[[2]]", len = NULL)
  check_pd_matrix(cov_t, "cov_t")
  concurrent <- design$control == "concurrent"
  if (concurrent && !inherits(data_c, "qdm_data_continuous2")) {
    stop("'data_c' must be created by data_continuous2().", call. = FALSE)
  }
  n_t <- as.integer(round(n_t))
  crit <- design$criterion
  control_post <- function(post_t) {
    if (concurrent) {
      control_posterior_continuous2(design, post_t, data_c$n, data_c$mean,
                                    data_c$cov)
    } else {
      control_posterior_continuous2(design, post_t)
    }
  }
  post0 <- posterior_continuous2(design$prior_t, n_t, matrix(0, 1L, 2L),
                                 cov_t, external = design$external_t)
  post_c0 <- control_post(post0)
  if (!is.null(post_c0)) {
    check_df_mm(post0$nu[1L] - 1, post_c0$nu[1L] - 1, design$method,
                "method = \"MC\"")
  }
  sets <- region_sets(crit)
  cond <- data.frame(n_t = n_t, mean_t1 = sort(unique(mean_t[[1L]])))
  prob <- function(idx, m2) {
    post_t <- posterior_continuous2(design$prior_t, n_t,
                                    cbind(cond$mean_t1[idx], m2), cov_t,
                                    external = design$external_t)
    post_c <- control_post(post_t)
    dist_c <- if (is.null(post_c)) {
      NULL
    } else {
      mean_dist_continuous2(post_c, crit$type, crit$m_c)
    }
    dist_t <- mean_dist_continuous2(post_t, crit$type, crit$m_t)
    p <- region_prob_continuous2(crit, dist_t, dist_c, c(sets$go, sets$nogo),
                                 design$method, design$nMC,
                                 fixed = design$fixed_c)
    split_sets(p, length(sets$go))
  }
  # Control means used for the observed differences
  control <- if (concurrent) {
    data_c$mean
  } else if (design$control == "fixed") {
    design$fixed_c
  } else if (design$control == "hypothetical") {
    design$hypothetical_c$mean
  } else {
    as.numeric(post_c0$mu[1L, ])
  }
  cond$effect1 <- cond$mean_t1 - control[1L]
  grid <- sort(unique(mean_t[[2L]]))
  span <- if (length(grid) > 1L) diff(range(grid)) else max(1, abs(grid))
  values <- sort(c(grid[1L] - c(1e4, 1e2) * span, grid,
                   grid[length(grid)] + c(1e2, 1e4) * span))
  ref <- format_num(control[2L], 4L)
  what <- switch(design$control,
                 concurrent = "control mean",
                 fixed = "fixed control mean",
                 hypothetical = "location of the hypothetical control",
                 external = "posterior mean of the control mean")
  list(cond = cond, values = rep(list(values), nrow(cond)), prob = prob,
       discrete = FALSE,
       derived = function(idx, m2) list(effect = m2 - control[2L]),
       axes = list(x = if (concurrent) "effect1" else "mean_t1",
                   y = if (concurrent) "effect" else "cut",
                   offset = if (concurrent) control[2L] else 0,
                   discrete_x = FALSE, discrete_y = FALSE,
                   xlab = if (concurrent) {
                     "Observed difference in means, endpoint 1"
                   } else {
                     "Observed treatment mean, endpoint 1"
                   },
                   ylab = if (concurrent) {
                     "Observed difference in means, endpoint 2"
                   } else {
                     "Observed treatment mean, endpoint 2"
                   },
                   sec_scale = NULL,
                   caption = if (concurrent) {
                     sprintf("Control data: %s", format(data_c, digits = 4L))
                   }),
       titles = c(cut = "Observed treatment mean of endpoint 2",
                  effect = sprintf(paste("Observed difference in means of",
                                         "endpoint 2 (%s %s)"), what, ref)),
       increasing = crit$direction[2L] == "greater", range = range(grid))
}
