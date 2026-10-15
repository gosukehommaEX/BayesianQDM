#' Decision boundary data for two continuous endpoints
#'
#' @param design An object of class \code{"qdm_design_continuous2"}.
#' @param n_t Number of subjects in the treatment group (a single number).
#' @param mean_t List of two numeric vectors of observed treatment means.
#' @param cov_t Observed treatment covariance matrix.
#' @param data_c Observed control data (concurrent control only).
#'
#' @return A list with elements \code{grid} (data frame), \code{p_go},
#'   \code{p_nogo} (vectors, or lists of two for the marginal rules) and
#'   \code{axes} (see \code{qdm_boundary()}).
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
  grid <- expand.grid(mean_t1 = sort(unique(mean_t[[1L]])),
                      mean_t2 = sort(unique(mean_t[[2L]])))
  n_t <- as.integer(round(n_t))
  post_t <- posterior_continuous2(design$prior_t, n_t,
                                  cbind(grid$mean_t1, grid$mean_t2), cov_t,
                                  external = design$external_t)
  post_c <- if (concurrent) {
    control_posterior_continuous2(design, post_t, data_c$n, data_c$mean,
                                  data_c$cov)
  } else {
    control_posterior_continuous2(design, post_t)
  }
  crit <- design$criterion
  dist_c <- NULL
  if (!is.null(post_c)) {
    check_df_mm(post_t$nu[1L] - 1, post_c$nu[1L] - 1, design$method,
                "method = \"MC\"")
    dist_c <- mean_dist_continuous2(post_c, crit$type, crit$m_c)
  }
  dist_t <- mean_dist_continuous2(post_t, crit$type, crit$m_t)
  sets <- region_sets(crit)
  p <- region_prob_continuous2(crit, dist_t, dist_c, c(sets$go, sets$nogo),
                               design$method, design$nMC,
                               fixed = design$fixed_c)
  s <- split_sets(p, length(sets$go))
  grid <- data.frame(n_t = n_t, mean_t1 = grid$mean_t1,
                     mean_t2 = grid$mean_t2)
  if (concurrent) {
    grid$diff1 <- grid$mean_t1 - data_c$mean[1L]
    grid$diff2 <- grid$mean_t2 - data_c$mean[2L]
    axes <- list(x = "diff1", y = "diff2",
                 xlab = "Observed difference in means, endpoint 1",
                 ylab = "Observed difference in means, endpoint 2",
                 caption = sprintf("Control data: %s",
                                   format(data_c, digits = 4L)))
  } else {
    axes <- list(x = "mean_t1", y = "mean_t2",
                 xlab = "Observed treatment mean, endpoint 1",
                 ylab = "Observed treatment mean, endpoint 2",
                 caption = NULL)
  }
  list(grid = grid, p_go = s$go, p_nogo = s$nogo, axes = axes)
}
