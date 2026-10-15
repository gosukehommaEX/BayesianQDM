#' Decision boundary specification for one continuous endpoint
#'
#' The cut variable is the observed treatment mean. Its grid is
#' \code{mean_t} or, by default, 201 equally spaced values over a range that
#' covers the thresholds of the criterion shifted by the location of the
#' control distribution, widened by ten times the scale of the difference of
#' the means (and mapped to the observed mean through the posterior location
#' of the treatment mean). Two values far beyond each end of the grid are
#' added, so that decision changes just outside the grid are found.
#'
#' @param design An object of class \code{"qdm_design_continuous"}.
#' @param n_t,sd_t Values of the conditions (at most one of them with more
#'   than one value).
#' @param mean_t Grid of the observed treatment mean, or \code{NULL}.
#' @param data_c Observed control data (concurrent control only).
#'
#' @return A list with the elements used by \code{qdm_boundary()} (see
#'   \code{boundary_binary()}) and \code{range}, the range of the grid.
#' @noRd
boundary_continuous <- function(design, n_t, mean_t, sd_t, data_c) {
  if (is.null(sd_t)) {
    stop("'sd_t' must be given for a continuous endpoint.", call. = FALSE)
  }
  check_numeric(n_t, "n_t", len = NULL, lower = 2, integer = TRUE)
  check_numeric(sd_t, "sd_t", len = NULL, lower = 0, lower_open = TRUE)
  if (!is.null(mean_t)) {
    check_numeric(mean_t, "mean_t", len = NULL)
  }
  n_t <- sort(unique(as.integer(round(n_t))))
  sd_t <- sort(unique(sd_t))
  if (length(n_t) > 1L && length(sd_t) > 1L) {
    stop("Only one of 'n_t' and 'sd_t' may have more than one value.",
         call. = FALSE)
  }
  concurrent <- design$control == "concurrent"
  if (concurrent && !inherits(data_c, "qdm_data_continuous")) {
    stop("'data_c' must be created by data_continuous().", call. = FALSE)
  }
  crit <- design$criterion
  cond <- expand.grid(n_t = n_t, sd_t = sd_t, KEEP.OUT.ATTRS = FALSE)
  control_post <- function(post_t) {
    if (concurrent) {
      control_posterior_continuous(design, post_t, data_c$n, data_c$mean,
                                   data_c$sd)
    } else {
      control_posterior_continuous(design, post_t)
    }
  }
  # Degrees of freedom do not depend on the observed mean and SD
  for (n in n_t) {
    post_t <- posterior_continuous(design$prior_t, n, 0, sd_t[1L],
                                   external = design$external_t)
    post_c <- control_post(post_t)
    if (!is.null(post_c)) {
      check_df_mm(post_t$nu, post_c$nu[1L], design$method)
    }
  }
  prob <- function(idx, m) {
    go <- numeric(length(idx))
    nogo <- numeric(length(idx))
    for (n in unique(cond$n_t[idx])) {
      rows <- which(cond$n_t[idx] == n)
      post_t <- posterior_continuous(design$prior_t, n, m[rows],
                                     cond$sd_t[idx[rows]],
                                     external = design$external_t)
      p <- region_prob_continuous(crit, post_t, control_post(post_t),
                                  design$method, design$nMC,
                                  fixed = design$fixed_c)
      go[rows] <- p$go
      nogo[rows] <- p$nogo
    }
    list(go = go, nogo = nogo)
  }
  # Location of the control distribution (the observed mean for a
  # concurrent control, used for the observed differences)
  if (concurrent) {
    control <- data_c$mean
    loc_c <- control_post(NULL)$mu
  } else {
    control <- control_location(design)
    loc_c <- control
  }
  if (is.null(mean_t)) {
    # Scale of the difference of the means
    post_t <- posterior_continuous(design$prior_t, min(n_t), 0, max(sd_t),
                                   external = design$external_t)
    post_c <- control_post(post_t)
    s2 <- post_t$sigma2 / post_t$kappa
    if (!is.null(post_c)) {
      s2 <- s2 + post_c$sigma2 / post_c$kappa
    }
    thr <- if (crit$type == "posterior") {
      c(crit$theta_tv, crit$theta_mav)
    } else {
      crit$theta_null
    }
    target <- loc_c + c(min(thr), max(thr)) + c(-10, 10) * sqrt(s2)
    # Observed mean that gives each posterior location (linear in the mean)
    ends <- unlist(lapply(n_t, function(n) {
      a <- posterior_continuous(design$prior_t, n, 0, 1,
                                external = design$external_t)$mu
      b <- posterior_continuous(design$prior_t, n, 1, 1,
                                external = design$external_t)$mu - a
      (target - a) / b
    }))
    mean_t <- seq(min(ends), max(ends), length.out = 201L)
  }
  grid <- sort(unique(mean_t))
  span <- if (length(grid) > 1L) diff(range(grid)) else max(1, abs(grid))
  values <- sort(c(grid[1L] - c(1e4, 1e2) * span, grid,
                   grid[length(grid)] + c(1e2, 1e4) * span))
  x <- if (length(sd_t) > 1L) "sd_t" else "n_t"
  ref <- format_num(control, 4L)
  effect_title <- switch(
    design$control,
    concurrent = sprintf(paste("Observed difference in means (control mean",
                               "%s)"), ref),
    fixed = sprintf("Observed difference from the fixed control mean (%s)",
                    ref),
    hypothetical = sprintf(paste("Observed difference from the location of",
                                 "the hypothetical control (%s)"), ref),
    external = sprintf(paste("Observed difference from the posterior mean",
                             "of the control mean (%s)"), ref)
  )
  list(cond = cond, values = rep(list(values), nrow(cond)), prob = prob,
       discrete = FALSE,
       derived = function(idx, m) list(effect = m - control),
       axes = list(x = x, y = if (concurrent) "effect" else "cut",
                   offset = if (concurrent) control else 0,
                   discrete_x = x == "n_t", discrete_y = FALSE,
                   xlab = if (x == "sd_t") {
                     "Observed treatment standard deviation"
                   } else {
                     "Treatment group size"
                   },
                   ylab = if (concurrent) {
                     "Observed difference in means"
                   } else {
                     "Observed treatment mean"
                   },
                   sec_scale = NULL,
                   caption = if (concurrent) {
                     sprintf("Control data: %s", format(data_c, digits = 4L))
                   }),
       titles = c(cut = "Observed treatment mean", effect = effect_title),
       increasing = crit$direction == "greater", range = range(grid))
}
