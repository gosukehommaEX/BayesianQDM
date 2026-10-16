#' Plot a decision boundary
#'
#' Draws the decision intervals computed by \code{\link{qdm_boundary}} as
#' regions filled by the decision (Go, NoGo, Gray, or Miss if it occurs):
#' the horizontal axis is the condition (for example, the number of control
#' responders or the treatment group size) and the vertical axis is the
#' observed quantity (numbers of responders or events, or the observed
#' treatment mean) or, for continuous endpoints with a concurrent control
#' group, the observed difference in means. The boundary values are
#' written on the Go and NoGo regions next to the boundaries, for example
#' the smallest number of responders that leads to a Go decision. With a
#' concurrent control group and a binary endpoint, the observed response
#' rate of the treatment group is shown on a secondary axis.
#'
#' The plot is a \code{ggplot} object, so it can be changed further by adding
#' \pkg{ggplot2} functions, for example \code{ggplot2::labs()},
#' \code{ggplot2::theme()} and \code{ggplot2::scale_fill_manual()}.
#'
#' @inheritParams plot.qdm_oc
#' @param x An object of class \code{"qdm_boundary"}.
#' @param labels \code{TRUE} (the default) to write the boundary values for
#'   every value of the horizontal axis (or, with more than 16 values, for
#'   about six of them, leaving out the first and the last), \code{FALSE}
#'   for none, or a numeric vector of the values of the horizontal axis at
#'   which they are written.
#' @param label_size Size of the boundary values, in millimeters, or
#'   \code{NULL} for the size of the axis tick labels,
#'   \code{0.8 * base_size / .pt}.
#' @param ... Not used; arguments given here are ignored with a warning.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' des <- design_binary(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_beta(0.5, 0.5), prior_c = prior_beta(0.5, 0.5),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                   gamma_go = 0.5, gamma_nogo = 0.8)
#' )
#' b <- qdm_boundary(des)
#' plot(b)
#' plot(b, labels = c(0, 4, 8), col_gray = "gray80", base_size = 14,
#'      legend_position = "right")
#' @export
plot.qdm_boundary <- function(x, col_go = "#658D1B", col_nogo = "#D91E49",
                              col_gray = "#939597", col_miss = "#E69F00",
                              labels = TRUE, label_size = NULL,
                              base_size = 11, legend_position = "bottom",
                              ...) {
  axes <- attr(x, "axes")
  if (is.null(axes)) {
    stop("'x' must be an unmodified result of qdm_boundary().",
         call. = FALSE)
  }
  warn_unused_dots(list(...), "plot.qdm_boundary")
  cols <- decision_colors(col_go = col_go, col_nogo = col_nogo,
                          col_gray = col_gray, col_miss = col_miss)
  style <- plot_style(label_size = label_size, base_size = base_size,
                      legend_position = legend_position)
  xv <- x[[axes$x]]
  ux <- sort(unique(xv))
  if (isTRUE(labels)) {
    lab_x <- if (length(ux) <= 16L) {
      ux
    } else {
      # About six values near pretty values, without the first and the last
      # value, whose labels would be cut at the edges of the panel
      thin <- unique(ux[vapply(pretty(ux, n = 6L), function(v) {
        which.min(abs(ux - v))
      }, integer(1))])
      setdiff(thin, ux[c(1L, length(ux))])
    }
  } else if (isFALSE(labels)) {
    lab_x <- numeric(0)
  } else {
    found <- is.numeric(labels) && length(labels) > 0L && !anyNA(labels) &&
      all(vapply(labels, function(v) {
        any(abs(ux - v) <= 1e-8 * max(1, abs(v)))
      }, logical(1)))
    if (!found) {
      stop(paste("'labels' must be TRUE, FALSE or values of the horizontal",
                 "axis of the plot."), call. = FALSE)
    }
    lab_x <- ux[vapply(ux, function(u) {
      any(abs(labels - u) <= 1e-8 * max(1, abs(u)))
    }, logical(1))]
  }
  # Horizontal extent of the regions of each value of the condition: up to
  # the midpoints between neighboring values (half a unit for consecutive
  # whole numbers)
  if (length(ux) == 1L) {
    left <- ux - 0.5
    right <- ux + 0.5
  } else {
    mids <- (ux[-1L] + ux[-length(ux)]) / 2
    k <- length(ux)
    left <- c(ux[1L] - (mids[1L] - ux[1L]), mids)
    right <- c(mids, ux[k] + (ux[k] - mids[k - 1L]))
  }
  pos <- match(xv, ux)
  if (axes$y == "effect") {
    lo <- pmin(x$effect_from, x$effect_to)
    hi <- pmax(x$effect_from, x$effect_to)
  } else {
    lo <- x$from
    hi <- x$to
  }
  half <- if (axes$discrete_y) 0.5 else 0
  d <- data.frame(xmin = left[pos], xmax = right[pos], ymin = lo - half,
                  ymax = hi + half, decision = x$decision)
  # Boundary values next to the boundaries of the Go and NoGo regions: the
  # lower end of a region with another region below, and the upper end of
  # a region with another region above
  ends <- lapply(seq_len(nrow(x)), function(i) {
    if (!(xv[i] %in% lab_x) || !(x$decision[i] %in% c("Go", "NoGo"))) {
      return(NULL)
    }
    same <- xv == xv[i]
    same[i] <- FALSE
    below <- any(same & hi <= lo[i])
    above <- any(same & lo >= hi[i])
    rbind(if (below && is.finite(lo[i])) {
      data.frame(x = xv[i], y = lo[i], vjust = if (half > 0) 0.5 else -0.4)
    }, if (above && is.finite(hi[i]) && (hi[i] != lo[i] || !below)) {
      data.frame(x = xv[i], y = hi[i], vjust = if (half > 0) 0.5 else 1.4)
    })
  })
  lab <- do.call(rbind, ends)
  if (!is.null(lab)) {
    lab$text <- vapply(lab$y, format, character(1), digits = 3L)
  }
  # Vertical range: the finite ends of the intervals with a margin
  fin <- c(lo, hi)[is.finite(c(lo, hi))]
  if (half == 0) {
    rng <- attr(x, "range") - axes$offset
    span <- diff(rng)
    fin <- fin[fin >= rng[1L] - span & fin <= rng[2L] + span]
    if (length(fin) == 0L) {
      fin <- rng
    }
    width <- max(diff(range(fin)), 0.1 * span, 1e-8)
    ylim <- range(fin) + c(-0.5, 0.5) * width
  } else {
    ylim <- range(c(lo, hi)) + c(-half, half)
  }
  p <- ggplot2::ggplot(d) +
    ggplot2::geom_rect(ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                    ymin = .data$ymin, ymax = .data$ymax,
                                    fill = .data$decision)) +
    ggplot2::scale_fill_manual(values = cols, drop = TRUE) +
    ggplot2::coord_cartesian(ylim = ylim, expand = FALSE) +
    ggplot2::labs(x = axes$xlab, y = axes$ylab, fill = "Decision",
                  caption = axes$caption) +
    style$theme +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  if (!is.null(lab)) {
    p <- p + ggplot2::geom_text(
      data = lab,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$text,
                   vjust = .data$vjust),
      colour = "white", fontface = "bold", size = style$label_size
    )
  }
  # Whole numbers only on the axes of counts and sample sizes
  int_breaks <- function(lim) {
    b <- pretty(lim)
    b[abs(b - round(b)) < 1e-8]
  }
  if (axes$discrete_x) {
    p <- p + ggplot2::scale_x_continuous(breaks = int_breaks)
  }
  sec <- if (is.null(axes$sec_scale)) {
    ggplot2::waiver()
  } else {
    s <- axes$sec_scale
    ggplot2::sec_axis(~ . / s, name = "Observed treatment response rate",
                      labels = function(v) paste0(round(100 * v), "%"))
  }
  # The secondary axis is used only with counts on the vertical axis (one
  # binary endpoint with a concurrent control group)
  if (axes$discrete_y) {
    p <- p + ggplot2::scale_y_continuous(breaks = int_breaks, sec.axis = sec)
  }
  p
}
