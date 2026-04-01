#' Plot Method for pbayesdecisionprob2bin Objects
#'
#' Displays operating characteristics for two-binary-endpoint results returned
#' by \code{\link{pbayesdecisionprob2bin}}.
#'
#' When the input scenarios form a regular grid over
#' \code{(pi_t1, pi_t2)} (i.e., every combination of the unique values of
#' \code{pi_t1} and \code{pi_t2} is present) and \code{rho_t} is constant,
#' the function produces a \strong{filled tile plot}: each panel (Go, Gray,
#' NoGo) is coloured by its own probability on a continuous gradient (white to
#' the panel colour), so intensity directly reflects the probability magnitude.
#' A solid white contour line is overlaid at the corresponding decision
#' threshold (\code{gamma_go} for the Go panel, \code{gamma_nogo} for the
#' NoGo panel, and their mean for the Gray panel) to mark the boundary where
#' the probability equals the threshold.  Otherwise the function falls back to
#' a \strong{scatter plot} in which point colour encodes the decision
#' probability on a continuous scale.
#'
#' When \code{which = "all"}, the three panels are arranged side-by-side using
#' \code{gridExtra::grid.arrange}, so each panel retains its own independent
#' colour scale.  This requires the \pkg{gridExtra} package.
#'
#' For \code{design = 'controlled'} or \code{design = 'external'}, both axes
#' are expressed as treatment-minus-control differences:
#' \eqn{\theta_1 = \pi_{t1} - \bar{\pi}_{c1}} and
#' \eqn{\theta_2 = \pi_{t2} - \bar{\pi}_{c2}},
#' where \eqn{\bar{\pi}_{c1}} and \eqn{\bar{\pi}_{c2}} are the means of the
#' supplied \code{pi_c1} and \code{pi_c2} vectors.
#' For \code{design = 'uncontrolled'}, the axes represent \eqn{\pi_{t1}} and
#' \eqn{\pi_{t2}} directly.
#'
#' Vertical and horizontal reference lines are drawn at the decision thresholds:
#' \itemize{
#'   \item When \code{prob = 'posterior'}: vertical lines at \eqn{\theta_{TV1}}
#'         and \eqn{\theta_{MAV1}} (x-axis) and horizontal lines at
#'         \eqn{\theta_{TV2}} and \eqn{\theta_{MAV2}} (y-axis).
#'   \item When \code{prob = 'predictive'}: a single vertical line at
#'         \eqn{\theta_{NULL1}} and a single horizontal line at
#'         \eqn{\theta_{NULL2}}.
#' }
#'
#' @param x An object of class \code{pbayesdecisionprob2bin}.
#' @param which A character string specifying which decision probability to
#'        plot.  Must be one of \code{"Go"}, \code{"Gray"}, \code{"NoGo"},
#'        \code{"all"}, or \code{"overlay"}.  Default is \code{"Go"}.
#' @param title A character string for the plot title.  Defaults to
#'        \code{NULL} (no title displayed).
#' @param xlab A character string or expression for the x-axis label.
#'        Defaults to \code{NULL}, which auto-generates a label based on
#'        \code{design}.
#' @param ylab A character string or expression for the y-axis label.
#'        Defaults to \code{NULL}, which auto-generates a label based on
#'        \code{design}.
#' @param col_go A character string specifying the high-end fill colour for the
#'        Go probability gradient.  Default is \code{"#658D1B"}.
#' @param col_nogo A character string specifying the high-end fill colour for
#'        the NoGo probability gradient.  Default is \code{"#D91E49"}.
#' @param col_gray A character string specifying the high-end fill colour for
#'        the Gray probability gradient.  Default is \code{"#939597"}.
#' @param base_size A positive numeric scalar specifying the base font size
#'        (in points) passed to \code{theme_bw()}.  Default is \code{28}.
#' @param ... Further arguments passed to or from other methods (ignored).
#'
#' @return Invisibly returns a \code{ggplot} object (single panel) or a
#'         \code{gtable} object (\code{which = "all"}).
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_vline geom_hline annotate
#'   scale_fill_gradient scale_fill_stepsn scale_fill_manual
#'   scale_x_continuous scale_y_continuous
#'   coord_cartesian labs theme_bw theme element_text element_line element_blank
#'   margin unit geom_point scale_color_gradient scale_color_stepsn geom_text
#'   guide_colorsteps
#' @importFrom grDevices col2rgb rgb
#' @importFrom gridExtra grid.arrange
#' @export
plot.pbayesdecisionprob2bin <- function(x,
                                        which     = "Go",
                                        title     = NULL,
                                        xlab      = NULL,
                                        ylab      = NULL,
                                        col_go    = "#658D1B",
                                        col_nogo  = "#D91E49",
                                        col_gray  = "#939597",
                                        base_size = 28,
                                        ...) {

  # --- Input validation ---
  if (!which %in% c("Go", "Gray", "NoGo", "all", "overlay")) {
    stop("'which' must be one of \"Go\", \"Gray\", \"NoGo\", \"all\", or \"overlay\"")
  }
  if (nrow(x) < 2L) {
    stop("'x' must contain at least 2 scenarios (rows) to produce a plot")
  }

  # --- Extract attributes ---
  design      <- attr(x, "design")
  prob        <- attr(x, "prob")
  theta_TV1   <- attr(x, "theta_TV1")
  theta_MAV1  <- attr(x, "theta_MAV1")
  theta_TV2   <- attr(x, "theta_TV2")
  theta_MAV2  <- attr(x, "theta_MAV2")
  theta_NULL1 <- attr(x, "theta_NULL1")
  theta_NULL2 <- attr(x, "theta_NULL2")
  gamma_go    <- attr(x, "gamma_go")
  gamma_nogo  <- attr(x, "gamma_nogo")

  # --- Build x/y axis variables ---
  if (design %in% c("controlled", "external")) {
    pi_c1_mean <- mean(x$pi_c1)
    pi_c2_mean <- mean(x$pi_c2)
    ax1     <- x$pi_t1 - pi_c1_mean
    ax2     <- x$pi_t2 - pi_c2_mean
    x_label <- if (is.null(xlab)) expression(theta[1] == pi[t1] - pi[c1]) else xlab
    y_label <- if (is.null(ylab)) expression(theta[2] == pi[t2] - pi[c2]) else ylab
  } else {
    pi_c1_mean <- NULL
    pi_c2_mean <- NULL
    ax1     <- x$pi_t1
    ax2     <- x$pi_t2
    x_label <- if (is.null(xlab)) expression(pi[t1]) else xlab
    y_label <- if (is.null(ylab)) expression(pi[t2]) else ylab
  }

  # --- Determine threshold positions on each axis ---
  # Threshold values are used directly as axis coordinates (no transformation).
  if (prob == "posterior") {
    vx_TV  <- theta_TV1;  vx_MAV <- theta_MAV1
    vy_TV  <- theta_TV2;  vy_MAV <- theta_MAV2
  } else {
    vx_NULL <- theta_NULL1
    vy_NULL <- theta_NULL2
  }

  # --- Detect grid layout and constant rho_t ---
  u1      <- sort(unique(ax1))
  u2      <- sort(unique(ax2))
  is_grid <- (length(u1) * length(u2) == nrow(x)) &&
    (length(unique(x$rho_t)) == 1L)

  # --- Helper: axis breaks ---
  axis_breaks <- function(vals) {
    b <- pretty(vals, n = 6L)
    eps <- 1e-9
    b[b >= (min(vals) - eps) & b <= (max(vals) + eps)]
  }

  # --- Helper: add threshold reference lines and labels ---
  add_thresholds <- function(p, x_rng, y_rng, bs) {
    off_x <- diff(x_rng) * 0.02
    off_y <- diff(y_rng) * 0.02
    if (prob == "posterior") {
      p <- p +
        ggplot2::geom_vline(xintercept = vx_MAV, color = "gray30",
                            linetype = "dotted", linewidth = 1.0) +
        ggplot2::geom_vline(xintercept = vx_TV,  color = "gray30",
                            linetype = "dotted", linewidth = 1.0) +
        ggplot2::geom_hline(yintercept = vy_MAV, color = "gray30",
                            linetype = "dotted", linewidth = 1.0) +
        ggplot2::geom_hline(yintercept = vy_TV,  color = "gray30",
                            linetype = "dotted", linewidth = 1.0) +
        ggplot2::annotate("text",
                          x = vx_MAV - off_x, y = y_rng[2],
                          label = "theta[MAV1]", parse = TRUE,
                          color = "gray30", hjust = 1, vjust = 1,
                          size = bs / 3) +
        ggplot2::annotate("text",
                          x = vx_TV  + off_x, y = y_rng[2],
                          label = "theta[TV1]",  parse = TRUE,
                          color = "gray30", hjust = 0, vjust = 1,
                          size = bs / 3) +
        ggplot2::annotate("text",
                          x = x_rng[1], y = vy_MAV - off_y,
                          label = "theta[MAV2]", parse = TRUE,
                          color = "gray30", hjust = 0, vjust = 1,
                          size = bs / 3) +
        ggplot2::annotate("text",
                          x = x_rng[1], y = vy_TV  + off_y,
                          label = "theta[TV2]",  parse = TRUE,
                          color = "gray30", hjust = 0, vjust = 0,
                          size = bs / 3)
    } else {
      p <- p +
        ggplot2::geom_vline(xintercept = vx_NULL, color = "gray30",
                            linetype = "dotted", linewidth = 1.0) +
        ggplot2::geom_hline(yintercept = vy_NULL, color = "gray30",
                            linetype = "dotted", linewidth = 1.0) +
        ggplot2::annotate("text",
                          x = vx_NULL + off_x, y = y_rng[2],
                          label = "theta[NULL1]", parse = TRUE,
                          color = "gray30", hjust = 0, vjust = 1,
                          size = bs / 3) +
        ggplot2::annotate("text",
                          x = x_rng[1], y = vy_NULL + off_y,
                          label = "theta[NULL2]", parse = TRUE,
                          color = "gray30", hjust = 0, vjust = 0,
                          size = bs / 3)
    }
    p
  }

  # --- Helper: common theme ---
  common_theme <- function(bs) {
    ggplot2::theme_bw(base_size = bs) +
      ggplot2::theme(
        legend.position   = "right",
        legend.margin     = ggplot2::margin(t = 0, r = 0, b = 0, l = 10),
        legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 5),
        legend.text       = ggplot2::element_text(size = bs * 0.74),
        legend.title      = ggplot2::element_text(size = bs * 0.74,
                                                  hjust = 0),
        legend.key.height = ggplot2::unit(bs * 0.04, "cm"),
        panel.grid.minor  = ggplot2::element_blank(),
        panel.grid.major  = ggplot2::element_blank(),
        panel.border      = ggplot2::element_blank(),
        axis.line         = ggplot2::element_line(color = "black",
                                                  linewidth = 0.8),
        strip.text        = ggplot2::element_text(size = bs * 0.64,
                                                  face = "bold")
      )
  }

  # ---------------------------------------------------------------------------
  # OVERLAY MODE
  # ---------------------------------------------------------------------------
  if (which == "overlay") {

    if (!is_grid) {
      stop("'which = \"overlay\"' requires a regular grid input.")
    }
    if (!requireNamespace("gridExtra", quietly = TRUE)) {
      stop("Package 'gridExtra' is required for which = \"overlay\".")
    }

    tile_w_ov <- if (length(unique(ax1)) > 1L) min(diff(sort(unique(ax1)))) / 2 else 0.05
    tile_h_ov <- if (length(unique(ax2)) > 1L) min(diff(sort(unique(ax2)))) / 2 else 0.05
    x_rng_ov  <- c(min(ax1) - tile_w_ov, max(ax1) + tile_w_ov)
    y_rng_ov  <- c(min(ax2) - tile_h_ov, max(ax2) + tile_h_ov)

    df_ov <- data.frame(ax1 = ax1, ax2 = ax2,
                        Go = x[["Go"]], Gray = x[["Gray"]], NoGo = x[["NoGo"]])
    df_ov$dominant <- apply(df_ov[, c("Go","Gray","NoGo")], 1L, function(r) {
      if (any(is.na(r))) return(NA_character_)
      c("Go","Gray","NoGo")[which.max(r)]
    })
    df_ov$max_prob <- apply(df_ov[, c("Go","Gray","NoGo")], 1L, function(r) {
      if (any(is.na(r))) return(NA_real_)
      max(r, na.rm = TRUE)
    })
    bin_labels <- c("< 0.50","0.50-0.60","0.60-0.70","0.70-0.80","> 0.80")
    df_ov$prob_bin <- cut(df_ov$max_prob,
                          breaks = c(0, 0.50, 0.60, 0.70, 0.80, 1.001),
                          labels = bin_labels, include.lowest = TRUE, right = FALSE)

    make_shades <- function(base_col, n = 5L) {
      alphas <- seq(0.20, 1.00, length.out = n)
      sapply(alphas, function(a) {
        v <- col2rgb(base_col) / 255
        rgb(1 - a*(1-v[1]), 1 - a*(1-v[2]), 1 - a*(1-v[3]))
      })
    }
    go_sh   <- make_shades(col_go);   names(go_sh)   <- bin_labels
    gray_sh <- make_shades(col_gray); names(gray_sh) <- bin_labels
    nogo_sh <- make_shades(col_nogo); names(nogo_sh) <- bin_labels

    df_ov$fill_col <- mapply(function(dom, bin) {
      if (is.na(dom) || is.na(bin)) return("white")
      b <- as.character(bin)
      if (dom == "Go")   return(unname(go_sh[b]))
      if (dom == "Gray") return(unname(gray_sh[b]))
      if (dom == "NoGo") return(unname(nogo_sh[b]))
      "white"
    }, df_ov$dominant, df_ov$prob_bin)

    p_main <- ggplot2::ggplot(df_ov, ggplot2::aes(x = ax1, y = ax2)) +
      ggplot2::geom_tile(ggplot2::aes(fill = I(fill_col)), color = "gray50",
                         linewidth = 0.3) +
      ggplot2::scale_x_continuous(breaks = axis_breaks(ax1), expand = c(0, 0)) +
      ggplot2::scale_y_continuous(breaks = axis_breaks(ax2), expand = c(0, 0)) +
      ggplot2::coord_cartesian(xlim = x_rng_ov, ylim = y_rng_ov) +
      ggplot2::labs(title = title, x = x_label, y = y_label) +
      common_theme(base_size) +
      ggplot2::theme(legend.position = "none",
                     plot.margin = ggplot2::margin(t = 10, r = 5, b = 5, l = 5))

    df_na <- df_ov[is.na(df_ov$dominant), , drop = FALSE]
    if (nrow(df_na) > 0L) {
      p_main <- p_main +
        ggplot2::geom_text(data = df_na,
                           mapping = ggplot2::aes(x = ax1, y = ax2, label = "NA"),
                           color = "gray50", size = base_size * 0.35, inherit.aes = FALSE)
    }
    p_main <- add_thresholds(p_main, x_rng_ov, y_rng_ov, base_size)

    gray_5 <- make_shades("#808080")
    leg_df <- data.frame(
      y   = c(7, 6, 5, 3, 2, 1, 0, -1),
      col = c(unname(go_sh[5L]), unname(gray_sh[5L]), unname(nogo_sh[5L]),
              unname(gray_5)),
      lab = c("Go", "Grey", "NoGo", bin_labels),
      stringsAsFactors = FALSE
    )
    p_leg <- ggplot2::ggplot(leg_df, ggplot2::aes(x = 0, y = y)) +
      ggplot2::geom_tile(ggplot2::aes(fill = I(col)), width = 0.8, height = 0.8) +
      ggplot2::geom_text(ggplot2::aes(x = 0.6, label = lab),
                         hjust = 0, size = base_size * 0.26) +
      ggplot2::annotate("text", x = -0.4, y = 8.0, label = "Decision",
                        hjust = 0, fontface = "bold", size = base_size * 0.35) +
      ggplot2::annotate("text", x = -0.4, y = 4.0, label = "Probability",
                        hjust = 0, fontface = "bold", size = base_size * 0.35) +
      ggplot2::scale_x_continuous(limits = c(-0.5, 4.0), expand = c(0, 0)) +
      ggplot2::scale_y_continuous(limits = c(-1.6, 8.5), expand = c(0, 0)) +
      ggplot2::coord_fixed(ratio = 1) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.margin = ggplot2::margin(t = 2, r = 2, b = 2, l = 2))

    out <- gridExtra::grid.arrange(
      p_main, p_leg, ncol = 2L, widths = c(6, 1)
    )
    return(invisible(out))
  }

  # ---------------------------------------------------------------------------
  # CONTOUR MODE
  # ---------------------------------------------------------------------------
  if (is_grid) {

    tile_w <- if (length(u1) > 1L) min(diff(u1)) / 2 else 0.05
    tile_h <- if (length(u2) > 1L) min(diff(u2)) / 2 else 0.05
    x_rng  <- c(min(ax1) - tile_w, max(ax1) + tile_w)
    y_rng  <- c(min(ax2) - tile_h, max(ax2) + tile_h)

    # Helper: build a single gradient tile panel
    # - prob_col  : column name in x ("Go", "Gray", "NoGo")
    # - high_col  : colour for high-probability end of the gradient
    # - gamma_thr : threshold value at which the white contour line is drawn
    make_tile_panel <- function(prob_col, high_col, gamma_thr,
                                panel_title, bs) {
      df <- data.frame(ax1      = ax1,
                       ax2      = ax2,
                       prob_val = x[[prob_col]])

      legend_label <- switch(prob_col,
                             Go   = "Pr(Go)",
                             NoGo = "Pr(NoGo)",
                             Gray = "Pr(Gray)")

      p <- ggplot2::ggplot(df, ggplot2::aes(x = ax1, y = ax2)) +
        ggplot2::geom_tile(ggplot2::aes(fill = prob_val),
                           color = "gray50", linewidth = 0.5) +
        ggplot2::scale_fill_stepsn(
          name   = legend_label,
          colors = c("white", high_col),
          breaks = seq(0.1, 0.9, by = 0.1),
          limits = c(0, 1),
          guide  = ggplot2::guide_colorsteps(
            title.position = "top",
            title.hjust    = 0,
            barheight      = ggplot2::unit(bs * 0.50, "cm"),
            barwidth       = ggplot2::unit(bs * 0.06, "cm"),
            direction      = "vertical",
            show.limits    = TRUE
          )
        ) +
        ggplot2::scale_x_continuous(breaks = axis_breaks(ax1),
                                    expand = c(0, 0)) +
        ggplot2::scale_y_continuous(breaks = axis_breaks(ax2),
                                    expand = c(0, 0)) +
        ggplot2::coord_cartesian(xlim = x_rng, ylim = y_rng) +
        ggplot2::labs(title = panel_title, x = x_label, y = y_label) +
        common_theme(bs)

      # Overlay white tiles and "NA" text on missing cells
      df_na_tile <- df[is.na(df$prob_val), , drop = FALSE]
      if (nrow(df_na_tile) > 0L) {
        p <- p +
          ggplot2::geom_tile(data = df_na_tile,
                             mapping = ggplot2::aes(x = ax1, y = ax2),
                             fill = "white", color = "gray50",
                             linewidth = 0.5, inherit.aes = FALSE) +
          ggplot2::geom_text(data = df_na_tile,
                             mapping = ggplot2::aes(x = ax1, y = ax2, label = "NA"),
                             color = "gray50", size = bs * 0.35,
                             inherit.aes = FALSE)
      }

      p <- add_thresholds(p, x_rng, y_rng, bs)
      p
    }

    if (which == "all") {
      # Build three independent panels with panel-specific gradients and
      # arrange them side-by-side; gridExtra is required
      if (!requireNamespace("gridExtra", quietly = TRUE)) {
        stop("Package 'gridExtra' is required for which = \"all\". ",
             "Please install it with install.packages(\"gridExtra\").")
      }
      p_go   <- make_tile_panel("Go",   col_go,   gamma_go,   "Go",   base_size)
      p_gray <- make_tile_panel("Gray", col_gray,
                                (gamma_go + gamma_nogo) / 2, "Gray", base_size)
      p_nogo <- make_tile_panel("NoGo", col_nogo, gamma_nogo, "NoGo", base_size)

      out <- gridExtra::grid.arrange(
        p_go, p_gray, p_nogo,
        nrow = 1L,
        top  = if (!is.null(title)) title else ""
      )
      return(invisible(out))

    } else {
      high_col  <- switch(which, Go = col_go, Gray = col_gray, NoGo = col_nogo)
      gamma_thr <- switch(which,
                          Go   = gamma_go,
                          NoGo = gamma_nogo,
                          Gray = (gamma_go + gamma_nogo) / 2)
      p <- make_tile_panel(which, high_col, gamma_thr, title, base_size)
    }

    # ---------------------------------------------------------------------------
    # SCATTER MODE (non-grid or multiple rho_t values)
    # ---------------------------------------------------------------------------
  } else {

    prob_col <- if (which == "all") {
      message("Non-grid input: 'which = \"all\"' shows Go probability only.")
      "Go"
    } else {
      which
    }

    high_col  <- switch(prob_col,
                        Go   = col_go,
                        NoGo = col_nogo,
                        Gray = col_gray)
    prob_vals <- x[[prob_col]]
    x_rng     <- range(ax1)
    y_rng     <- range(ax2)

    df_scatter <- data.frame(ax1 = ax1, ax2 = ax2, prob_val = prob_vals)

    legend_label <- switch(prob_col,
                           Go   = "Pr(Go)",
                           NoGo = "Pr(NoGo)",
                           Gray = "Pr(Gray)")

    p <- ggplot2::ggplot(df_scatter,
                         ggplot2::aes(x = ax1, y = ax2, color = prob_val)) +
      ggplot2::geom_point(size = 4) +
      ggplot2::scale_color_stepsn(
        name   = legend_label,
        colors = c("white", high_col),
        breaks = seq(0.1, 0.9, by = 0.1),
        limits = c(0, 1),
        guide  = ggplot2::guide_colorsteps(
          title.position  = "top",
          title.hjust     = 0,
          barheight       = ggplot2::unit(base_size * 0.50, "cm"),
          barwidth        = ggplot2::unit(base_size * 0.06, "cm"),
          direction       = "vertical",
          ticks           = TRUE,
          ticks.colour    = "white",
          ticks.linewidth = 1.5
        )
      ) +
      ggplot2::scale_x_continuous(breaks = axis_breaks(ax1),
                                  expand = c(0.05, 0)) +
      ggplot2::scale_y_continuous(breaks = axis_breaks(ax2),
                                  expand = c(0.05, 0)) +
      ggplot2::labs(title = title, x = x_label, y = y_label) +
      common_theme(base_size)

    df_na_sc <- df_scatter[is.na(df_scatter$prob_val), , drop = FALSE]
    if (nrow(df_na_sc) > 0L) {
      p <- p +
        ggplot2::geom_text(data = df_na_sc,
                           mapping = ggplot2::aes(x = ax1, y = ax2, label = "NA"),
                           color = "gray50", size = base_size * 0.35,
                           inherit.aes = FALSE)
    }

    p <- add_thresholds(p, x_rng, y_rng, base_size)
  }

  print(p)
  invisible(p)
}
