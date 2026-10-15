#' Shades of a color by probability
#'
#' Interpolates linearly in RGB between white (probability 0) and a color
#' (probability 1), channel by channel, and rounds to an integer. Missing
#' probabilities give white.
#'
#' @param prob Numeric vector of probabilities.
#' @param color Colors, of length 1 or the length of \code{prob}.
#'
#' @return A character vector of hexadecimal color codes.
#' @noRd
shade_color <- function(prob, color) {
  rgb <- t(grDevices::col2rgb(rep_len(color, length(prob))))
  p <- pmin(pmax(prob, 0), 1)
  v <- round(255 + p * (rgb - 255))
  out <- sprintf("#%02X%02X%02X", as.integer(v[, 1L]), as.integer(v[, 2L]),
                 as.integer(v[, 3L]))
  out[is.na(prob)] <- "#FFFFFF"
  out
}
