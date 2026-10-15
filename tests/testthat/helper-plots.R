# Parts of a ggplot object. With ggplot2 4.0 and later the plot, its theme
# elements and their settings are S7 objects whose properties are stored as
# attributes; with earlier versions they are lists.
plot_theme <- function(p) {
  th <- attr(p, "theme")
  if (is.null(th)) p$theme else th
}

plot_data <- function(p) {
  d <- attr(p, "data")
  if (is.null(d)) p$data else d
}

element_size <- function(el) {
  size <- attr(el, "size")
  if (is.null(size)) el$size else size
}
