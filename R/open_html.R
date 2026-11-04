#' Open an HTML file
#'
#' Opens the file in the viewer of RStudio if it is available and the file is
#' in the temporary folder of the R session (the viewer shows only those
#' files), and in the web browser otherwise. The tests replace this function,
#' because it opens a window.
#'
#' @param file Name of the HTML file.
#'
#' @return \code{file}, invisibly.
#' @noRd
open_html <- function(file) {
  # nocov start
  viewer <- getOption("viewer")
  path <- normalizePath(file, winslash = "/", mustWork = FALSE)
  temp <- normalizePath(tempdir(), winslash = "/", mustWork = FALSE)
  if (is.function(viewer) && startsWith(path, temp)) {
    viewer(file)
  } else {
    utils::browseURL(file)
  }
  invisible(file)
  # nocov end
}
