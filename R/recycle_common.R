#' Recycle arguments to a common length
#'
#' Recycles a named list of vectors to the length of the longest element.
#' Each element must have length 1 or the common length.
#'
#' @param args A named list of vectors.
#'
#' @return A named list of vectors of the common length.
#' @noRd
recycle_common <- function(args) {
  lens <- lengths(args)
  n <- max(lens)
  bad <- names(args)[lens != 1L & lens != n]
  if (length(bad) > 0L) {
    stop(sprintf("%s must have length 1 or %d.",
                 paste0("'", bad, "'", collapse = ", "), n),
         call. = FALSE)
  }
  lapply(args, rep_len, length.out = n)
}
