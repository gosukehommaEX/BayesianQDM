#' Value of an argument as R code
#'
#' @param x A value of an argument of a \code{design_*()} function: a prior
#'   distribution, a decision criterion, data, external data, a hypothetical
#'   control, the conduct of a time-to-event study, a matrix, or a numeric,
#'   character or logical vector.
#' @param keep_type Whether integers keep the suffix \code{L} (see
#'   \code{code_num()}).
#'
#' @return A character string of R code or a call from \code{code_node()}.
#' @noRd
code_value <- function(x, keep_type = TRUE) {
  # Elements that the functions of the package convert to integers
  whole <- c("n", "y", "counts", "events")
  fields <- function(fun, nms) {
    code_node(fun, stats::setNames(lapply(nms, function(nm) {
      code_value(x[[nm]], keep_type = !(nm %in% whole))
    }), nms))
  }
  if (inherits(x, "qdm_prior")) {
    return(fields(sub("^qdm_", "", class(x)[1L]), names(unclass(x))))
  }
  if (inherits(x, "qdm_data_binary2")) {
    return(fields("data_binary2", "counts"))
  }
  if (inherits(x, "qdm_data")) {
    return(fields(sub("^qdm_", "", class(x)[1L]), names(unclass(x))))
  }
  if (inherits(x, "qdm_borrow")) {
    return(fields("borrow", c("data", "a0")))
  }
  if (inherits(x, "qdm_hypothetical")) {
    return(if (x$type == "data") {
      fields("hypothetical", "data")
    } else {
      fields("hypothetical", c("mean", "r"))
    })
  }
  if (inherits(x, "qdm_criterion")) {
    return(code_criterion(x))
  }
  if (inherits(x, "qdm_conduct")) {
    return(code_conduct(x))
  }
  if (is.matrix(x)) {
    return(code_node("matrix", list(code_num(as.vector(x), keep_type),
                                    nrow = code_num(nrow(x), FALSE))))
  }
  if (is.character(x)) {
    s <- vapply(x, encodeString, character(1), quote = "\"", USE.NAMES = FALSE)
  } else if (is.logical(x)) {
    s <- as.character(x)
  } else if (is.numeric(x)) {
    return(code_num(x, keep_type))
  } else {
    stop("R code cannot be written for this value.", call. = FALSE)
  }
  if (length(s) == 1L) s else paste0("c(", paste(s, collapse = ", "), ")")
}
