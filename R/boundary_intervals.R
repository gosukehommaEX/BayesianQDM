#' Decision intervals of a decision boundary
#'
#' Computes the decision at each value of the grid of the observed quantity
#' (the cut variable) for each condition, and merges consecutive values with
#' the same decision into intervals. For a discrete cut variable (numbers of
#' responders or events), the intervals are ranges of grid values. For a
#' continuous cut variable (observed means), the points at which the Go or
#' NoGo criterion changes between two consecutive grid values are located by
#' bisection, the intervals are separated by these points, and the first and
#' last intervals extend to minus and plus infinity (decision changes
#' outside the grid are not detected, and neither are two changes of the
#' same criterion between two consecutive grid values).
#'
#' @param spec List with elements \code{values} (list of sorted numeric
#'   vectors, the grid of each condition), \code{prob} (function of condition
#'   indices and values that returns a list with the probabilities
#'   \code{go} and \code{nogo}, as numeric vectors or, for two endpoints with
#'   a marginal rule, lists of two numeric vectors) and \code{discrete}
#'   (logical).
#' @param criterion The criterion of the design.
#'
#' @return A data frame with the columns \code{row} (index of the
#'   condition), \code{decision} (character), \code{from} and \code{to}.
#' @noRd
boundary_intervals <- function(spec, criterion) {
  met <- function(idx, x) {
    p <- spec$prob(idx, x)
    list(go = meets_criterion(p$go, criterion$gamma_go, criterion$rule),
         nogo = meets_criterion(p$nogo, criterion$gamma_nogo,
                                criterion$rule))
  }
  idx <- rep(seq_along(spec$values), lengths(spec$values))
  x <- unlist(spec$values, use.names = FALSE)
  m <- met(idx, x)
  dec <- decision_label(m$go, m$nogo)
  n <- length(x)
  new_run <- c(TRUE, idx[-1L] != idx[-n] | dec[-1L] != dec[-n])
  if (spec$discrete) {
    start <- which(new_run)
    end <- c(start[-1L] - 1L, n)
    return(data.frame(row = idx[start], decision = dec[start],
                      from = x[start], to = x[end]))
  }
  # Cells (j, j + 1) of the same condition in which a criterion changes
  same <- idx[-1L] == idx[-n]
  cell_go <- which(same & m$go[-1L] != m$go[-n])
  cell_nogo <- which(same & m$nogo[-1L] != m$nogo[-n])
  pred_go <- function(i, v) met(i, v)$go
  pred_nogo <- function(i, v) met(i, v)$nogo
  flips <- rbind(
    data.frame(row = idx[cell_go], at = bisect_flip(pred_go, idx[cell_go],
                                                    x[cell_go],
                                                    x[cell_go + 1L]),
               which = rep("go", length(cell_go))),
    data.frame(row = idx[cell_nogo],
               at = bisect_flip(pred_nogo, idx[cell_nogo], x[cell_nogo],
                                x[cell_nogo + 1L]),
               which = rep("nogo", length(cell_nogo)))
  )
  first <- which(!duplicated(idx))
  out <- lapply(seq_along(spec$values), function(r) {
    f <- flips[flips$row == r, , drop = FALSE]
    f <- f[order(f$at), , drop = FALSE]
    go <- m$go[first[r]]
    nogo <- m$nogo[first[r]]
    decs <- decision_label(go, nogo)
    for (k in seq_len(nrow(f))) {
      if (f$which[k] == "go") go <- !go else nogo <- !nogo
      decs <- c(decs, decision_label(go, nogo))
    }
    from <- c(-Inf, f$at)
    to <- c(f$at, Inf)
    # Drop intervals of zero width (two changes at the same point) and merge
    # neighbors with the same decision
    keep <- to - from > 1e-12 * pmax(1, abs(from), abs(to), na.rm = TRUE) |
      is.infinite(from) | is.infinite(to)
    from <- from[keep]
    to <- to[keep]
    decs <- decs[keep]
    run <- c(TRUE, decs[-1L] != decs[-length(decs)])
    start <- which(run)
    end <- c(start[-1L] - 1L, length(decs))
    data.frame(row = rep(r, length(start)), decision = decs[start],
               from = from[start], to = to[end])
  })
  do.call(rbind, out)
}
