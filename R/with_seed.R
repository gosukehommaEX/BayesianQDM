#' Evaluate an expression with a temporary random seed
#'
#' If \code{seed} is \code{NULL}, \code{expr} is evaluated with the current
#' random number stream. Otherwise the random number generator is seeded with
#' \code{seed}, \code{expr} is evaluated, and the previous state of the
#' generator (\code{.Random.seed} in the global environment) is restored on
#' exit.
#'
#' @param seed A single number or \code{NULL}.
#' @param expr Expression to evaluate. It is evaluated lazily, after seeding.
#'
#' @return The value of \code{expr}.
#' @noRd
with_seed <- function(seed, expr) {
  if (is.null(seed)) {
    return(expr)
  }
  check_numeric(seed, "seed", integer = TRUE)
  env <- globalenv()
  has_seed <- exists(".Random.seed", envir = env, inherits = FALSE)
  if (has_seed) {
    old_seed <- get(".Random.seed", envir = env, inherits = FALSE)
  }
  on.exit({
    if (has_seed) {
      assign(".Random.seed", old_seed, envir = env)
    } else if (exists(".Random.seed", envir = env, inherits = FALSE)) {
      rm(".Random.seed", envir = env)
    }
  }, add = TRUE)
  set.seed(seed)
  expr
}
