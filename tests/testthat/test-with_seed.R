test_that("with_seed() gives reproducible results", {
  x1 <- with_seed(123, stats::runif(3))
  x2 <- with_seed(123, stats::runif(3))
  expect_identical(x1, x2)
})

test_that("with_seed() restores the random number state", {
  set.seed(42)
  state <- get(".Random.seed", envir = globalenv())
  with_seed(1, stats::runif(5))
  expect_identical(get(".Random.seed", envir = globalenv()), state)
})

test_that("with_seed() with NULL uses the current random number stream", {
  set.seed(7)
  x1 <- with_seed(NULL, stats::runif(2))
  set.seed(7)
  x2 <- stats::runif(2)
  expect_identical(x1, x2)
})

test_that("with_seed() removes .Random.seed if none existed before", {
  env <- globalenv()
  stats::runif(1)
  saved <- get(".Random.seed", envir = env)
  rm(".Random.seed", envir = env)
  with_seed(1, stats::runif(1))
  expect_false(exists(".Random.seed", envir = env, inherits = FALSE))
  assign(".Random.seed", saved, envir = env)
})

test_that("with_seed() validates the seed", {
  expect_error(with_seed(1.5, 1), "'seed' must contain whole numbers only.")
})
