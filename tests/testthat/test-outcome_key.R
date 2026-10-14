test_that("outcome_key() distinguishes count vectors", {
  x <- outcomes_binary2(6)
  expect_identical(anyDuplicated(outcome_key(x)), 0L)
  expect_equal(outcome_key(matrix(c(1L, 2L, 3L, 4L), 1L)),
               2 + 5 * 3 + 5 ^ 2 * 4)
  # The base grows with the counts (a fixed base of 1001 would give the
  # same key to these two rows)
  y <- rbind(c(0L, 1001L, 0L, 0L), c(1000L, 0L, 1L, 0L))
  expect_false(outcome_key(y)[1] == outcome_key(y)[2])
  expect_length(outcome_key(matrix(integer(0), 0L, 4L)), 0L)
})
