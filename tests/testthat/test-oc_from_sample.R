test_that("oc_from_sample() sums the probabilities of the outcomes", {
  s <- list(p_go = c(0.9, 0.6, 0.2, 0.05),
            p_nogo = c(0.01, 0.3, 0.7, 0.9),
            weights = cbind(c(0.1, 0.2, 0.3, 0.4), c(0.25, 0.25, 0.25, 0.25)))
  out <- oc_from_sample(s, gamma_go = 0.5, gamma_nogo = 0.8)
  expect_identical(colnames(out), c("Go", "NoGo", "Gray", "Miss"))
  expect_equal(out[1, ], c(Go = 0.3, NoGo = 0.4, Gray = 0.3, Miss = 0))
  expect_equal(out[2, ], c(Go = 0.5, NoGo = 0.25, Gray = 0.25, Miss = 0))
  out2 <- oc_from_sample(s, gamma_go = 0.5, gamma_nogo = 0.25)
  expect_equal(out2[1, ], c(Go = 0.1, NoGo = 0.7, Gray = 0, Miss = 0.2))
})

test_that("oc_from_sample() accepts probabilities that differ by scenario", {
  s <- list(p_go = cbind(c(0.9, 0.1), c(0.1, 0.9)),
            p_nogo = cbind(c(0.05, 0.9), c(0.9, 0.05)),
            weights = cbind(c(0.5, 0.5), c(0.2, 0.8)))
  out <- oc_from_sample(s, 0.5, 0.5)
  expect_equal(out[, "Go"], c(0.5, 0.8))
  expect_equal(out[, "NoGo"], c(0.5, 0.2))
})

test_that("oc_from_sample() handles pairs of outcomes and marginal rules", {
  # Pairs of 2 treatment and 1 control outcome, 2 scenarios; marginal rule
  s <- list(p_go = list(matrix(c(0.9, 0.2), 2L), matrix(c(0.8, 0.9), 2L)),
            p_nogo = list(matrix(c(0.1, 0.7), 2L), matrix(c(0.1, 0.1), 2L)),
            weights = NULL,
            weights_t = cbind(c(0.4, 0.6), c(NA, NA)),
            weights_c = matrix(1, 1L, 2L))
  out <- oc_from_sample(s, c(0.5, 0.5), c(0.5, 0.5), "or")
  expect_equal(out[1, ], c(Go = 0.4, NoGo = 0, Gray = 0, Miss = 0.6))
  expect_true(all(is.na(out[2, ])))
  and <- oc_from_sample(s, c(0.5, 0.5), c(0.5, 0.5), "and")
  expect_equal(and[1, ], c(Go = 0.4, NoGo = 0, Gray = 0.6, Miss = 0))
})
