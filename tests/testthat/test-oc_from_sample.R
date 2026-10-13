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
