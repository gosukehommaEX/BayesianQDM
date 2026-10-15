# Python_reference/ref_plot_overlay.py: the decision with the highest
# probability (ties within 1e-10 go to the first of Go, Gray, NoGo, Miss)
# and the step of its probability
test_that("dominant_decision() chooses the decision with the highest probability", {
  x <- data.frame(Go = c(0.62, 0.10, 0.05, 0.4, 0.2, 0.2, NA, 0.72, 0.45),
                  Gray = c(0.30, 0.45, 0.15, 0.4, 0.3, 0.3, NA, 0.2,
                           0.45 + 1e-12),
                  NoGo = c(0.08, 0.45, 0.80, 0.2, 0.5 - 1e-11, 0.499, NA,
                           0.08, 0.1))
  d <- dominant_decision(x, c("Go", "Gray", "NoGo"))
  expect_identical(names(d), c("decision", "prob", "step"))
  expect_identical(d$decision, c("Go", "Gray", "NoGo", "Go", "NoGo", "NoGo",
                                 NA, "Go", "Go"))
  expect_identical(d$prob, c(0.62, 0.45, 0.80, 0.4, 0.5 - 1e-11, 0.499, NA,
                             0.72, 0.45))
  expect_identical(d$step, c(3L, 1L, 5L, 1L, 2L, 1L, NA, 4L, 1L))
})

test_that("dominant_decision() includes Miss decisions and every step", {
  x <- data.frame(Go = 0.3, Gray = 0.1, NoGo = 0.05, Miss = 0.55)
  d <- dominant_decision(x, c("Go", "Gray", "NoGo", "Miss"))
  expect_identical(d$decision, "Miss")
  expect_identical(d$step, 2L)
  # Without Miss among the decisions, Go has the highest probability
  expect_identical(dominant_decision(x, c("Go", "Gray", "NoGo"))$decision,
                   "Go")
  p <- c(0, 0.4999, 0.5, 0.6, 0.7, 0.79, 0.8, 1)
  y <- data.frame(Go = p, Gray = 0, NoGo = 0)
  expect_identical(dominant_decision(y, c("Go", "Gray", "NoGo"))$step,
                   c(1L, 1L, 2L, 3L, 4L, 4L, 5L, 5L))
})
