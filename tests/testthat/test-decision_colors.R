test_that("decision_colors() names the colors by decision", {
  expect_identical(
    decision_colors(col_go = "#658D1B", col_nogo = "red", col_gray = "gray",
                    col_miss = "#E69F00"),
    c(Go = "#658D1B", NoGo = "red", Gray = "gray", Miss = "#E69F00")
  )
  expect_identical(decision_colors(col_nogo = "blue", col_go = "black"),
                   c(NoGo = "blue", Go = "black"))
})

test_that("decision_colors() names the invalid argument", {
  expect_error(decision_colors(col_go = "red", col_gray = "grey500"),
               "'col_gray' must be a single color name", fixed = TRUE)
})
