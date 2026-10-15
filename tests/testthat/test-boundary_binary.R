test_that("boundary_binary() matches the fixed-control references", {
  # Reference: Python_reference/ref_fixed_tte_qdm.py, designs D8 and D9
  b <- boundary_binary(ref_design("D8"), 20)
  expect_identical(b$grid$y_t, 0:20)
  expect_identical(b$grid$n_t, rep(20L, 21L))
  expect_equal(b$p_go[12L], 0.915075681321, tolerance = 1e-10)
  expect_equal(b$p_nogo[12L], 0.00168707910052, tolerance = 1e-10)
  expect_equal(b$p_go, stats::pbeta(0.4, 1 + 0:20, 21 - 0:20,
                                    lower.tail = FALSE))
  ref_d9 <- c(0.0350877192982, 0.125902992776, 0.267801857585,
              0.438080495356, 0.608359133127, 0.755060728745,
              0.865086925458, 0.93653250774, 0.975827577995,
              0.993292053664, 0.998968008256)
  b9 <- boundary_binary(ref_design("D9"), 10)
  expect_equal(b9$p_go, ref_d9, tolerance = 1e-10)
  expect_equal(b9$p_go + b9$p_nogo, rep(1, 11L))
})

test_that("boundary_binary() stacks several sample sizes without a control", {
  b <- boundary_binary(ref_design("D8"), c(12, 10, 10))
  expect_identical(b$grid$n_t, rep(c(10L, 12L), c(11L, 13L)))
  expect_identical(b$grid$y_t, c(0:10, 0:12))
  expect_identical(b$axes$x, "n_t")
  expect_identical(b$axes$y, "y_t")
  expect_null(b$axes$caption)
  expect_equal(b$p_go[11L + 1:13],
               stats::pbeta(0.4, 1 + 0:12, 13 - 0:12, lower.tail = FALSE))
})

test_that("boundary_binary() covers all outcomes with a concurrent control", {
  des <- ref_design("D1")
  b <- boundary_binary(des, 12)
  expect_identical(nrow(b$grid), 13L * 13L)
  expect_identical(names(b$grid), c("n_t", "y_t", "n_c", "y_c"))
  expect_identical(b$axes$x, "y_c")
  expect_identical(b$axes$y, "y_t")
  expect_identical(b$axes$xlab, "Responders in the control group (n = 12)")
  expect_identical(b$axes$ylab,
                   "Responders in the treatment group (n = 12)")
  row <- which(b$grid$y_t == 9 & b$grid$y_c == 4)
  dec <- decide_binary(des, data_binary(12, 9), data_binary(12, 4))
  expect_equal(b$p_go[row], dec$prob_go)
  expect_equal(b$p_nogo[row], dec$prob_nogo)
  expect_equal(b$p_go[row], pbetadiff(0.2, 10, 4, 5, 9, lower.tail = FALSE))
  expect_error(boundary_binary(des, c(10, 12)),
               "'n_t' must be a single number for a design with a concurrent")
})

test_that("boundary_binary() uses hypothetical and external controls", {
  for (name in c("D3", "D4")) {
    des <- ref_design(name)
    b <- boundary_binary(des, 12)
    dec <- decide_binary(des, data_binary(12, 7), NULL)
    expect_equal(b$p_go[8L], dec$prob_go, info = name)
    expect_equal(b$p_nogo[8L], dec$prob_nogo, info = name)
  }
  expect_error(boundary_binary(ref_design("D3"), 0),
               "'n_t' must be at least 1.")
})
