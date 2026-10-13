test_that("cellprob_binary2() matches the closed-form probabilities", {
  # Reference values computed in Python from p11 = pi1 pi2 + rho sd1 sd2
  expect_equal(cellprob_binary2(0.3, 0.4, 0.3),
               c("00" = 0.48734983296193085, "01" = 0.2126501670380691,
                 "10" = 0.11265016703806907, "11" = 0.18734983296193092),
               tolerance = 1e-12)
  expect_equal(cellprob_binary2(0.3, 0.4, -0.2),
               c("00" = 0.37510011135871263, "01" = 0.3248998886412873,
                 "10" = 0.2248998886412873, "11" = 0.0751001113587127),
               tolerance = 1e-12)
})

test_that("cellprob_binary2() gives independence for rho = 0", {
  p <- cellprob_binary2(0.25, 0.35, 0)
  expect_equal(unname(p["11"]), 0.25 * 0.35, tolerance = 1e-14)
  expect_equal(unname(p["10"] + p["11"]), 0.25, tolerance = 1e-14)
  expect_equal(unname(p["01"] + p["11"]), 0.35, tolerance = 1e-14)
  expect_equal(sum(p), 1, tolerance = 1e-14)
})

test_that("cellprob_binary2() rejects infeasible correlations", {
  # For pi1 = 0.1 and pi2 = 0.9 the largest feasible correlation is 1/9
  expect_silent(cellprob_binary2(0.1, 0.9, 1 / 9))
  expect_error(cellprob_binary2(0.1, 0.9, 0.5), "is infeasible")
  expect_error(cellprob_binary2(0, 0.5, 0), "'pi1' must be greater than 0")
})
