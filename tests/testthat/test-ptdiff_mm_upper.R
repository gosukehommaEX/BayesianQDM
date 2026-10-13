test_that("ptdiff_mm_upper() follows Proposition 1 of Yamaguchi et al.", {
  # With equal scales s and degrees of freedom nu, nu* = 2 nu - 4 and
  # scale* = s sqrt(2 nu (nu* - 2) / ((nu - 2) nu*)).
  for (nu in c(5, 10, 30)) {
    nu_star <- 2 * nu - 4
    scale_star <- sqrt(2 * nu * (nu_star - 2) / ((nu - 2) * nu_star))
    expect_equal(ptdiff_mm_upper(1.3, 0, 1, nu, 0, 1, nu),
                 stats::pt(1.3 / scale_star, nu_star, lower.tail = FALSE),
                 tolerance = 1e-12, info = paste("nu =", nu))
  }
})

test_that("ptdiff_mm_upper() reduces to the normal case for infinite df", {
  expect_equal(ptdiff_mm_upper(1, 2, 1, Inf, 0, 1.5, Inf),
               stats::pnorm(1, 2, sqrt(1 + 1.5 ^ 2), lower.tail = FALSE),
               tolerance = 1e-12)
})
