test_that("reflect_criterion() leaves a criterion with larger effects favorable unchanged", {
  cr <- crit_post()
  r <- reflect_criterion(cr)
  expect_identical(r$criterion, cr)
  expect_false(r$flip)
})

test_that("reflect_criterion() negates the thresholds of the reflected endpoints", {
  cr <- criterion_posterior(-0.1, 0.05, 0.6, 0.7, direction = "less")
  r <- reflect_criterion(cr)
  expect_true(r$flip)
  expect_identical(r$criterion$theta_tv, 0.1)
  expect_identical(r$criterion$theta_mav, -0.05)
  expect_identical(r$criterion$direction, "greater")
  expect_identical(r$criterion$gamma_go, 0.6)
  cp <- criterion_predictive(c(0, 0.1), m_t = 20, gamma_go = 0.5,
                             gamma_nogo = 0.5,
                             direction = c("greater", "less"))
  rp <- reflect_criterion(cp)
  expect_identical(rp$flip, c(FALSE, TRUE))
  expect_identical(rp$criterion$theta_null, c(0, -0.1))
  expect_identical(rp$criterion$direction, c("greater", "greater"))
})
