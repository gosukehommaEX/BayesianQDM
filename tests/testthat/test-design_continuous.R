test_that("design_continuous() stores the specification", {
  des <- ref_design_cont("C2")
  expect_s3_class(des, c("qdm_design_continuous", "qdm_design"),
                  exact = TRUE)
  expect_identical(des$endpoint, "continuous")
  expect_identical(des$n_t, 10L)
  expect_identical(des$n_c, 8L)
  expect_identical(des$control, "concurrent")
  expect_identical(des$method, "NI")
  expect_identical(des$nMC, 10000L)
  expect_s3_class(des$prior_c, "qdm_prior_nix")
  expect_equal(des$external_c$a0, 0.4)
  expect_null(des$hypothetical_c)
  expect_identical(ref_design_cont("C3")$control, "hypothetical")
  expect_null(ref_design_cont("C3")$prior_c)
  expect_identical(ref_design_cont("C4")$control, "external")
  mc <- ref_design_cont("C1", method = "MC", nMC = 500)
  expect_identical(mc$method, "MC")
  expect_identical(mc$nMC, 500L)
  expect_identical(ref_design_cont("C1", method = "MM")$method, "MM")
})

test_that("design_continuous() validates the sample sizes and priors", {
  vague <- prior_vague()
  expect_error(design_continuous(1, 12, vague, vague, crit_cont()),
               "'n_t' must be at least 2.")
  expect_error(design_continuous(12, 1, vague, vague, crit_cont()),
               "'n_c' must be 0 or a whole number of at least 2.")
  expect_error(design_continuous(12, 2.5, vague, vague, crit_cont()),
               "'n_c' must contain whole numbers only.")
  expect_error(design_continuous(12, 12, prior_beta(1, 1), vague,
                                 crit_cont()),
               "'prior_t' must be created by prior_nix\\(\\) or prior_vague")
  expect_error(design_continuous(12, 12, vague, criterion = crit_cont()),
               "'prior_c' must be created by prior_nix\\(\\) or prior_vague")
  expect_error(design_continuous(12, 0, vague, vague, crit_cont(),
                                 hypothetical_c = hypothetical(mean = 0,
                                                               r = 1)),
               "'prior_c' must be NULL for a hypothetical control")
})

test_that("design_continuous() validates the other arguments", {
  vague <- prior_vague()
  expect_error(design_continuous(12, 12, vague, vague,
                                 criterion_posterior(c(1, 1), c(0, 0), 0.5,
                                                     0.5)),
               "'criterion' must be created by criterion_posterior\\(\\) or")
  expect_error(design_continuous(12, 12, vague, vague, crit_cont(),
                                 external_t = borrow(data_binary(10, 3),
                                                     0.5)),
               "'external_t' must be NULL or created by borrow\\(\\) with")
  expect_error(design_continuous(12, 0, vague, criterion = crit_cont(),
                                 hypothetical_c = hypothetical(
                                   data_binary(10, 3))),
               "'hypothetical_c' must be NULL or created by hypothetical")
  expect_error(design_continuous(12, 0, vague, criterion = crit_cont(),
                                 hypothetical_c = hypothetical(
                                   mean = c(0, 1), r = 1)),
               "'hypothetical_c' must be NULL or created by hypothetical")
  expect_error(design_continuous(12, 12, vague, criterion = crit_cont(),
                                 hypothetical_c = hypothetical(mean = 0,
                                                               r = 1)),
               "'hypothetical_c' can be used only when 'n_c' is 0.")
  expect_error(design_continuous(12, 0, vague, vague, crit_cont()),
               "When 'n_c' is 0, specify 'external_c'")
  expect_error(design_continuous(12, 12, vague, vague, crit_cont(),
                                 method = "XX"),
               "'arg' should be one of")
  expect_error(design_continuous(12, 12, vague, vague, crit_cont(),
                                 method = "MC", nMC = 0),
               "'nMC' must be at least 1 and at most 2147483647.")
})

test_that("design_continuous() requires a proper control distribution", {
  # With the vague prior and external control data only, the control mean
  # has a0 * n_e - 1 degrees of freedom
  vague <- prior_vague()
  expect_error(design_continuous(12, 0, vague, vague, crit_cont(),
                                 external_c = borrow(data_continuous(4, 0,
                                                                     1),
                                                     0.25)),
               "'a0' times the number of external subjects must be greater")
  des <- design_continuous(12, 0, vague, vague, crit_cont(),
                           external_c = borrow(data_continuous(5, 0, 1),
                                               0.25))
  expect_identical(des$control, "external")
})

test_that("design_continuous() checks the degrees of freedom for MM", {
  vague <- prior_vague()
  expect_error(design_continuous(5, 12, vague, vague, crit_cont(),
                                 method = "MM"),
               "\\(here 4 and 11\\)")
  expect_error(design_continuous(12, 0, vague, vague, crit_cont(),
                                 external_c = borrow(data_continuous(10, 0,
                                                                     1),
                                                     0.5),
                                 method = "MM"),
               "\\(here 11 and 4\\)")
  # A NIX prior adds nu0 degrees of freedom
  des <- design_continuous(5, 5, prior_nix(0, 1, 1, 1), prior_nix(0, 1, 1, 1),
                           crit_cont(), method = "MM")
  expect_identical(des$method, "MM")
  # The other methods do not need the fourth moments
  expect_identical(design_continuous(2, 2, vague, vague, crit_cont())$n_t,
                   2L)
})
