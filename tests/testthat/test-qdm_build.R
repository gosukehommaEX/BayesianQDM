# The references are the designs of the helper files, created by calling the
# design_*() functions directly.
test_that("qdm_build() gives the designs of design_binary()", {
  unif <- prior_beta(1, 1)
  spec <- qdm_endpoint("binary") |>
    qdm_prior(unif, unif) |>
    qdm_criterion(crit_post())
  expect_identical(spec |> qdm_sample_size(12, 12) |> qdm_build(),
                   ref_design("D1"))
  expect_identical(
    spec |>
      qdm_sample_size(10, 8) |>
      qdm_prior(unif, prior_beta(2, 3)) |>
      qdm_external(external_t = borrow(data_binary(10, 6), 0.5),
                   external_c = borrow(data_binary(20, 8), 0.5)) |>
      qdm_build(),
    ref_design("D2")
  )
  expect_identical(
    spec |>
      qdm_sample_size(12, 0) |>
      qdm_control(hypothetical_c = hypothetical(data_binary(12, 4))) |>
      qdm_build(),
    ref_design("D3")
  )
  expect_identical(
    spec |>
      qdm_sample_size(12, 0) |>
      qdm_external(external_c = borrow(data_binary(20, 6), 0.5)) |>
      qdm_build(),
    ref_design("D4")
  )
  expect_identical(
    spec |>
      qdm_sample_size(10, 10) |>
      qdm_criterion(criterion_predictive(0.05, m_t = 30, m_c = 20,
                                         gamma_go = 0.9, gamma_nogo = 0.7)) |>
      qdm_build(),
    ref_design("D5")
  )
  expect_identical(
    qdm_endpoint("binary") |>
      qdm_sample_size(20, 0) |>
      qdm_prior(unif) |>
      qdm_criterion(crit_post(0.8, 0.8, tv = 0.2, mav = 0.05)) |>
      qdm_control(fixed_c = 0.2) |>
      qdm_build(),
    ref_design("D8")
  )
})

test_that("qdm_build() gives the designs of design_continuous()", {
  vague <- prior_vague()
  expect_identical(
    qdm_endpoint("continuous") |>
      qdm_sample_size(10, 8) |>
      qdm_prior(prior_nix(5, 2, 4, 20), prior_nix(-5, 1, 3, 25)) |>
      qdm_criterion(crit_cont()) |>
      qdm_external(external_t = borrow(data_continuous(6, 20, 30), 0.5),
                   external_c = borrow(data_continuous(30, -12, 20), 0.4)) |>
      qdm_computation(method = "MM") |>
      qdm_build(),
    ref_design_cont("C2", method = "MM")
  )
  expect_identical(
    qdm_endpoint("continuous") |>
      qdm_sample_size(12, 0) |>
      qdm_prior(vague) |>
      qdm_criterion(crit_cont()) |>
      qdm_control(hypothetical_c = hypothetical(mean = -10, r = 1.5)) |>
      qdm_build(),
    ref_design_cont("C3")
  )
  expect_identical(
    qdm_endpoint("continuous") |>
      qdm_sample_size(12, 0) |>
      qdm_prior(vague) |>
      qdm_criterion(crit_cont()) |>
      qdm_control(fixed_c = -10) |>
      qdm_computation(method = "MC", nMC = 2000L) |>
      qdm_build(),
    ref_design_cont("C7", method = "MC", nMC = 2000L)
  )
})

test_that("qdm_build() gives the designs of design_binary2()", {
  expect_identical(
    qdm_endpoint("binary2") |>
      qdm_sample_size(8, 8) |>
      qdm_prior(dir_q(), dir_q()) |>
      qdm_criterion(criterion_posterior(c(0.15, 0.25), c(0.10, 0.15),
                                        gamma_go = 0.453, gamma_nogo = 0.229,
                                        go_regions = 1:2)) |>
      qdm_external(external_t = borrow(data_binary2(c(2, 1, 1, 2)), 0.7),
                   external_c = borrow(data_binary2(c(14, 2, 3, 1)), 0.5)) |>
      qdm_computation(nMC = 20000L) |>
      qdm_build(),
    ref_design_bin2("B6")
  )
  expect_identical(
    qdm_endpoint("binary2") |>
      qdm_sample_size(3, 0) |>
      qdm_prior(dir_q(), dir_q()) |>
      qdm_criterion(criterion_predictive(c(0.1, 0), m_t = 4, m_c = 3,
                                         gamma_go = 0.4, gamma_nogo = 0.3)) |>
      qdm_control(
        hypothetical_c = hypothetical(data_binary2(c(2, 1, 1, 0)))
      ) |>
      qdm_computation(nMC = 500L, enumerate = FALSE) |>
      qdm_build(),
    ref_design_bin2("B7", nMC = 500L, enumerate = FALSE)
  )
  expect_identical(
    qdm_endpoint("binary2") |>
      qdm_sample_size(3, 0) |>
      qdm_prior(dir_q()) |>
      qdm_criterion(criterion_predictive(c(0.1, 0), m_t = 4, gamma_go = 0.4,
                                         gamma_nogo = 0.3)) |>
      qdm_control(fixed_c = c(0.25, 0.5)) |>
      qdm_computation(nMC = 20000L) |>
      qdm_build(),
    ref_design_bin2("BF1")
  )
})

test_that("qdm_build() gives the designs of design_continuous2()", {
  expect_identical(
    qdm_endpoint("continuous2") |>
      qdm_sample_size(12, 12) |>
      qdm_prior(prior_niw(c(1.984, 1.191), 0.1, 4, lambda_cs4),
                prior_niw(c(-0.5, -1.5), 0.1, 4, lambda_cs4)) |>
      qdm_criterion(criterion_predictive(c(1.984, 1.191), m_t = 50, m_c = 50,
                                         gamma_go = 0.6, gamma_nogo = 0.5)) |>
      qdm_external(
        external_t = borrow(data_continuous2(8, c(2, 1),
                                             matrix(c(25, 3, 3, 9), 2L)),
                            0.5),
        external_c = borrow(data_continuous2(20, c(-0.5, -1.5),
                                             matrix(c(24, 2, 2, 10), 2L)),
                            0.4)
      ) |>
      qdm_computation(method = "MC", nMC = 500L) |>
      qdm_build(),
    ref_design_cont2("K3", method = "MC", nMC = 500L)
  )
  expect_identical(
    qdm_endpoint("continuous2") |>
      qdm_sample_size(15, 0) |>
      qdm_prior(prior_vague()) |>
      qdm_criterion(crit_cont2(0.4, 0.3, go_regions = 1:2)) |>
      qdm_control(hypothetical_c = hypothetical(mean = c(-10, -1),
                                                r = 1.5)) |>
      qdm_build(),
    ref_design_cont2("K4")
  )
  expect_identical(
    qdm_endpoint("continuous2") |>
      qdm_sample_size(10, 0) |>
      qdm_prior(prior_vague()) |>
      qdm_criterion(crit_cont2()) |>
      qdm_control(fixed_c = c(-10, -1)) |>
      qdm_build(),
    ref_design_cont2("KF1")
  )
})

test_that("qdm_build() gives the designs of design_tte() (milestone)", {
  g <- prior_gamma(0.1, 1)
  spec <- qdm_endpoint("tte") |>
    qdm_prior(g, g) |>
    qdm_criterion(crit_tte()) |>
    qdm_tte(milestone = 3, conduct = conduct_example)
  expect_identical(spec |> qdm_sample_size(20, 20) |> qdm_build(),
                   ref_design_tte("T2"))
  expect_identical(
    spec |>
      qdm_sample_size(30, 0) |>
      qdm_control(hypothetical_c = hypothetical(data_tte(8, 60))) |>
      qdm_build(),
    ref_design_tte("T5")
  )
  expect_identical(
    spec |>
      qdm_sample_size(30, 0) |>
      qdm_external(external_t = borrow(data_tte(4, 30), 0.5),
                   external_c = borrow(data_tte(20, 100), 0.5)) |>
      qdm_build(),
    ref_design_tte("T6")
  )
  expect_identical(
    spec |>
      qdm_sample_size(30, 0) |>
      qdm_prior(g) |>
      qdm_tte(milestone = 3, conduct = conduct_complete) |>
      qdm_control(fixed_c = 0.72) |>
      qdm_build(),
    ref_design_tte("T1")
  )
  expect_identical(
    spec |>
      qdm_sample_size(20, 20) |>
      qdm_criterion(criterion_predictive(0.05, m_t = 10, m_c = 8,
                                         gamma_go = 0.8, gamma_nogo = 0.3)) |>
      qdm_tte(milestone = 3, conduct = conduct_complete) |>
      qdm_computation(nMC = 500L) |>
      qdm_build(),
    ref_design_tte("T4", nMC = 500L)
  )
})

test_that("qdm_build() gives the designs of design_tte() (hazard ratio)", {
  g <- prior_gamma(0.1, 1)
  spec <- qdm_endpoint("tte") |>
    qdm_criterion(crit_hr()) |>
    qdm_tte(effect = "hazard_ratio", conduct = conduct_complete)
  expect_identical(
    spec |> qdm_sample_size(30, 30) |> qdm_prior(g, g) |> qdm_build(),
    ref_design_hr("H1")
  )
  expect_identical(
    spec |>
      qdm_sample_size(30, 0) |>
      qdm_prior(g) |>
      qdm_control(fixed_c = 6) |>
      qdm_build(),
    ref_design_hr("H2")
  )
})

test_that("qdm_build() passes the conduct of the future study", {
  g <- prior_gamma(0.1, 1)
  crit <- criterion_predictive(0.05, m_t = 10, m_c = 8, gamma_go = 0.8,
                               gamma_nogo = 0.3)
  des <- qdm_endpoint("tte") |>
    qdm_sample_size(20, 20) |>
    qdm_prior(g, g) |>
    qdm_criterion(crit) |>
    qdm_tte(milestone = 3, conduct = conduct_example,
            conduct_future = conduct_complete) |>
    qdm_computation(nMC = 500L) |>
    qdm_build()
  expect_identical(des, design_tte(20, 20, g, g, crit, milestone = 3,
                                   conduct = conduct_example,
                                   conduct_future = conduct_complete,
                                   nMC = 500L))
  expect_identical(des$conduct, conduct_example)
  expect_identical(des$conduct_future, conduct_complete)
})

test_that("the comparisons with the reference designs detect a difference", {
  spec <- qdm_endpoint("binary") |>
    qdm_sample_size(13, 12) |>
    qdm_prior(prior_beta(1, 1), prior_beta(1, 1)) |>
    qdm_criterion(crit_post())
  expect_false(identical(qdm_build(spec), ref_design("D1")))
  expect_false(identical(qdm_build(spec |> qdm_sample_size(12, 12) |>
                                     qdm_criterion(crit_post(0.6, 0.8))),
                         ref_design("D1")))
  expect_identical(qdm_build(spec |> qdm_sample_size(12, 12)),
                   ref_design("D1"))
})

test_that("qdm_build() reports the required arguments that are not set", {
  expect_error(qdm_build(qdm_endpoint("binary")),
               paste("The design specification is incomplete: set the sample",
                     "sizes (qdm_sample_size()), the prior distributions",
                     "(qdm_prior()), the decision criterion",
                     "(qdm_criterion())."),
               fixed = TRUE)
  expect_error(qdm_endpoint("binary") |>
                 qdm_sample_size(12, 12) |>
                 qdm_prior(prior_beta(1, 1)) |>
                 qdm_build(),
               paste("The design specification is incomplete: set the",
                     "decision criterion (qdm_criterion())."),
               fixed = TRUE)
})

test_that("qdm_build() gives the errors of the design_*() function", {
  spec <- qdm_endpoint("binary") |>
    qdm_sample_size(0, 10) |>
    qdm_prior(prior_beta(1, 1), prior_beta(1, 1)) |>
    qdm_criterion(crit_post())
  expect_error(qdm_build(spec), "'n_t' must be at least 1.", fixed = TRUE)
  expect_error(spec |>
                 qdm_sample_size(12, 12) |>
                 qdm_prior(prior_nix(0, 1, 3, 1), prior_beta(1, 1)) |>
                 qdm_build(),
               "'prior_t' must be created by prior_beta().", fixed = TRUE)
  expect_error(qdm_build(list()),
               "'x' must be a design specification", fixed = TRUE)
  # A design gives the same design
  expect_identical(qdm_build(ref_design("D1")), ref_design("D1"))
})
