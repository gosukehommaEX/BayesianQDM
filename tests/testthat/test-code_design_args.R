test_that("code_design_args() leaves out the arguments equal to defaults", {
  args <- code_design_args(design_spec(ref_design_cont("C1")))
  expect_identical(names(args),
                   c("n_t", "n_c", "prior_t", "prior_c", "criterion"))
  args <- code_design_args(design_spec(ref_design_cont("C1", method = "MC",
                                                       nMC = 500L)))
  expect_identical(args$method, "MC")
  expect_identical(args$nMC, 500L)
  expect_identical(names(code_design_args(design_spec(ref_design_hr("H1")))),
                   c("n_t", "n_c", "prior_t", "prior_c", "criterion",
                     "effect", "conduct"))
  # The default of conduct_future is the value of conduct
  expect_false("conduct_future" %in%
                 names(code_design_args(design_spec(ref_design_tte("T4")))))
  spec <- qdm_endpoint("tte") |>
    qdm_tte(milestone = 3, conduct = conduct_example,
            conduct_future = conduct_complete)
  expect_identical(code_design_args(spec)$conduct_future, conduct_complete)
  # Arguments of a specification that are not yet converted are kept
  spec2 <- qdm_endpoint("binary2") |> qdm_computation(nMC = 10000)
  expect_identical(code_design_args(spec2)$nMC, 10000)
})
