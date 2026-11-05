test_that("show_code() gives code that recreates every design", {
  designs <- all_ref_designs()
  for (nm in names(designs)) {
    des <- designs[[nm]]
    for (style in c("pipe", "design")) {
      code <- show_code(des, style = style, name = NULL)
      expect_identical(eval(parse(text = code)), des,
                       label = paste(nm, style))
    }
  }
})

test_that("show_code() writes the pipe of the design steps", {
  # Expected lines: Python_reference/port_code_format.out
  code <- show_code(ref_design("D2"))
  expect_s3_class(code, "qdm_code")
  expect_identical(unclass(code), c(
    "design <- qdm_endpoint(\"binary\") |>",
    "  qdm_sample_size(n_t = 10, n_c = 8) |>",
    "  qdm_prior(prior_t = prior_beta(a = 1, b = 1),",
    "            prior_c = prior_beta(a = 2, b = 3)) |>",
    "  qdm_criterion(criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,",
    "                                    gamma_go = 0.5, gamma_nogo = 0.8)) |>",
    "  qdm_external(external_t = borrow(data = data_binary(n = 10, y = 6),",
    "                                   a0 = 0.5),",
    "               external_c = borrow(data = data_binary(n = 20, y = 8),",
    "                                   a0 = 0.5)) |>",
    "  qdm_build()"
  ))
  expect_true(all(nchar(code) <= 80L))
  expect_identical(unclass(show_code(ref_design_tte("T2"), name = NULL)), c(
    "qdm_endpoint(\"tte\") |>",
    "  qdm_sample_size(n_t = 20, n_c = 20) |>",
    "  qdm_prior(prior_t = prior_gamma(shape = 0.1, rate = 1),",
    "            prior_c = prior_gamma(shape = 0.1, rate = 1)) |>",
    "  qdm_criterion(criterion_posterior(theta_tv = 0.15, theta_mav = 0,",
    "                                    gamma_go = 0.6, gamma_nogo = 0.8)) |>",
    "  qdm_tte(milestone = 3,",
    "          conduct = conduct_tte(accrual_time = c(0, 5, 10, 15),",
    "                                accrual_rate = c(0.4, 2, 4.5),",
    paste0("                                dropout_hazard = 0.00226, ",
           "cut_time = 15.44)) |>"),
    "  qdm_build()"
  ))
})

test_that("show_code() writes the call of the design function", {
  # Expected lines: Python_reference/port_code_format.out
  expect_identical(unclass(show_code(ref_design("D2"), style = "design")), c(
    "design <- design_binary(",
    "  n_t = 10, n_c = 8, prior_t = prior_beta(a = 1, b = 1),",
    "  prior_c = prior_beta(a = 2, b = 3),",
    "  criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,",
    "                                  gamma_go = 0.5, gamma_nogo = 0.8),",
    "  external_t = borrow(data = data_binary(n = 10, y = 6), a0 = 0.5),",
    "  external_c = borrow(data = data_binary(n = 20, y = 8), a0 = 0.5)",
    ")"
  ))
  # Expected lines: Python_reference/port_code_format.out
  des <- design_tte(30, 0, prior_gamma(1, 1), criterion = crit_tte(),
                    milestone = 3, fixed_c = 0.7)
  expect_identical(unclass(show_code(des, style = "design", name = NULL)), c(
    "design_tte(",
    "  n_t = 30, n_c = 0, prior_t = prior_gamma(shape = 1, rate = 1),",
    "  criterion = criterion_posterior(theta_tv = 0.15, theta_mav = 0,",
    "                                  gamma_go = 0.6, gamma_nogo = 0.8),",
    "  milestone = 3, fixed_c = 0.7",
    ")"
  ))
  # The name starts the first line; a short call is on one line
  des2 <- design_binary(5, 0, prior_beta(1, 1), criterion = crit_post(),
                        fixed_c = 0.2)
  expect_identical(show_code(des2, style = "design", name = "d")[1L],
                   "d <- design_binary(")
  spec <- qdm_endpoint("binary") |> qdm_sample_size(5, 5)
  expect_identical(unclass(show_code(spec, style = "design")),
                   "design <- design_binary(n_t = 5, n_c = 5)")
})

test_that("show_code() leaves out the arguments equal to their defaults", {
  code <- paste(show_code(ref_design_cont("C1")), collapse = "\n")
  expect_false(grepl("method|nMC|qdm_computation", code))
  code_mm <- show_code(ref_design_cont("C1", method = "MM"))
  expect_true(any(grepl("qdm_computation(method = \"MM\")", code_mm,
                        fixed = TRUE)))
  # The future conduct is the conduct of the study
  code_t4 <- paste(show_code(ref_design_tte("T4")), collapse = "\n")
  expect_false(grepl("conduct_future", code_t4))
  expect_true(grepl("cut_time = 1000000", code_t4, fixed = TRUE))
})

test_that("show_code() writes specifications and the designs of results", {
  code <- show_code(spec_d1(), name = NULL)
  expect_false(any(grepl("qdm_build", code)))
  spec <- eval(parse(text = code))
  expect_s3_class(spec, "qdm_spec")
  expect_identical(qdm_build(spec), ref_design("D1"))
  des <- ref_design("D1")
  expect_identical(show_code(qdm_oc(des, scen_c)), show_code(des))
  expect_identical(show_code(qdm_boundary(des)), show_code(des))
  expect_identical(show_code(qdm_decide(des, data_binary(12, 9),
                                        data_binary(12, 4))),
                   show_code(des))
  cal <- run_cal("D1")
  expect_identical(show_code(cal), show_code(cal$design))
  expect_true(any(grepl("gamma_go = 0.7", show_code(cal), fixed = TRUE)))
  res <- find_n_d1(n_t = seq(6, 20, by = 2))
  expect_identical(show_code(res), show_code(res$design))
})

test_that("show_code() assigns the code to a name", {
  des <- ref_design("D1")
  expect_true(startsWith(show_code(des)[1L], "design <- "))
  expect_true(startsWith(show_code(des, name = "pilot")[1L], "pilot <- "))
  expect_false(grepl("<-", show_code(des, name = NULL)[1L], fixed = TRUE))
  expect_error(show_code(des, name = "two words"),
               "'name' must be NULL or a syntactic name", fixed = TRUE)
  expect_error(show_code(des, name = c("a", "b")),
               "'name' must be NULL or a syntactic name", fixed = TRUE)
  expect_error(show_code(des, style = "other"), "'arg' should be one of",
               fixed = TRUE)
})
