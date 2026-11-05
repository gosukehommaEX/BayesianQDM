test_that("design_spec() gives the specification that builds the design", {
  designs <- all_ref_designs()
  expect_length(designs, 47L)
  for (nm in names(designs)) {
    des <- designs[[nm]]
    spec <- design_spec(des)
    expect_s3_class(spec, "qdm_spec")
    expect_identical(spec$endpoint, des$endpoint)
    expect_identical(qdm_build(spec), des, label = nm)
  }
})

test_that("design_spec() sets the arguments that are not NULL", {
  spec <- design_spec(ref_design("D2"))
  expect_setequal(names(spec$args),
                  c("n_t", "n_c", "prior_t", "prior_c", "criterion",
                    "external_t", "external_c"))
  expect_identical(spec$args$external_c, borrow(data_binary(20, 8), 0.5))
  expect_false("hypothetical_c" %in% names(design_spec(ref_design("D1"))$args))
  # The future conduct is kept only for a predictive criterion
  expect_false("conduct_future" %in%
                 names(design_spec(ref_design_tte("T2"))$args))
  expect_identical(design_spec(ref_design_tte("T4"))$args$conduct_future,
                   conduct_complete)
  # A design that differs in one argument does not give the reference
  other <- ref_design("D2")
  other$external_c <- borrow(data_binary(20, 8), 0.4)
  expect_false(identical(qdm_build(design_spec(other)), ref_design("D2")))
})
