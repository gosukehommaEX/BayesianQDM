test_that("report_symbol() gives the symbols of the columns of the scenarios", {
  expect_identical(report_symbol(c("pi_t", "sd_c", "rho_t", "median_t")),
                   c("\\pi_t", "\\sigma_c", "\\rho_t", "M_t"))
  expect_identical(report_symbol(c("mu_t1", "label")),
                   c("\\mu_{t1}", NA_character_))
})

test_that("report_symbol() has a symbol for every column of every endpoint", {
  for (des in list(ref_design("D1"), ref_design_cont("C1"),
                   ref_design_bin2("B4"), ref_design_cont2("K1"),
                   ref_design_tte("T2"))) {
    cols <- scenario_columns(des)$names
    expect_false(anyNA(report_symbol(cols)), info = des$endpoint)
  }
})
