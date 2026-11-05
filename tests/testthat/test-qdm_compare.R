test_that("qdm_compare() gives the results of qdm_oc() for each design", {
  d1 <- ref_design("D1")
  d2 <- ref_design("D2")
  cmp <- qdm_compare(list(a = d1, b = d2), scen_c)
  expect_s3_class(cmp, c("qdm_comparison", "data.frame"), exact = TRUE)
  expect_identical(levels(cmp$design), c("a", "b"))
  expect_identical(names(cmp), c("design", "pi_t", "pi_c", "Go", "NoGo",
                                 "Gray", "Miss"))
  designs <- list(a = d1, b = d2)
  for (lab in names(designs)) {
    oc <- qdm_oc(designs[[lab]], scen_c)
    rows <- cmp[cmp$design == lab, ]
    for (col in c("pi_t", "pi_c", "Go", "NoGo", "Gray", "Miss")) {
      expect_identical(rows[[col]], oc[[col]], label = paste(lab, col))
    }
  }
  expect_identical(attr(cmp, "designs"), designs)
  expect_identical(attr(cmp, "method"),
                   c(a = attr(qdm_oc(d1, scen_c), "method"),
                     b = attr(qdm_oc(d2, scen_c), "method")))
})

test_that("qdm_compare() names the designs and builds specifications", {
  cmp <- qdm_compare(list(spec_d1(), b = ref_design("D2"), ref_design("D1")),
                     scen_c)
  expect_identical(levels(cmp$design), c("Design 1", "b", "Design 3"))
  expect_identical(attr(cmp, "designs")[["Design 1"]], ref_design("D1"))
  expect_identical(nrow(cmp), 12L)
  expect_identical(levels(qdm_compare(list(ref_design("D1")),
                                      scen_c)$design), "Design 1")
})

test_that("qdm_compare() uses the same seed for every design", {
  d <- ref_design_cont("C1", method = "MM")
  d2 <- qdm_sample_size(d, 20, 20)
  cmp <- qdm_compare(list(d, d2), scen_cont_c, nsim = 200L, seed = 3)
  oc2 <- qdm_oc(d2, scen_cont_c, nsim = 200L, seed = 3)
  expect_identical(cmp$Go[cmp$design == "Design 2"], oc2$Go)
  expect_identical(cmp$NoGo[cmp$design == "Design 2"], oc2$NoGo)
  expect_identical(unname(attr(cmp, "method")[2L]), attr(oc2, "method"))
})

test_that("qdm_compare() names the design in the warnings", {
  expect_warning(qdm_compare(list(miss = ref_design("D6")), scen_c),
                 "^miss: The Go and NoGo criteria are both met")
  expect_no_warning(qdm_compare(list(ref_design("D1")), scen_c))
})

test_that("qdm_compare() rejects invalid arguments", {
  msg <- "'designs' must be a list of designs or design specifications."
  expect_error(qdm_compare(ref_design("D1"), scen_c), msg, fixed = TRUE)
  expect_error(qdm_compare(list(), scen_c), msg, fixed = TRUE)
  expect_error(qdm_compare(list(ref_design("D1"), 1), scen_c),
               "Each element of 'designs' must be a design", fixed = TRUE)
  expect_error(qdm_compare(list(a = ref_design("D1"), a = ref_design("D2")),
                           scen_c),
               "The names of the designs must be unique.", fixed = TRUE)
  expect_error(qdm_compare(list(ref_design("D1"), ref_design_cont("C1")),
                           scen_c),
               "The designs must have the same endpoint.", fixed = TRUE)
  expect_error(qdm_compare(list(ref_design("D1")), cbind(scen_c, design = 1)),
               "'scenarios' must not have a column named design.",
               fixed = TRUE)
})
