test_that("report_study() describes the endpoints and the sample sizes", {
  expect_identical(
    report_study(ref_design("D1")),
    list(list(type = "heading", text = "Study design"),
         list(type = "par",
              text = paste("The proof-of-concept study has one binary",
                           "endpoint and will enroll 12 subjects in the",
                           "treatment group and 12 subjects in the control",
                           "group.")))
  )
  expect_identical(
    report_study(ref_design_cont2("K4"))[[2L]][["text"]],
    paste("The proof-of-concept study has two continuous endpoints and will",
          "enroll 15 subjects in the treatment group, without a concurrent",
          "control group.")
  )
})

test_that("report_study() names each type of endpoint", {
  text <- function(des) report_study(des)[[2L]][["text"]]
  expect_match(text(ref_design_cont("C1")), "has one continuous endpoint",
               fixed = TRUE)
  expect_match(text(ref_design_bin2("B6")), "has two binary endpoints",
               fixed = TRUE)
  expect_match(text(ref_design_tte("T2")), "has one time-to-event endpoint",
               fixed = TRUE)
})
