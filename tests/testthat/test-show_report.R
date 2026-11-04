report_note <- function() {
  sprintf("Report written by BayesianQDM %s.",
          utils::packageVersion("BayesianQDM"))
}

test_that("show_report() writes a text report and returns the file name", {
  f <- tempfile(fileext = ".txt")
  on.exit(unlink(f), add = TRUE)
  oc <- qdm_oc(ref_design("D1"), scen_c)
  expect_invisible(show_report(oc, file_type = "text", file = f,
                               open = FALSE))
  path <- spec_d1() |>
    qdm_oc(scen_c) |>
    show_report(file_type = "text", file = f, open = FALSE)
  expect_identical(path, f)
  expect_identical(readLines(f),
                   report_text("Decision criteria and operating characteristics",
                               report_blocks(oc, 3L), report_note()))
})

test_that("show_report() writes an HTML report", {
  f <- tempfile(fileext = ".html")
  on.exit(unlink(f), add = TRUE)
  show_report(ref_design("D1"), file = f, open = FALSE)
  expect_identical(readLines(f),
                   report_html("Decision criteria",
                               report_blocks(ref_design("D1"), 3L),
                               report_note()))
})

test_that("show_report() uses a temporary file by default", {
  html <- show_report(ref_design("D1"), open = FALSE)
  text <- show_report(ref_design("D1"), file_type = "text", open = FALSE)
  on.exit(unlink(c(html, text)), add = TRUE)
  expect_true(file.exists(html))
  expect_true(file.exists(text))
  expect_true(startsWith(basename(html), "qdm_report_"))
  expect_true(endsWith(html, ".html"))
  expect_true(endsWith(text, ".txt"))
  expect_identical(normalizePath(dirname(html)), normalizePath(tempdir()))
})

test_that("show_report() uses the title and the number of decimal places", {
  f <- tempfile(fileext = ".txt")
  on.exit(unlink(f), add = TRUE)
  oc <- qdm_oc(ref_design("D1"), scen_c)
  show_report(oc, file_type = "text", file = f, title = "My title",
              digits = 2L, open = FALSE)
  lines <- readLines(f)
  expect_identical(lines[1L], "# My title")
  expect_identical(lines, report_text("My title", report_blocks(oc, 2L),
                                      report_note()))
  expect_false(identical(lines, report_text("My title", report_blocks(oc, 3L),
                                            report_note())))
})

test_that("show_report() prints a text report and opens an HTML report", {
  f <- tempfile(fileext = ".txt")
  on.exit(unlink(f), add = TRUE)
  expect_output(show_report(ref_design("D1"), file_type = "text", file = f,
                            open = TRUE),
                "# Decision criteria", fixed = TRUE)
  opened <- NULL
  local_mocked_bindings(open_html = function(file) {
    opened <<- file
    invisible(file)
  })
  g <- tempfile(fileext = ".html")
  on.exit(unlink(g), add = TRUE)
  expect_silent(show_report(ref_design("D1"), file = g, open = TRUE))
  expect_identical(opened, g)
})

test_that("show_report() checks its arguments", {
  des <- ref_design("D1")
  expect_error(show_report(des, file_type = "pdf", open = FALSE),
               "should be one of")
  expect_error(show_report(des, file = 1, open = FALSE),
               "'file' must be NULL or the name of a file.", fixed = TRUE)
  expect_error(show_report(des, file = "", open = FALSE),
               "'file' must be NULL or the name of a file.", fixed = TRUE)
  expect_error(show_report(des, title = c("a", "b"), open = FALSE),
               "'title' must be NULL or a character string.", fixed = TRUE)
  expect_error(show_report(des, digits = 0, open = FALSE),
               "'digits' must be at least 1 and at most 10.", fixed = TRUE)
  expect_error(show_report(des, digits = 1.5, open = FALSE),
               "'digits' must contain whole numbers only.", fixed = TRUE)
  expect_error(show_report(des, open = NA), "'open' must be TRUE or FALSE.",
               fixed = TRUE)
  expect_error(show_report(list(), open = FALSE),
               "'x' must be a design, a design specification", fixed = TRUE)
})
