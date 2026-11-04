test_that("html_escape() replaces the characters of HTML markup", {
  expect_identical(html_escape("a < b & c > d"), "a &lt; b &amp; c &gt; d")
  expect_identical(html_escape(c("&lt;", "x")), c("&amp;lt;", "x"))
  expect_identical(html_escape("$\\theta \\le 0.1$"), "$\\theta \\le 0.1$")
})
