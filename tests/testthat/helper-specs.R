# Design specification built step by step that gives the design D1 of
# helper-designs.R, used by the tests of the functions that accept a
# specification in place of a design.
spec_d1 <- function() {
  qdm_endpoint("binary") |>
    qdm_sample_size(12, 12) |>
    qdm_prior(prior_beta(1, 1), prior_beta(1, 1)) |>
    qdm_criterion(crit_post())
}
