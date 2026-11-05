# Validation of qdm_compare() against Kang, Yamaguchi and Han (2026),
# Supplemental Tables 3 and 4: the operating characteristics of a single
# binary endpoint for the sample sizes 6/6, 12/6, 12/12 and 24/24 under the
# same scenarios (posterior criterion in Table 3, posterior predictive
# criterion in Table 4). Each table is computed by one call of
# qdm_compare() with the four designs.
#
# Verdicts are those of validate_kang2026_binary() (PASS, EXPLAINED with the
# reason "kang2026_mc", FAIL), and a row is also FAIL when the probabilities
# of qdm_compare() differ from those of qdm_oc() for the same design and
# scenario (validate_kang2026_binary()), which must be identical.
#
# This file defines validate_kang2026_compare(), same_probs() and
# kang2026_compare_self_check(); run_validation.R calls them after
# validate_kang2026_binary.R.

# Whether the probabilities of two sets of rows are identical
same_probs <- function(a, b) {
  vapply(seq_len(nrow(a)), function(i) {
    identical(unlist(a[i, c("Go", "NoGo", "Gray", "Miss")], use.names = FALSE),
              unlist(b[i, c("Go", "NoGo", "Gray", "Miss")], use.names = FALSE))
  }, logical(1))
}

validate_kang2026_compare <- function(ref_file, res) {
  ref <- utils::read.csv(ref_file, stringsAsFactors = FALSE)
  out <- list()
  for (tb in c(3L, 4L)) {
    d <- ref[ref$table == tb, ]
    key <- paste(d$n_t, d$n_c, d$m_t, d$m_c, sep = "|")
    first <- d[!duplicated(key), ]
    labels <- sprintf("n_t = %d, n_c = %d", first$n_t, first$n_c)
    designs <- lapply(seq_len(nrow(first)),
                      function(i) kang2026_design(first[i, ]))
    names(designs) <- labels
    scen <- unique(d[, c("pi_t", "pi_c")])
    rownames(scen) <- NULL
    cmp <- qdm_compare(designs, scen)
    design <- labels[match(key, key[!duplicated(key)])]
    idx <- match(paste(design, d$pi_t, d$pi_c),
                 paste(cmp$design, cmp$pi_t, cmp$pi_c))
    m <- match(paste(tb, d$row, d$criterion),
               paste(res$table, res$row, res$criterion))
    r <- data.frame(
      table = tb, row = d$row, criterion = d$criterion, design = design,
      effect = d$effect, Go = cmp$Go[idx], NoGo = cmp$NoGo[idx],
      Gray = cmp$Gray[idx], Miss = cmp$Miss[idx], pub_Go = d$Go,
      pub_NoGo = d$NoGo, pub_Gray = d$Gray,
      sensitivity = res$sensitivity[m]
    )
    r$same_as_qdm_oc <- same_probs(r, res[m, ])
    v <- kang2026_verdict(r)
    v$verdict[!r$same_as_qdm_oc] <- "FAIL"
    v$reason[!r$same_as_qdm_oc] <- "differs_from_qdm_oc"
    out[[length(out) + 1L]] <- cbind(r, v)
  }
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

# Self-checks: shifting the published Go probabilities by 0.02 must give FAIL
# for all rows of the posterior criterion, and a change of 1e-12 in one
# probability of qdm_compare() must be detected as a difference from
# qdm_oc().
kang2026_compare_self_check <- function(res_cmp, res) {
  shift <- kang2026_self_check(res_cmp)
  m <- match(paste(res_cmp$table, res_cmp$row, res_cmp$criterion),
             paste(res$table, res$row, res$criterion))
  changed <- res_cmp
  changed$Go[1L] <- changed$Go[1L] + 1e-12
  detected <- !same_probs(changed, res[m, ])[1L] &&
    all(same_probs(res_cmp, res[m, ]))
  list(n_posterior = shift$n_posterior, n_fail = shift$n_fail,
       ok_shift = shift$ok, ok_identity = detected,
       ok = shift$ok && detected)
}
