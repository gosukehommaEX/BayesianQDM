# All the designs of the helper files, named as in the helpers, for the
# tests that loop over every kind of design.
all_ref_designs <- function() {
  c(lapply(stats::setNames(nm = paste0("D", 1:9)), ref_design),
    lapply(stats::setNames(nm = paste0("C", 1:8)), ref_design_cont),
    lapply(stats::setNames(nm = c(paste0("B", 1:8), "BF1", "BF2")),
           ref_design_bin2),
    lapply(stats::setNames(nm = c(paste0("K", 1:6), "KF1", "KF2")),
           ref_design_cont2),
    lapply(stats::setNames(nm = paste0("T", 1:6)), ref_design_tte),
    lapply(stats::setNames(nm = c(paste0("H", 1:4), "P1", "P2")),
           ref_design_hr))
}
