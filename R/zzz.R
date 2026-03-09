# Suppress R CMD check NOTEs for ggplot2 aes() column name references.
# These variables are data frame column names used inside aes() and are not
# true global variables; utils::globalVariables() informs the static analyser.
utils::globalVariables(c(
  "Decision",   # plot.pbayesdecisionprob1bin, plot.pbayesdecisionprob1cont
  "prob",       # plot.getgamma1bin, plot.getgamma1cont,
  # plot.getgamma2bin, plot.getgamma2cont
  "prob_val",   # plot.pbayesdecisionprob2bin, plot.pbayesdecisionprob2cont
  "pt_grp",     # plot.getgamma1bin, plot.getgamma1cont,
  # plot.getgamma2bin, plot.getgamma2cont
  "Scenario"    # plot.getgamma1bin, plot.getgamma1cont,
  # plot.getgamma2bin, plot.getgamma2cont
))
