# Avoid NOTES about undefined globals from data.table's NSE
utils::globalVariables(c(
  '.', "r", "tau_r", "tau_wt", "x", "bs_id", "xbin", "ape",
  "signal", "prior", "posterior", "Y", "Y0", "dX", "dY",
  "alpha", "alpha_rank", "Z", "rowid", "LLS.INT.R", "lb", "ub", "se", "cluster_id", "bsweights"
))
