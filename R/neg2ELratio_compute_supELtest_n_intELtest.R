
neg2ELratio_compute_supELtest_n_intELtest <- function(fit, fit_time_restrict_boot, lowerbindx_boot, upperbindx_boot, Td_sort_boot, barNt, teststat_pre){
  # function for computing supELtest and intELtest statistics based on pointwise EL statistics
  inttest_pre_dF     <- teststat_pre*diff(c(0, 1 - fit$surv[fit_time_restrict_boot]))[lowerbindx_boot:upperbindx_boot]
  inttest_pre_dbarNt <- teststat_pre*diff(c(0, barNt))[lowerbindx_boot:upperbindx_boot]
  inttest_pre_dt     <- teststat_pre[-(upperbindx_boot - lowerbindx_boot + 1)]*diff(Td_sort_boot)
  
  suptest        <- max(teststat_pre)
  inttest_dF     <- sum(inttest_pre_dF)
  inttest_dbarNt <- sum(inttest_pre_dbarNt)
  inttest_dt     <- sum(inttest_pre_dt)
  return(list(suptest = suptest, inttest_dF = inttest_dF, inttest_dbarNt = inttest_dbarNt, inttest_dt = inttest_dt))
}
