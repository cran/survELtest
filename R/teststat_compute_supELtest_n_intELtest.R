
teststat_compute_supELtest_n_intELtest <- function(sided, k, M_vec, iter, iter_without_last, fitls, lowerbindx_boot, upperbindx_boot, Td_sort_boot, wt_dbarNt, wt_dt, wt_dF){
  # function for computing supELtest and intELtest statistics based on pointwise EL statistics
  ### computing quantities related to the NPLR statistics:
  teststat_pre = 1:length(Td_sort_boot) * 0
  error_vec = 1:length(Td_sort_boot) * 0
  for (j in 1:(upperbindx_boot - lowerbindx_boot + 1)) {
    t = Td_sort_boot[j]  # the given t >= 0 in localized EL statistic
    neg2logRt = teststat_neg2logR(sided, k, M_vec, iter, iter_without_last, fitls, t, EL_CBcrit = 0)
    error_vec[j] = neg2logRt$still_error
    if (error_vec[j] == 1) next
    teststat_pre[j] = neg2logRt$out  # the -2 log NPLR statistic in Section 2.2 of the paper at the given t >=0 specified above
  }
  
  inttest_pre_dF = teststat_pre * wt_dF[lowerbindx_boot:upperbindx_boot]  # -2logR(t) d\hat{F}_0(t) in (3.1) of the paper
  inttest_pre_dbarNt = teststat_pre * wt_dbarNt[lowerbindx_boot:upperbindx_boot]
  inttest_pre_dt     = (teststat_pre * wt_dt[lowerbindx_boot:upperbindx_boot])[-(upperbindx_boot - lowerbindx_boot + 1)]
  
  suptest        = max(teststat_pre) # same as test = max(teststat_pre) in original teststat
  inttest_dF     = sum(inttest_pre_dF)
  inttest_dbarNt = sum(inttest_pre_dbarNt)
  inttest_dt     = sum(inttest_pre_dt)
  
  return(list(teststat_pre = teststat_pre, suptest = suptest, inttest_dF = inttest_dF, inttest_dbarNt = inttest_dbarNt, inttest_dt = inttest_dt))
}

