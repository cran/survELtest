#' @importFrom stats quantile
neg2ELratio_nocrossings <- function(alpha, fit, fit1, fit2, fit_time_restrict_boot, nn, lowerbindx_boot, upperbindx_boot, sum_DWknGImuw_1_big, sum_DWknGImuw_2_big){
  # function for computation for nocrossings
  test_nocross = 0
  
  HW_CBdistr = apply(abs(sum_DWknGImuw_2_big[,lowerbindx_boot:upperbindx_boot]-sum_DWknGImuw_1_big[,lowerbindx_boot:upperbindx_boot]),1,max)
  HW_CBcrit = as.vector(quantile(HW_CBdistr,1-alpha))
  S1_hat = (((c(1,fit1$surv)[cumsum(c(0,fit$time) %in% c(0,fit1$time))])[-1])[fit_time_restrict_boot])[lowerbindx_boot:upperbindx_boot] #except for survival, other "expansions" can't use this
  S2_hat = (((c(1,fit2$surv)[cumsum(c(0,fit$time) %in% c(0,fit2$time))])[-1])[fit_time_restrict_boot])[lowerbindx_boot:upperbindx_boot]
  diff = log(S1_hat)-log(S2_hat)
  HW_CB_ubs = diff+HW_CBcrit/sqrt(sum(nn))
  HW_CB_lbs = diff-HW_CBcrit/sqrt(sum(nn))
  eq_S_all = sum(HW_CB_lbs <= 0) == (upperbindx_boot - lowerbindx_boot + 1) & sum(HW_CB_ubs >= 0) == (upperbindx_boot - lowerbindx_boot + 1)
  g_S_all = sum(HW_CB_ubs >= 0) == (upperbindx_boot - lowerbindx_boot + 1) & sum(HW_CB_lbs <= 0) != (upperbindx_boot - lowerbindx_boot + 1)
  
  if (g_S_all | eq_S_all) {
    test_nocross = 1
  }
  return(list(test_nocross = test_nocross))
}
