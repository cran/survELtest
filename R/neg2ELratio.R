#' @importFrom survival Surv
neg2ELratio <- function(formula, data, group_order, t1, t2, sided, nboot, alpha, details.return, seed, nlimit) {
  set.seed(seed)
  
  list_t1t2 <- neg2ELratio_t1t2(formula, group_order, t1, t2, nboot)
  list_bootstrap <- neg2ELratio_bootstrap(nboot, nlimit, list_t1t2$fit, list_t1t2$fit1, list_t1t2$fit2, list_t1t2$NBOOT, list_t1t2$fit_time_restrict_boot, list_t1t2$mm, list_t1t2$nn)
  list_bootstrap_result <- neg2ELratio_bootstrap_results(sided, alpha, list_t1t2$fit, list_t1t2$fit1, list_t1t2$fit2, list_t1t2$NBOOT, list_t1t2$fit_time_restrict_boot, list_t1t2$mm, list_t1t2$nn, list_t1t2$lowerbindx_boot,  list_t1t2$upperbindx_boot, list_bootstrap$U_boot_H1)
  gc()
  list_EL <- neg2ELratio_EL(sided, list_t1t2$fit1, list_t1t2$fit2, list_t1t2$T_11, list_t1t2$T_21, list_t1t2$upperbindx_boot, list_t1t2$lowerbindx_boot, list_t1t2$Td_sort_boot, list_bootstrap_result$EL_CBcrit)
  list_compute_supELtest_n_intELtest <- neg2ELratio_compute_supELtest_n_intELtest(list_t1t2$fit, list_t1t2$fit_time_restrict_boot, list_t1t2$lowerbindx_boot, list_t1t2$upperbindx_boot, list_t1t2$Td_sort_boot, list_bootstrap_result$barNt, list_EL$teststat_pre)
  list_nocrossings <- neg2ELratio_nocrossings(alpha, list_t1t2$fit, list_t1t2$fit1, list_t1t2$fit2, list_t1t2$fit_time_restrict_boot, list_t1t2$nn, list_t1t2$lowerbindx_boot, list_t1t2$upperbindx_boot, list_bootstrap$sum_DWknGImuw_1_big, list_bootstrap$sum_DWknGImuw_2_big)
  
  if (details.return == TRUE) {
    return (list(
      test_nocross = list_nocrossings$test_nocross,
      neg2ELratio_at_ts = list_EL$teststat_pre,
      wt_dbarNt = list_bootstrap_result$wt_dbarNt,
      wt_dt = list_bootstrap_result$wt_dt,
      wt_dF = list_bootstrap_result$wt_dF,
      neg2ELratio_bootstrap_at_ts = list_bootstrap_result$Up2_boot_H1_1sided,
      wt_dbarNt_boot = list_bootstrap_result$wt_dbarNt_boot,
      wt_dt_boot = list_bootstrap_result$wt_dt_boot,
      wt_dF_boot = list_bootstrap_result$wt_dF_boot,
      lowerbindx_boot = list_t1t2$lowerbindx_boot,
      upperbindx_boot = list_t1t2$upperbindx_boot,
      ts = list_t1t2$Td_sort_boot,
      bootstrap_ts = list_t1t2$fit$time[list_t1t2$fit_time_restrict_boot],

      EL_SOcrit           = list_bootstrap_result$EL_SOcrit,
      int_dFEL_SOcrit     = list_bootstrap_result$int_dFEL_SOcrit,
      int_dbarNtEL_SOcrit = list_bootstrap_result$int_dbarNtEL_SOcrit,
      int_dtEL_SOcrit     = list_bootstrap_result$int_dtEL_SOcrit,
      
      suptest        = list_compute_supELtest_n_intELtest$suptest,
      inttest_dF     = list_compute_supELtest_n_intELtest$inttest_dF,
      inttest_dbarNt = list_compute_supELtest_n_intELtest$inttest_dbarNt,
      inttest_dt     = list_compute_supELtest_n_intELtest$inttest_dt,
      
      p_value_suptest        = mean(list_bootstrap_result$sup_boot_H1        > list_compute_supELtest_n_intELtest$suptest),
      p_value_inttest_dF     = mean(list_bootstrap_result$int_dF_boot_H1     > list_compute_supELtest_n_intELtest$inttest_dF),
      p_value_inttest_dbarNt = mean(list_bootstrap_result$int_dbarNt_boot_H1 > list_compute_supELtest_n_intELtest$inttest_dbarNt),
      p_value_inttest_dt     = mean(list_bootstrap_result$int_dt_boot_H1     > list_compute_supELtest_n_intELtest$inttest_dt)
    ))
  } else return (list(neg2ELratio_at_ts = list_EL$teststat_pre))
}
