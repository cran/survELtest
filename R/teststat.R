teststat = function(formula, data, group_order, t1, t2, sided, nboot, alpha, details.return, seed, nlimit) {
  set.seed(seed)
  
  list_t1t2 = teststat_t1t2(formula, group_order, t1, t2, nboot)
  list_bootstrap = teststat_bootstrap(nlimit, list_t1t2$k, list_t1t2$iter, list_t1t2$fit, list_t1t2$fitls, list_t1t2$NBOOT, list_t1t2$mm, list_t1t2$nn, list_t1t2$fit_time_restrict_boot)
  list_bootstrap_results = teststat_bootstrap_results(sided, alpha, list_t1t2$k, list_t1t2$M_vec, list_t1t2$iter, list_t1t2$fit, list_t1t2$fitls, list_t1t2$NBOOT, list_t1t2$mm, list_t1t2$nn, list_t1t2$fit_time_restrict_boot, list_t1t2$lowerbindx_boot, list_t1t2$upperbindx_boot, list_t1t2$Td_sort_boot, list_bootstrap$sum_DWknGImuw_big)
  gc()
  list_compute_supELtest_n_intELtest = teststat_compute_supELtest_n_intELtest(sided, list_t1t2$k, list_t1t2$M_vec, list_t1t2$iter, list_t1t2$iter_without_last, list_t1t2$fitls, list_t1t2$lowerbindx_boot, list_t1t2$upperbindx_boot, list_t1t2$Td_sort_boot, list_bootstrap_results$wt_dbarNt, list_bootstrap_results$wt_dt, list_bootstrap_results$wt_dF)
  list_nocrossings = teststat_nocrossings(alpha, list_t1t2$k, list_t1t2$M_vec, list_t1t2$iter, list_t1t2$iter_without_last, list_t1t2$fit, list_t1t2$fitls, list_t1t2$nn, list_t1t2$fit_time_restrict_boot, list_t1t2$lowerbindx_boot, list_t1t2$upperbindx_boot, list_bootstrap$sum_DWknGImuw_big)
  
  if (details.return == TRUE) {
    return (list(
      test_nocross = list_nocrossings$test_nocross,
      neg2ELratio_at_ts = list_compute_supELtest_n_intELtest$teststat_pre,
      wt_dbarNt = list_bootstrap_results$wt_dbarNt,
      wt_dt = list_bootstrap_results$wt_dt,
      wt_dF = list_bootstrap_results$wt_dF,
      neg2ELratio_bootstrap_at_ts = list_bootstrap_results$Up2_boot_H1_1sided,
      wt_dbarNt_boot = list_bootstrap_results$wt_dbarNt_boot,
      wt_dt_boot = list_bootstrap_results$wt_dt_boot,
      wt_dF_boot = list_bootstrap_results$wt_dF_boot,
      lowerbindx_boot = list_t1t2$lowerbindx_boot,
      upperbindx_boot = list_t1t2$upperbindx_boot,
      ts = list_t1t2$Td_sort_boot,
      bootstrap_ts = list_t1t2$fit$time[list_t1t2$fit_time_restrict_boot],
      
      EL_SOcrit           = list_bootstrap_results$EL_SOcrit,
      int_dFEL_SOcrit     = list_bootstrap_results$int_dFEL_SOcrit,
      int_dbarNtEL_SOcrit = list_bootstrap_results$int_dbarNtEL_SOcrit,
      int_dtEL_SOcrit     = list_bootstrap_results$int_dtEL_SOcrit,
      
      suptest        = list_compute_supELtest_n_intELtest$suptest,
      inttest_dF     = list_compute_supELtest_n_intELtest$inttest_dF,
      inttest_dbarNt = list_compute_supELtest_n_intELtest$inttest_dbarNt,
      inttest_dt     = list_compute_supELtest_n_intELtest$inttest_dt,
      
      p_value_suptest        = mean(list_bootstrap_results$sup_boot_H1        > list_compute_supELtest_n_intELtest$suptest), # same as out_sup_pval = mean(sup_boot_H1 >= max(teststat_pre)) in original teststat
      p_value_inttest_dF     = mean(list_bootstrap_results$int_dF_boot_H1     > list_compute_supELtest_n_intELtest$inttest_dF), # same as  = mean(int_dF_boot_H1 >= sum(inttest_pre_dF)) in original teststat
      p_value_inttest_dbarNt = mean(list_bootstrap_results$int_dbarNt_boot_H1 > list_compute_supELtest_n_intELtest$inttest_dbarNt),
      p_value_inttest_dt     = mean(list_bootstrap_results$int_dt_boot_H1     > list_compute_supELtest_n_intELtest$inttest_dt)
    ))
  } else return (list(neg2ELratio_at_ts = list_compute_supELtest_n_intELtest$teststat_pre))
}
