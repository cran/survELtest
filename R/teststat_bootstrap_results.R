#' @importFrom stats quantile
#' @importFrom Iso "pava"
#' @importFrom plyr "alply"
teststat_bootstrap_results <- function(sided, alpha, k, M_vec, iter, fit, fitls, NBOOT, mm, nn, fit_time_restrict_boot, lowerbindx_boot, upperbindx_boot, Td_sort_boot, sum_DWknGImuw_big){
  # function for computing bootstrap results for supELtest and intELtest
  ### bootstrap components in the limiting distribution in Theorem 1:
  theta_hat_js = theta_hat_j_k(Td_sort_boot, fit, fitls, M_vec)
  sigma2_hat_overpjs = sigma2_hat_overpj_k(Td_sort_boot, fit, fitls)
  Dsqsigmadp_big = lapply(iter, FUN = function(iter){
    matrix(rep(1 / sqrt(sigma2_hat_overpjs[, iter]), times = NBOOT), byrow = TRUE, nrow = NBOOT)
  })
  Ujps = array(0, c(k, NBOOT, upperbindx_boot - lowerbindx_boot + 1))
  for (i in iter){
    Ujps[i, , ] = product_mat(as.matrix(-sum_DWknGImuw_big[[i]][, lowerbindx_boot:upperbindx_boot]), Dsqsigmadp_big[[i]])
  }
  wjs = array(0, c(k, NBOOT, upperbindx_boot - lowerbindx_boot + 1))
  for (j in iter){
    wjs[j, , ] = matrix(rep(division00(1, theta_hat_js[, j] ), times = NBOOT), byrow = TRUE, nrow = NBOOT)
  }
  wjs_denom = 0
  for (j in iter){
    wjs_denom = wjs_denom + wjs[j,,]
  }
  for (j in iter){
    wjs[j, , ] =  wjs[j, , ] / wjs_denom
  }
  if (sided == 1){
    pava_time1 = Sys.time()
    pava_result = array(unlist(mapply(pava, alply(division00(Ujps, sqrt(wjs)), c(2, 3)), alply(wjs, c(2, 3)), decreasing = TRUE)), c(k, NBOOT, upperbindx_boot - lowerbindx_boot + 1))
    pava_time2 = Sys.time()
  } else {
    pava_result = division00(Ujps, sqrt(wjs))
  }
  avg_Ujps = apply(sqrt(wjs) * Ujps, c(2, 3), sum)
  avg_Ujps_karray = array(rep(avg_Ujps, each = k), c(k, NBOOT, upperbindx_boot - lowerbindx_boot + 1))
  Up2_boot_H1_1sided = apply(wjs * (pava_result - avg_Ujps_karray) ^ 2, c(2, 3), sum)  # bootstrap SSB(t) in Theorem 1
  
  ### bootstrap sup_{t \in [t_1, t_2]} SSB(t) in Theorem 1:
  sup_boot_H1 = apply(as.matrix(Up2_boot_H1_1sided), 1, max)  # bootstrap sup_{t \in [t_1, t_2]} SSB(t) in Theorem 1
  EL_SOcrit = as.vector(quantile(sup_boot_H1, probs = 1 - alpha))  # quantiles based on bootstrapped values of sup_{t \in [t_1, t_2]} SSB(t)
  
  nobd_pooled         = length(rep(fit$time, times = fit$n.event))
  data_big_pooled     = matrix(rep(rep(fit$time, times = fit$n.event), times = mm), byrow = TRUE, nrow = mm)
  fit_time_big_pooled = matrix(rep(fit$time[fit_time_restrict_boot], times=nobd_pooled), byrow = FALSE, ncol = nobd_pooled)
  Indt_big_pooled     = (data_big_pooled <= fit_time_big_pooled)
  barNt        = apply(Indt_big_pooled, 1, sum) / sum(nn)
  wt_dbarNt    = diff(c(0, barNt))
  
  barNt_big      = matrix(rep(barNt, time = NBOOT), byrow = TRUE, nrow = NBOOT)
  wt_dbarNt_boot = t(apply(cbind(0, barNt_big), 1, diff))
  
  wt_dt      = diff(c(fit$time[fit_time_restrict_boot], Inf))  # we don't want the last one
  t_big      = matrix(rep(fit$time[fit_time_restrict_boot], time = NBOOT), byrow = TRUE, nrow = NBOOT)
  wt_dt_boot = t(apply(cbind(t_big, Inf), 1, diff))
  
  wt_dF      = diff(c(0, 1 - fit$surv[fit_time_restrict_boot]))
  F_big      = matrix(rep(1 - fit$surv[fit_time_restrict_boot], time = NBOOT), byrow = TRUE, nrow = NBOOT)
  wt_dF_boot = t(apply(cbind(0, F_big), 1, diff))
  Up2_boot_H1_1sided_times_dF = Up2_boot_H1_1sided * t(apply(cbind(0,  F_big), 1, diff))[, lowerbindx_boot:upperbindx_boot]
  int_dF_boot_H1 = apply(as.matrix(Up2_boot_H1_1sided_times_dF), 1, sum)  # bootstrap \int_{t_1}^{t_2} SSB(t) dF_0(t) in Theorem 1
  int_dFEL_SOcrit = as.vector(quantile(int_dF_boot_H1, probs = 1 - alpha))  # quantiles based on bootstrapped values of \int_{t_1}^{t_2} SSB(t) dF_0(t) in Theorem 1
  
  Up2_boot_H1_1sided_times_dbarNt = Up2_boot_H1_1sided*wt_dbarNt_boot[, lowerbindx_boot:upperbindx_boot]
  int_dbarNt_boot_H1  = apply(as.matrix(Up2_boot_H1_1sided_times_dbarNt), 1, sum)
  int_dbarNtEL_SOcrit = as.vector(quantile(int_dbarNt_boot_H1, 1 - alpha))
  
  Up2_boot_H1_1sided_times_dt = as.matrix(Up2_boot_H1_1sided[, -(upperbindx_boot - lowerbindx_boot + 1)]) * as.matrix(wt_dt_boot[, lowerbindx_boot:(upperbindx_boot - 1)])
  int_dt_boot_H1  = apply(as.matrix(Up2_boot_H1_1sided_times_dt), 1, sum)
  int_dtEL_SOcrit = as.vector(quantile(int_dt_boot_H1, 1 - alpha))
  return(list(Up2_boot_H1_1sided = Up2_boot_H1_1sided, sup_boot_H1 = sup_boot_H1, EL_SOcrit = EL_SOcrit, wt_dbarNt = wt_dbarNt, wt_dbarNt_boot = wt_dbarNt_boot,
              wt_dt = wt_dt, wt_dt_boot = wt_dt_boot, wt_dF = wt_dF, wt_dF_boot = wt_dF_boot, int_dF_boot_H1 = int_dF_boot_H1, int_dFEL_SOcrit = int_dFEL_SOcrit,
              int_dbarNt_boot_H1 = int_dbarNt_boot_H1, int_dbarNtEL_SOcrit = int_dbarNtEL_SOcrit, int_dt_boot_H1 = int_dt_boot_H1, int_dtEL_SOcrit = int_dtEL_SOcrit))
}
