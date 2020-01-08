#' @importFrom stats quantile
neg2ELratio_bootstrap_results <- function(sided, alpha, fit, fit1, fit2, NBOOT, fit_time_restrict_boot, mm, nn, lowerbindx_boot,  upperbindx_boot, U_boot_H1){
  # function for computing bootstrap results for supELtest and intELtest
  nboot      <- NBOOT
  EL_CBdistr <- apply((U_boot_H1^2),1,max)
  EL_CBcrit  <- as.vector(quantile(EL_CBdistr,1-alpha))
  
  if (sided == 1) {
    Up2_boot_H1_1sided <- (U_boot_H1 > 0) * (U_boot_H1^2)
  } else {
    Up2_boot_H1_1sided <- (U_boot_H1^2)
  }
  sup_boot_H1 <- apply(Up2_boot_H1_1sided[, lowerbindx_boot:upperbindx_boot], 1, max)
  EL_SOcrit   <- as.vector(quantile(sup_boot_H1, 1-alpha))
  
  nobd         <- length(rep(fit$time, times = fit$n.event))
  data_big     <- matrix(rep(rep(fit$time, times = fit$n.event), times = mm), byrow = TRUE, nrow = mm)
  fit_time_big <- matrix(rep(fit$time[fit_time_restrict_boot], times=nobd), byrow = FALSE, ncol = nobd)
  Indt_big     <- (data_big <= fit_time_big)
  barNt        <- apply(Indt_big,1, sum)/sum(nn)
  wt_dbarNt    <- diff(c(0, barNt))
  
  barNt_big      <- matrix(rep(barNt, time = nboot), byrow = TRUE, nrow = nboot)
  wt_dbarNt_boot <- t(apply(cbind(0, barNt_big), 1, diff))
  
  wt_dt      <- diff(c(fit$time[fit_time_restrict_boot], Inf))
  t_big      <- matrix(rep(fit$time[fit_time_restrict_boot], time = nboot), byrow = TRUE, nrow = nboot)
  wt_dt_boot <- t(apply(cbind(t_big, Inf), 1, diff))
  wt_dF      <- diff(c(0, 1 - fit$surv[fit_time_restrict_boot]))
  F_big      <- matrix(rep(1 - fit$surv[fit_time_restrict_boot], time = nboot), byrow = TRUE, nrow = nboot)
  wt_dF_boot <- t(apply(cbind(0, F_big), 1, diff))
  
  Up2_boot_H1_1sided_times_dF <- Up2_boot_H1_1sided*t(apply(cbind(0, F_big), 1, diff))
  int_dF_boot_H1  <- apply(as.matrix(Up2_boot_H1_1sided_times_dF[, lowerbindx_boot:upperbindx_boot]), 1, sum)
  int_dFEL_SOcrit <- as.vector(quantile(int_dF_boot_H1, 1-alpha))
  Up2_boot_H1_1sided_times_dbarNt <- Up2_boot_H1_1sided*t(apply(cbind(0, barNt_big), 1, diff))
  int_dbarNt_boot_H1  <- apply(as.matrix(Up2_boot_H1_1sided_times_dbarNt[, lowerbindx_boot:upperbindx_boot]), 1, sum)
  int_dbarNtEL_SOcrit <- as.vector(quantile(int_dbarNt_boot_H1, 1-alpha))
  
  Up2_boot_H1_1sided_times_dt <- Up2_boot_H1_1sided*t(apply(cbind(t_big, Inf), 1, diff))
  int_dt_boot_H1  <- apply(as.matrix(Up2_boot_H1_1sided_times_dt[, lowerbindx_boot:(upperbindx_boot - 1)]), 1, sum)
  int_dtEL_SOcrit <- as.vector(quantile(int_dt_boot_H1, 1-alpha))
  return(list(EL_CBcrit = EL_CBcrit, Up2_boot_H1_1sided = Up2_boot_H1_1sided, sup_boot_H1 = sup_boot_H1, EL_SOcrit = EL_SOcrit, barNt = barNt,
              wt_dbarNt = wt_dbarNt, wt_dbarNt_boot = wt_dbarNt_boot, wt_dt = wt_dt, wt_dt_boot = wt_dt_boot,
              wt_dF = wt_dF, wt_dF_boot = wt_dF_boot, int_dF_boot_H1 = int_dF_boot_H1, int_dFEL_SOcrit = int_dFEL_SOcrit, int_dbarNt_boot_H1 = int_dbarNt_boot_H1,
              int_dbarNtEL_SOcrit = int_dbarNtEL_SOcrit, int_dt_boot_H1 = int_dt_boot_H1, int_dtEL_SOcrit = int_dtEL_SOcrit))
}
