#' @importFrom stats update
#' @importFrom stats terms
#' @importFrom survival Surv
#' @importFrom survival survfit
neg2ELratio_t1t2 <- function(formula, group_order, t1, t2, nboot) {
  # function for preparing data and restrict to [t1, t2]
  fit  <- survfit(update(formula,.~1))
  fit_group = survfit(formula)
  group = eval(parse(text = labels(terms(formula))))
  fit1 = fit_group[group = group_order[[1]]]
  fit2 = fit_group[group = group_order[[2]]]
  NBOOT <- nboot
  
  fit_time_restrict_boot <- (fit$n.event != 0 & fit$time <= t2 & fit$time >= t1)
  
  mm <- length(fit$time[fit_time_restrict_boot])
  S1_hat_Wei <- ((c(1, fit1$surv)[cumsum(c(0, fit$time) %in% c(0, fit1$time))])[-1])[fit_time_restrict_boot]
  S2_hat_Wei <- ((c(1, fit2$surv)[cumsum(c(0, fit$time) %in% c(0, fit2$time))])[-1])[fit_time_restrict_boot]
  
  nn   <- c(fit1$n, fit2$n)
  T_11 <- min((fit$time[fit_time_restrict_boot])[S1_hat_Wei < 1])
  T_21 <- min((fit$time[fit_time_restrict_boot])[S2_hat_Wei < 1])
  T_1m <- max((fit$time[fit_time_restrict_boot])[S1_hat_Wei > 0])
  T_2m <- max((fit$time[fit_time_restrict_boot])[S2_hat_Wei > 0])
  
  lowerb <- max(T_11, T_21)
  upperb <- min(T_1m, T_2m)
  
  lowerbindx_boot <- which.min(abs(fit$time[fit_time_restrict_boot] - lowerb))
  upperbindx_boot <- which.min(abs(fit$time[fit_time_restrict_boot] - upperb))
  
  if (lowerbindx_boot >= upperbindx_boot) {
    warning("lowerbindx_boot >= upperbindx_boot in computing the integral and sup statistic.\nEither your sample size is too small or the overlapping region of the two samples is empty or just one observed event time.")
    return (NULL)
  }
  Td_sort_boot <- (fit$time[fit_time_restrict_boot])[lowerbindx_boot:upperbindx_boot]
  return(list(fit = fit, fit1 = fit1, fit2 = fit2, NBOOT = NBOOT, fit_time_restrict_boot = fit_time_restrict_boot,
              mm = mm, nn = nn, T_11 = T_11, T_21 = T_21, lowerbindx_boot = lowerbindx_boot,
              upperbindx_boot = upperbindx_boot, Td_sort_boot = Td_sort_boot))
}
