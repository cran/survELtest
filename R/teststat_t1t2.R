#' @importFrom stats terms
#' @importFrom survival Surv
#' @importFrom survival survfit
teststat_t1t2 <- function(formula, group_order, t1, t2, nboot){
  # function for preparing data and restrict to [t1, t2]
  group = eval(parse(text = labels(terms(formula))))
  k = length(unique(group))
  M_vec = rep(1, k)
  iter = 1 : k
  iter_without_last = 1:(k-1)
  fit  <- survfit(update(formula,.~1))
  fit_group = survfit(formula)
  fitls = lapply(iter, FUN = function(iter) {
    fit_group[group = group_order[[iter]]]
  })
  NBOOT = nboot
  nn = c()
  for (i in iter){
    nn = c(nn,fitls[[i]]$n)
  }
  
  fit_time_restrict_boot = (fit$n.event != 0 & fit$time <= t2 & fit$time >= t1)  # ordered observed uncensored times
  
  mm = length(fit$time[fit_time_restrict_boot])
  S_hat_Wei = lapply(iter , FUN = function(iter){
    ((c(1, fitls[[iter]]$surv)[cumsum(c(0, fit$time) %in% c(0, fitls[[iter]]$time))])[-1])[fit_time_restrict_boot]
  })
  
  ### finding t \in [t_1, t_2] that are ordered observed uncensored times
  T_i1 = unlist(lapply(iter, FUN= function(iter){
    min((fit$time[fit_time_restrict_boot])[S_hat_Wei[[iter]] < 1])  # data driven t_1, as in Remark 3 of Theorem 1
  }))
  T_ms = unlist(lapply(iter, FUN = function(iter){
    max((fit$time[fit_time_restrict_boot])[S_hat_Wei[[iter]] > 0])  # data driven t_2, as in Remark 3 of Theorem 1
  } ))
  
  lowerb = max(T_i1)
  upperb = min(T_ms)
  lowerbindx_boot = which.min(abs(fit$time[fit_time_restrict_boot] - lowerb))
  upperbindx_boot = which.min(abs(fit$time[fit_time_restrict_boot] - upperb))
  
  Td_sort_boot = (fit$time[fit_time_restrict_boot])[lowerbindx_boot:upperbindx_boot]  # t \in [t_1, t_2] that are ordered observed uncensored times
  return(list(k = k, M_vec = M_vec, iter = iter, iter_without_last = iter_without_last, fit = fit, fitls = fitls, NBOOT = NBOOT,
              mm = mm, nn = nn, fit_time_restrict_boot = fit_time_restrict_boot, lowerbindx_boot = lowerbindx_boot,
              upperbindx_boot = upperbindx_boot, Td_sort_boot = Td_sort_boot))
}
