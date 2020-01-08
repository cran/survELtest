#' @importFrom stats quantile
teststat_nocrossings <- function(alpha, k, M_vec, iter, iter_without_last, fit, fitls, nn, fit_time_restrict_boot, lowerbindx_boot, upperbindx_boot, sum_DWknGImuw_big){
  # function for computation for nocrossings
  ### computing quantities related to the first step of our composite procedure:
  test_nocross = 0
  # pmax_elmts is a two-dimension list used in next parameter HW_CBdistr with apply(abs(M_vec[i] * sum_DWknGImuw_{i}_big - M_vec[i]
  pmax_elmts = lapply(iter_without_last, FUN=function(iter_without_last){
    apply(abs(M_vec[iter_without_last+1] * as.matrix(sum_DWknGImuw_big[[iter_without_last+1]][, lowerbindx_boot:upperbindx_boot]) - M_vec[iter_without_last] * as.matrix(sum_DWknGImuw_big[[iter_without_last]][, lowerbindx_boot:upperbindx_boot])), 1, max)
  })
  HW_CBdistr = pmax_elmts[[1]]
  if (k >2){
    for (i in 2:(k-1)){
      HW_CBdistr = pmax(HW_CBdistr, pmax_elmts[[i]])
    }
  }
  HW_CBcrit = as.vector(quantile(HW_CBdistr, 1-alpha))
  
  ### Kaplan--Meier estimates for each treatment group and related quantities for the first step of our composite procedure:
  S_hat = lapply(iter, FUN= function(iter){
    (((c(1, fitls[[iter]]$surv)[cumsum(c(0, fit$time) %in% c(0, fitls[[iter]]$time))])[-1])[fit_time_restrict_boot])[lowerbindx_boot:upperbindx_boot]
  })
  diff = lapply(iter_without_last , FUN = function(iter_without_last){
    M_vec[iter_without_last] * log(S_hat[[iter_without_last]]) - M_vec[iter_without_last+1] * log(S_hat[[iter_without_last+1]])
  })
  
  HW_CB_ubs_pairwise = lapply(iter_without_last, FUN = function(iter_without_last){
    diff[[iter_without_last]] + HW_CBcrit / sqrt(sum(nn))
  })
  HW_CB_lbs_pairwise = lapply(iter_without_last, FUN =function(iter_without_last){
    diff[[iter_without_last]] - HW_CBcrit / sqrt(sum(nn))
  })
  
  g_S_all = sum(HW_CB_ubs_pairwise[[1]] >= 0) == (upperbindx_boot - lowerbindx_boot + 1) & sum(HW_CB_lbs_pairwise[[1]] <= 0) != (upperbindx_boot - lowerbindx_boot + 1)
  if (k >2){
    for (i in 2:(k-1)){
      g_S_all = g_S_all & sum(HW_CB_ubs_pairwise[[i]] >= 0) == (upperbindx_boot - lowerbindx_boot + 1) & sum(HW_CB_lbs_pairwise[[i]] <= 0) != (upperbindx_boot - lowerbindx_boot + 1)
    }
  }
  eq_S_all = sum(HW_CB_lbs_pairwise[[1]] <= 0) == (upperbindx_boot - lowerbindx_boot + 1) & sum(HW_CB_ubs_pairwise[[1]] >= 0) == (upperbindx_boot - lowerbindx_boot + 1)
  if (k > 2){
    for (i in 2: (k-1)){
      eq_S_all = eq_S_all & sum(HW_CB_lbs_pairwise[[i]] <= 0) == (upperbindx_boot - lowerbindx_boot + 1) & sum(HW_CB_ubs_pairwise[[i]] >= 0) == (upperbindx_boot - lowerbindx_boot + 1)
    }
  }
  if (g_S_all | eq_S_all) {
    test_nocross = 1
  }
  
  return(list(test_nocross = test_nocross))
}

