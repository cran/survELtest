#' @importFrom stats median
#' @importFrom methods is
teststat_neg2logR = function(sided, k, M_vec, iter, iter_without_last, fitls, t, EL_CBcrit = 0) {
  # function for computing pointwise EL statistics
  ### functions related to the nonparametric likelihood ratio statistics:
  d = lapply(fitls, FUN = function (f) {
    f$n.event[f$time <= t & f$n.event != 0]
  })
  r = lapply(fitls, FUN = function (f) {
    f$n.risk[f$time <= t & f$n.event != 0]
  })
  A = lapply(iter, FUN = function (iter) {
    r[[iter]]-d[[iter]]
  })
  D = c()
  for (i in iter){
    D[i] = max(d[[i]] - r[[i]]) / M_vec[i] + 0.0001
  }
  cumsum_neg_D = -cumsum(D[k:2])[(k - 1):1]
  init_lambdas = rep(0, times = k)
  for (i in iter_without_last){
    init_lambdas[i + 1] = median(c(D[i] + init_lambdas[i], cumsum_neg_D[i]))
  }
  init_lambdas = init_lambdas[2:k]
  num_neg_loglik   = neg_log_likelihood_le_t_eq_h_k(t, fitls, M_vec, init_lambdas = init_lambdas)
  num_warn         = tryCatch(neg_log_likelihood_le_t_eq_h_k(t, fitls, M_vec, init_lambdas = init_lambdas),
                              error = function(e) e, warning = function(w) w)
  if (sided == 1) {
    denom_neg_loglik = neg_log_likelihood_le_t_ineq_h_k(t, fitls, M_vec,init_lambdas = init_lambdas)
    denom_warn       = tryCatch(neg_log_likelihood_le_t_ineq_h_k(t, fitls, M_vec, init_lambdas = init_lambdas), error = function(e) e, warning = function(w) w)
  } else {
    denom_neg_loglik = neg_log_likelihood_le_t_glob_h_k(t, fitls, M_vec,init_lambdas = init_lambdas)
    denom_warn       = tryCatch(neg_log_likelihood_le_t_glob_h_k(t, fitls, M_vec, init_lambdas = init_lambdas), error = function(e) e, warning = function(w) w)
  }
  init_dir1 = 1
  bdry_r    = cumsum_neg_D[1] - 0.0001
  bdry_l    = D[1] + 0.0001
  init_wid1 = 0.1 * (bdry_r - bdry_l)
  while (sum(abs(num_neg_loglik$resultEEs) >= 0.0001)        != 0 |
         sum(denom_neg_loglik$resultEEs    >= 10^-8)         != 0 |
         is(num_warn, "warning") * is(denom_warn, "warning") != 0)
  {
    init_lambdas[1] = init_lambdas[1] + init_wid1 * init_dir1 * (-1)^init_dir1
    if (init_lambdas[1] > bdry_r | init_lambdas[1] < bdry_l) break
    for (i in 2:(k - 1)) {
      init_lambdas[i] = median(c(D[i] + init_lambdas[i - 1], cumsum_neg_D[i]))
    }
    num_neg_loglik   = neg_log_likelihood_le_t_eq_h_k  (t, fitls, M_vec, init_lambdas = init_lambdas)
    num_warn         = tryCatch(neg_log_likelihood_le_t_eq_h_k(t, fitls, M_vec, init_lambdas = init_lambdas),
                                error = function(e) e, warning = function(w) w)
    if (sided == 1) {
      denom_neg_loglik = neg_log_likelihood_le_t_ineq_h_k(t, fitls, M_vec,init_lambdas = init_lambdas)
      denom_warn       = tryCatch(neg_log_likelihood_le_t_ineq_h_k(t, fitls, M_vec, init_lambdas = init_lambdas), error = function(e) e, warning = function(w) w)
    } else {
      denom_neg_loglik = neg_log_likelihood_le_t_glob_h_k(t, fitls, M_vec,init_lambdas = init_lambdas)
      denom_warn       = tryCatch(neg_log_likelihood_le_t_glob_h_k(t, fitls, M_vec, init_lambdas = init_lambdas), error = function(e) e, warning = function(w) w)
    }
    init_dir1 = init_dir1 + 1
  }
  
  still_error = (sum(abs(num_neg_loglik$resultEEs) >= 0.0001)        != 0 |
                   sum(denom_neg_loglik$resultEEs    >= 10^-8)         != 0 |
                   is(num_warn, "warning") * is(denom_warn, "warning") != 0)
  test1t      = 2 * (num_neg_loglik$objective - denom_neg_loglik$objective)
  return(list(out           = test1t - EL_CBcrit,
              still_error   = still_error,
              num_lambdas   = num_neg_loglik$lambdas,
              denom_lambdas = denom_neg_loglik$lambdas))
}
