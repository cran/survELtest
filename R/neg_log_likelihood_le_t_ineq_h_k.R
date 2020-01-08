#' @importFrom nloptr nloptr
neg_log_likelihood_le_t_ineq_h_k = function(t, fitls, M_vec,init_lambdas = rep(0, times = length(fitls)-1), maxeval = 10000) {
    k = length(fitls)
    iter = 1:k
    iter_without_last = 1:(k-1)
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
    init_weights = c()
    for (i in 1:k){
        D[i] = max(d[[i]] - r[[i]]) / M_vec[i] + 0.0001
        if (i == 1 ){
            init_weights[i] = list(d[[i]] / (r[[i]] + (M_vec[i] * init_lambdas[i])))
        }else if (i == k){
            init_weights[i] = list(d[[i]] / (r[[i]] - (M_vec[i] * init_lambdas[i-1])))
        }else{
            init_weights[i] = list(d[[i]] / (r[[i]] + (M_vec[i] * (init_lambdas[i] - init_lambdas[i-1]))))
        }
    }
    init_weights = as.numeric(unlist(init_weights))
    lgth_d = unlist(lapply(iter, FUN =function(iter){
        length(d[[iter]])
    }))
    local_opts = list("algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-10)
    opts_used = list("algorithm" = "NLOPT_LD_AUGLAG", "xtol_rel" = 1.0e-10, "maxeval" = maxeval, "local_opts" = local_opts)
    nlopt = nloptr(x0 = init_weights, eval_f = fnc_objective_h, eval_grad_f = gnc_objective_h,
    lb = rep(0.0, times = cumsum(lgth_d[1:k])[k]),
    ub = rep(1.0, times = cumsum(lgth_d[1:k])[k]),
    eval_g_ineq = constraints_h, eval_jac_g_ineq = jac_constraints_h, opts = opts_used,
    k = k, iter = iter, iter_without_last = iter_without_last, d = d, r = r, A = A, D = D, M_vec = M_vec)
    nlopt$Ds = D
    nlopt$resultEEs = constraints_h(nlopt$solution, k, iter, iter_without_last, d, r, A, D, M_vec)
    lambda = c()
    lambda[1] = (d[[1]][1] / nlopt$solution[1] - r[[1]][1]) / M_vec[1]
    if (k > 2){
        for (j in 2:(k-1)){
            lambda[j] = lambda[j-1] - ( r[[j]][1] - d[[j]][1] / nlopt$solution[1 + cumsum(lgth_d[1:k])[j] - lgth_d[j]]) / M_vec[j]
        }
    }
    nlopt$lambdas = lambda
    return (nlopt)
}
