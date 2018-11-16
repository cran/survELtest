#' @importFrom nloptr nloptr

neg_log_likelihood_le_t_eq_h_k = function(t, fitls, M_vec,init_lambdas = rep(0, times = length(fitls)-1), maxeval = 10000) {
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
    fnc_objective_h = function(h) {
        lgth_d = unlist(lapply(iter, FUN =function(iter){
            length(d[[iter]])
        }))
        hh = lapply(iter, FUN =function(iter){
            h[ (1 + cumsum(lgth_d[1:k])[iter] - lgth_d[iter])  : cumsum(lgth_d[1:k])[iter] ]
        })
        sum_elmts = unlist(lapply(iter, FUN =function(iter){
            c(sum(d[[iter]]*log(hh[[iter]])),sum(product(A[[iter]], log(1-hh[[iter]]))))
        }))
        out = sum(sum_elmts)
        return (-out)
    }
    
    gnc_objective_h = function(h) {
        lgth_d = unlist(lapply(iter, FUN =function(iter){
            length(d[[iter]])
        }))
        hh = lapply(iter, FUN =function(iter){
            h[ (1 + cumsum(lgth_d[1:k])[iter] - lgth_d[iter])  : cumsum(lgth_d[1:k])[iter] ]
        })
        out = unlist(lapply(iter, FUN =function(iter){
            division00(d[[iter]], hh[[iter]]) - division00(r[[iter]] - d[[iter]], 1 - hh[[iter]])
        }))
        return(-out)
    }
    
    constraints_h = function(h) {
        lgth_d = unlist(lapply(iter, FUN =function(iter){
            length(d[[iter]])
        }))
        hh = lapply(iter, FUN =function(iter){
            h[ (1 + cumsum(lgth_d[1:k])[iter] - lgth_d[iter])  : cumsum(lgth_d[1:k])[iter] ]
        })
        out = unlist(lapply(iter_without_last, FUN =function(iter_without_last){
            log(division00(prod((1 - hh[[iter_without_last+1]]) ^ M_vec[iter_without_last+1]), prod((1 - hh[[iter_without_last]]) ^ M_vec[iter_without_last])))
        }))
        return(out)
    }
    
    jac_constraints_h = function(h) {
        lgth_d = unlist(lapply(iter, FUN =function(iter){
            length(d[[iter]])
        }))
        hh = lapply(iter, FUN =function(iter){
            h[ (1 + cumsum(lgth_d[1:k])[iter] - lgth_d[iter])  : cumsum(lgth_d[1:k])[iter] ]
        })
        out2D = lapply(iter_without_last,FUN =function(iter_without_last){
            out = c()
            if (iter_without_last != 1){
                out = c(out, rep(0, length = cumsum(lgth_d[1:k])[iter_without_last-1]))
            }
            out = c(out,division00(M_vec[iter_without_last], 1 - hh[[iter_without_last]]), division00(-M_vec[iter_without_last+1], 1 - hh[[iter_without_last+1]]))
            if (iter_without_last != k-1){
                out = c(out,rep(0, length = cumsum(lgth_d[1:k])[k] - cumsum(lgth_d[1:k])[iter_without_last+1]))
            }
            return(out)
        })
        rslt = c()
        for (i in iter_without_last){
            rslt = rbind(rslt, out2D[[i]])
        }
        return(rslt)
    }
    
    lgth_d = unlist(lapply(iter, FUN =function(iter){
        length(d[[iter]])
    }))
    local_opts = list("algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-10)
    opts_used = list("algorithm" = "NLOPT_LD_AUGLAG", "xtol_rel" = 1.0e-10, "maxeval" = maxeval, "local_opts" = local_opts)
    nlopt = nloptr(x0 = init_weights, eval_f = fnc_objective_h, eval_grad_f = gnc_objective_h,
    lb = rep(0.0, times = cumsum(lgth_d[1:k])[k]),
    ub = rep(1.0, times = cumsum(lgth_d[1:k])[k]),
    eval_g_eq = constraints_h, eval_jac_g_eq = jac_constraints_h, opts = opts_used)
    nlopt$Ds = D
    nlopt$resultEEs = constraints_h(nlopt$solution)
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
