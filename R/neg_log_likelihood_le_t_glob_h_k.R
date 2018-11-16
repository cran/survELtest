
neg_log_likelihood_le_t_glob_h_k = function(t, fitls, M_vec,init_lambdas = rep(0, times = length(fitls)-1), maxeval = 10000) {
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
    nlopt = list()
    nlopt$objective = -sum(sapply(iter, FUN = function (iter) {
        sum(d[[iter]] * log(d[[iter]] / r[[iter]])) + sum(A[[iter]] * log(1 - d[[iter]] / r[[iter]]))
    }))
    nlopt$resultEEs = 1:(k-1) * 0
    nlopt$lambdas = 1:(k-1) * 0
    return (nlopt)
}

