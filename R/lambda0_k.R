#' @importFrom stats uniroot

lambda0_k = function(t, fit1, fit2, tilde_theta = 1, M_vec2) {
    out = 1:length(t) * 0
    for (i in 1:length(t)) {
        if (sum(fit1$time[fit1$n.event!=0] <= t[i]) == 0 | sum(fit2$time[fit2$n.event != 0] <= t[i]) == 0) {
            out[i] = NA
        }  else {
            D1 = max((fit1$n.event - fit1$n.risk)[fit1$time <= t[i] & fit1$n.event != 0]) / M_vec2[1]
            D2 = max((fit2$n.event - fit2$n.risk)[fit2$time <= t[i] & fit2$n.event != 0]) / M_vec2[2]
            if (D1 != (-D2)) {
                out[i] = uniroot(a_1_k, interval = c(D1 + 0.0001, -D2 - 0.0001), tol = 0.0001, fit1 = fit1, fit2 = fit2, t = t[i], tilde_theta = tilde_theta, M_vec2 = M_vec2)$root
            } else {
                out[i] = D1
            }
        }
    }
    return(out)
}
