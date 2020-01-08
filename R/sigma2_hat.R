
sigma2_hat <- function(t, fit, fit1, fit2) {
    n   <- fit$n
    out <- 1:length(t)*0
    for (i in 1:length(t)) {
        d1     <- fit1$n.event[fit1$time <= t[i] & fit1$n.event != 0]
        r1     <- fit1$n.risk [fit1$time <= t[i] & fit1$n.event != 0]
        d2     <- fit2$n.event[fit2$time <= t[i] & fit2$n.event != 0]
        r2     <- fit2$n.risk [fit2$time <= t[i] & fit2$n.event != 0]
        out[i] <- n*(sum(division00(d1, r1*(r1 - d1))) + sum(division00(d2, r2*(r2 - d2))))
    }
    return (out)
}
