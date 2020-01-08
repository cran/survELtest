#' @importFrom stats uniroot
lambda0 <- function(t, fit1, fit2, tilde_theta = 1) {
    out <- 1:length(t)*0
    for (i in 1:length(t)) {
        out[i] <- if (sum(fit1$time[fit1$n.event != 0] <= t[i]) == 0 |
        sum(fit2$time[fit2$n.event != 0] <= t[i]) == 0) {
            NA
        } else {
            D1 <- max((fit1$n.event - fit1$n.risk)[fit1$time <= t[i] & fit1$n.event != 0])
            D2 <- max((fit2$n.event - fit2$n.risk)[fit2$time <= t[i] & fit2$n.event != 0])
            if (D1 != (-D2)) {
                uniroot(a_1, interval = c(D1 + 0.0001, -D2 - 0.0001), tol = 0.0001,
                fit1 = fit1, fit2 = fit2, t = t[i], tilde_theta = tilde_theta)$root
            } else D1
        }
    }
    return (out)
}
