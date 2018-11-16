
a_1_k = function(lambda, fit1, fit2, t, tilde_theta = 1, M_vec2) {
    if (sum(fit1$n.risk == fit1$n.event) == 0) {
        h1 = division00(fit1$n.event, (fit1$n.risk + M_vec2[1] * lambda))
    } else {
        h1 = division00(fit1$n.event, (fit1$n.risk + M_vec2[1] * lambda)) * as.numeric(fit1$n.risk != fit1$n.event) + division00(fit1$n.event, (fit1$n.event + M_vec2[1] * lambda)) * as.numeric(fit1$n.risk == fit1$n.event)
    }
    if (sum(fit2$n.risk == fit2$n.event) == 0) {
        h2 = division00(fit2$n.event, (fit2$n.risk - M_vec2[2] * lambda))
    } else {
        h2 = division00(fit2$n.event, (fit2$n.risk - M_vec2[2] * lambda)) * as.numeric(fit2$n.risk != fit2$n.event) + division00(fit2$n.event, (fit2$n.event - M_vec2[2] * lambda)) * as.numeric(fit2$n.risk == fit2$n.event)
    }
    
    num = ((1 - h1) ^ M_vec2[1])[fit1$time <= t]
    denom = ((1 - h2) ^ M_vec2[2])[fit2$time <= t]
    
    return(division00(prod(num), prod(denom)) - tilde_theta)
}  
