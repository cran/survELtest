
a_1 <- function(lambda, fit1, fit2, t, tilde_theta = 1) {
  h1 <- if (sum(fit1$n.risk == fit1$n.event) == 0) {
    division00(fit1$n.event, (fit1$n.risk + lambda))
  } else {
    division00(fit1$n.event, (fit1$n.risk  + lambda))*as.numeric(fit1$n.risk != fit1$n.event) + 
    division00(fit1$n.event, (fit1$n.event + lambda))*as.numeric(fit1$n.risk == fit1$n.event)
  }
  h2 <- if (sum(fit2$n.risk == fit2$n.event) == 0) {
    division00(fit2$n.event, (fit2$n.risk - lambda))
  } else {
    division00(fit2$n.event, (fit2$n.risk  - lambda))*as.numeric(fit2$n.risk != fit2$n.event) + 
    division00(fit2$n.event, (fit2$n.event - lambda))*as.numeric(fit2$n.risk == fit2$n.event)
  }
  num   <- (1 - h1)[fit1$time <= t]
  denom <- (1 - h2)[fit2$time <= t]
  return (division00(prod(num), prod(denom)) - tilde_theta)
}
