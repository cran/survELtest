
neg2ELratio_EL <- function(sided, fit1, fit2, T_11, T_21, upperbindx_boot, lowerbindx_boot, Td_sort_boot, EL_CBcrit){
  # function for computing pointwise EL statistics
  teststat_pre <- 1:length(Td_sort_boot)*0
  for (j in 1:(upperbindx_boot - lowerbindx_boot + 1)) {
    t <- Td_sort_boot[j]
    lambda0_hat <- lambda0(t, fit1, fit2, tilde_theta = 1)
    if (lambda0_hat < 0 | sided != 1) {
      d1 <- fit1$n.event[fit1$time <= t & fit1$n.event != 0]
      r1 <- fit1$n.risk [fit1$time <= t & fit1$n.event != 0]
      d2 <- fit2$n.event[fit2$time <= t & fit2$n.event != 0]
      r2 <- fit2$n.risk [fit2$time <= t & fit2$n.event != 0]
      A1 <- r1 - d1
      A2 <- r2 - d2
      B1 <- d1/(r1 + lambda0_hat)*as.numeric(r1 != d1) + d1/(d1 + lambda0_hat)*as.numeric(r1 == d1)
      B2 <- d2/(r2 - lambda0_hat)*as.numeric(r2 != d2) + d2/(d2 - lambda0_hat)*as.numeric(r2 == d2)
      EL_CBcrit <- 0
      teststat_pre[j] <- -2*sum(d1*log(B1)) - 2*sum(product(A1, log(1 - B1))) -
        2*sum(d2*log(B2)) - 2*sum(product(A2, log(1 - B2))) +
        2*sum(d1*log(d1/r1)) + 2*sum(d2*log(d2/r2)) +
        2*sum(product(A1, log(1 - d1/r1))) + 2*sum(product(A2, log(1 - d2/r2))) -
        EL_CBcrit
    } else if (t < T_11 & t >= T_21) {
      teststat_pre[j]=0
    }
  }
  return(list(teststat_pre = teststat_pre))
}
