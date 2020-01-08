#' @importFrom stats rnorm
neg2ELratio_bootstrap <- function(nboot, nlimit, fit, fit1, fit2, NBOOT, fit_time_restrict_boot, mm, nn){
  # function for bootstrap
  nobd1 <- length(rep(fit1$time, times = fit1$n.event))
  nobd2 <- length(rep(fit2$time, times = fit2$n.event))
  
  data1_big      <- matrix(rep(rep(fit1$time, times = fit1$n.event), times = mm), byrow = TRUE, nrow = mm)
  fit_time_1_big <- matrix(rep(fit$time[fit_time_restrict_boot], times = nobd1), byrow = FALSE, ncol = nobd1)
  Ind1t_big      <- (data1_big <= fit_time_1_big)
  Ind1x          <- rep(fit1$n.risk, times = fit1$n.event)
  Ind1x_big      <- matrix(rep(Ind1x, times = mm), byrow = TRUE, nrow = mm)
  
  data2_big      <- matrix(rep(rep(fit2$time, times = fit2$n.event), times = mm), byrow = TRUE, nrow = mm)
  fit_time_2_big <- matrix(rep(fit$time[fit_time_restrict_boot], times = nobd2), byrow = FALSE, ncol = nobd2)
  Ind2t_big      <- (data2_big <= fit_time_2_big)
  Ind2x          <- rep(fit2$n.risk, times = fit2$n.event)
  Ind2x_big      <- matrix(rep(Ind2x, times = mm), byrow = TRUE, nrow = mm)
  
  nsplit <- ceiling(mm/nlimit)
  if (nsplit > 1) {
    U_boot_H1 <- matrix(rep(NA, NBOOT*mm), NBOOT, mm)
    for (i in 1:nsplit) {
      if (i == nsplit) {
        nboot <- NBOOT - floor(NBOOT/nsplit)*(nsplit - 1)
        Gs_1  <- matrix(rnorm(nobd1*nboot), nrow = nboot, ncol = nobd1)
        Imuw_BIG_1 <- array(rep(Ind1t_big/Ind1x_big, each = nboot), c(nboot, mm, nobd1))
        Gs_1_BIG   <- array(matrix(rep(t(Gs_1), each = mm), byrow = TRUE, nrow = nboot), c(nboot, mm, nobd1))
        sum_DWknGImuw_1_big <- sqrt(sum(nn))*apply(Imuw_BIG_1*Gs_1_BIG, c(1, 2), sum)
        Gs_2       <- matrix(rnorm(nobd2*nboot), nrow = nboot, ncol = nobd2)
        Imuw_BIG_2 <- array(rep(Ind2t_big/Ind2x_big, each = nboot), c(nboot, mm, nobd2))
        Gs_2_BIG   <- array(matrix(rep(t(Gs_2), each = mm), byrow = TRUE, nrow = nboot), c(nboot, mm, nobd2))
        sum_DWknGImuw_2_big <- sqrt(sum(nn))*apply(Imuw_BIG_2*Gs_2_BIG, c(1, 2), sum)
        DsqSig_big <- matrix(rep(1/sqrt(sigma2_hat(fit$time[fit_time_restrict_boot], fit, fit1, fit2)), times = nboot), byrow = TRUE, nrow = nboot)
        U_boot_H1[(NBOOT - nboot + 1):NBOOT, ] <- product_mat(- sum_DWknGImuw_1_big + sum_DWknGImuw_2_big, DsqSig_big)
      } else {
        nboot <- floor(NBOOT/nsplit)
        Gs_1  <- matrix(rnorm(nobd1*nboot), nrow = nboot, ncol = nobd1)
        Imuw_BIG_1 <- array(rep(Ind1t_big/Ind1x_big, each = nboot), c(nboot, mm, nobd1))
        Gs_1_BIG   <- array(matrix(rep(t(Gs_1), each = mm), byrow = TRUE, nrow = nboot), c(nboot, mm, nobd1))
        sum_DWknGImuw_1_big <- sqrt(sum(nn))*apply(Imuw_BIG_1*Gs_1_BIG, c(1, 2), sum)
        Gs_2       <- matrix(rnorm(nobd2*nboot), nrow = nboot, ncol = nobd2)
        Imuw_BIG_2 <- array(rep(Ind2t_big/Ind2x_big, each = nboot), c(nboot, mm, nobd2))
        Gs_2_BIG   <- array(matrix(rep(t(Gs_2), each = mm), byrow = TRUE, nrow = nboot), c(nboot, mm, nobd2))
        sum_DWknGImuw_2_big <- sqrt(sum(nn))*apply(Imuw_BIG_2*Gs_2_BIG, c(1, 2), sum)
        DsqSig_big <- matrix(rep(1/sqrt(sigma2_hat(fit$time[fit_time_restrict_boot], fit, fit1, fit2)), times = nboot), byrow = TRUE, nrow = nboot)
        U_boot_H1[((i - 1)*nboot + 1):(i*nboot), ] <- product_mat(- sum_DWknGImuw_1_big + sum_DWknGImuw_2_big, DsqSig_big)
      }
    }
  } else {
    Gs_1       <- matrix(rnorm(nobd1*nboot), nrow = nboot, ncol = nobd1)
    Imuw_BIG_1 <- array(rep(Ind1t_big/Ind1x_big, each = nboot), c(nboot, mm, nobd1))
    Gs_1_BIG   <- array(matrix(rep(t(Gs_1), each = mm), byrow = TRUE, nrow = nboot), c(nboot, mm, nobd1))
    sum_DWknGImuw_1_big <- sqrt(sum(nn))*apply(Imuw_BIG_1*Gs_1_BIG, c(1, 2), sum)
    Gs_2       <- matrix(rnorm(nobd2*nboot), nrow = nboot, ncol = nobd2)
    Imuw_BIG_2 <- array(rep(Ind2t_big/Ind2x_big, each = nboot), c(nboot, mm, nobd2))
    Gs_2_BIG   <- array(matrix(rep(t(Gs_2), each = mm), byrow = TRUE, nrow = nboot), c(nboot, mm, nobd2))
    sum_DWknGImuw_2_big <- sqrt(sum(nn))*apply(Imuw_BIG_2*Gs_2_BIG, c(1, 2), sum)
    DsqSig_big <- matrix(rep(1/sqrt(sigma2_hat(fit$time[fit_time_restrict_boot], fit, fit1, fit2)), times = nboot), byrow = TRUE, nrow = nboot)
    U_boot_H1  <- product_mat(- sum_DWknGImuw_1_big + sum_DWknGImuw_2_big, DsqSig_big)
    }
  return(list(sum_DWknGImuw_1_big = sum_DWknGImuw_1_big, sum_DWknGImuw_2_big = sum_DWknGImuw_2_big, U_boot_H1 = U_boot_H1))
}
