#' @importFrom survival Surv
#' @importFrom survival survfit
#' @importFrom stats rnorm
#' @importFrom stats quantile

neg2ELratio <- function(data, group_order, t1, t2, sided, nboot, alpha, details.return, seed, nlimit) {
    set.seed(seed)
    
    dat  <- Surv(data[, 1], data[, 2])
    fit  <- survfit(dat ~ 1)
    fit1 <- survfit(dat[data[, 3] == group_order[1]] ~ 1)
    fit2 <- survfit(dat[data[, 3] == group_order[2]] ~ 1)
    
    # fit1  <- survfit(dat[data[, 3] == 1] ~ 1) #
    # fit2  <- survfit(dat[data[, 3] == 2] ~ 1) #
    NBOOT <- nboot #
    
    nn   <- c(fit1$n, fit2$n)
    T_11 <- min(fit1$time[fit1$n.event != 0])
    T_21 <- min(fit2$time[fit2$n.event != 0])
    T_1m <- max(fit1$time[fit1$n.event != 0])
    T_2m <- max(fit2$time[fit2$n.event != 0])
    fit_time_restrict_boot <- (fit$n.event != 0 & fit$time <= t2 & fit$time >= t1)
    
    # fit_time_restrict_boot <- (fit$n.event != 0) #
    
    mm <- length(fit$time[fit_time_restrict_boot])
    
    nsplit <- ceiling(mm/nlimit) #
    lowerb <- max(T_11, T_21)
    upperb <- min(T_1m, T_2m)
    
    lowerbindx_boot <- which.min(abs(fit$time[fit_time_restrict_boot] - lowerb))
    upperbindx_boot <- which.min(abs(fit$time[fit_time_restrict_boot] - upperb))
    
    if (lowerbindx_boot >= upperbindx_boot) {
        warning("lowerbindx_boot >= upperbindx_boot in computing the integral and sup statistic.
        Either your sample size is too small
        or the overlapping region of the two samples is empty or just one observed event time.")
        return (NULL)
    }
    
    # lowerbindx_boot <- (1:mm)[abs(fit$time[fit_time_restrict_boot] - lowerb) == min(abs(fit$time[fit_time_restrict_boot] - lowerb))] #
    # upperbindx_boot <- (1:mm)[abs(fit$time[fit_time_restrict_boot] - upperb) == min(abs(fit$time[fit_time_restrict_boot] - upperb))] #
    
    Td_sort_boot <- (fit$time[fit_time_restrict_boot])[lowerbindx_boot:upperbindx_boot]
    nobd1 <- length(rep(fit1$time, times = fit1$n.event))
    nobd2 <- length(rep(fit2$time, times = fit2$n.event))
    
    data1_big      <- matrix(rep(rep(fit1$time, times = fit1$n.event), times = mm), byrow = TRUE, nrow = mm)
    fit_time_1_big <- matrix(rep(fit$time[fit_time_restrict_boot], times = nobd1), byrow = FALSE, ncol = nobd1)
    Ind1t_big      <- (data1_big <= fit_time_1_big)
    Ind1x          <- rep(fit1$n.risk, times = fit1$n.event)
    Ind1x_big      <- matrix(rep(Ind1x, times = mm), byrow = TRUE, nrow = mm)
    
    # Gs_1           <- matrix(rnorm(nobd1*nboot), nrow = nboot, ncol = nobd1)
    # Imuw_BIG_1     <- array(rep(Ind1t_big/Ind1x_big, each = nboot), c(nboot, mm, nobd1))
    # Gs_1_BIG       <- array(matrix(rep(t(Gs_1), each = mm), byrow = TRUE, nrow = nboot), c(nboot, mm, nobd1))
    # sum_DWknGImuw_1_big <- sqrt(sum(nn))*apply(Imuw_BIG_1*Gs_1_BIG, c(1, 2), sum)
    
    data2_big      <- matrix(rep(rep(fit2$time, times = fit2$n.event), times = mm), byrow = TRUE, nrow = mm)
    fit_time_2_big <- matrix(rep(fit$time[fit_time_restrict_boot], times = nobd2), byrow = FALSE, ncol = nobd2)
    Ind2t_big      <- (data2_big <= fit_time_2_big)
    Ind2x          <- rep(fit2$n.risk, times = fit2$n.event)
    Ind2x_big      <- matrix(rep(Ind2x, times = mm), byrow = TRUE, nrow = mm)
    
    # Gs_2           <- matrix(rnorm(nobd2*nboot), nrow = nboot, ncol = nobd2)
    # Imuw_BIG_2     <- array(rep(Ind2t_big/Ind2x_big, each = nboot), c(nboot, mm, nobd2))
    # Gs_2_BIG       <- array(matrix(rep(t(Gs_2), each = mm), byrow = TRUE, nrow = nboot), c(nboot, mm, nobd2))
    # sum_DWknGImuw_2_big <- sqrt(sum(nn))*apply(Imuw_BIG_2*Gs_2_BIG, c(1, 2), sum)
    
    ########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
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
                U_boot_H1[(NBOOT - nboot + 1):NBOOT, ] <- product_mat(sum_DWknGImuw_1_big + sum_DWknGImuw_2_big, DsqSig_big)
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
                U_boot_H1[((i - 1)*nboot + 1):(i*nboot), ] <- product_mat(sum_DWknGImuw_1_big + sum_DWknGImuw_2_big, DsqSig_big)
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
        U_boot_H1  <- product_mat(sum_DWknGImuw_1_big + sum_DWknGImuw_2_big, DsqSig_big)
    }
    
    ########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
    S1_hat_Wei <- ((c(1, fit1$surv)[cumsum(c(0, fit$time) %in% c(0, fit1$time))])[-1])[fit_time_restrict_boot]
    S2_hat_Wei <- ((c(1, fit2$surv)[cumsum(c(0, fit$time) %in% c(0, fit2$time))])[-1])[fit_time_restrict_boot]
    
    S2_hat_big <- matrix(rep(S2_hat_Wei, times = nboot), byrow = TRUE, nrow = nboot)
    S1_hat_big <- matrix(rep(S1_hat_Wei, times = nboot), byrow = TRUE, nrow = nboot)
    DsqSig_S1mS2_big <- matrix(rep(1/sqrt(sigma2_hat_S1mS2(fit$time[fit_time_restrict_boot], fit, fit1, fit2)), times = nboot), byrow = TRUE, nrow = nboot)
    ourZt_boot       <- product_mat(sum_DWknGImuw_2_big/sqrt(sum(nn))*S2_hat_big - sum_DWknGImuw_1_big/sqrt(sum(nn))*S1_hat_big, DsqSig_S1mS2_big)
    
    nboot      <- NBOOT #
    EL_CBdistr <- apply((U_boot_H1^2),1,max) #
    EL_CBcrit  <- as.vector(quantile(EL_CBdistr,1-alpha)) #
    
    if (sided == 1) {
        Up2_boot_H1_1sided <- (U_boot_H1  > 0) * (U_boot_H1^2)
        Up2_Zt_boot        <- (ourZt_boot > 0) * (ourZt_boot^2)
    } else {
        Up2_boot_H1_1sided <- (U_boot_H1^2)
        Up2_Zt_boot        <- (ourZt_boot^2)
    }
    sup_boot_H1 <- apply(Up2_boot_H1_1sided[, lowerbindx_boot:upperbindx_boot], 1, max) #
    EL_SOcrit   <- as.vector(quantile(sup_boot_H1, 1-alpha)) #
    
    nobd         <- length(rep(fit$time, times = fit$n.event))
    data_big     <- matrix(rep(rep(fit$time, times = fit$n.event), times = mm), byrow = TRUE, nrow = mm)
    fit_time_big <- matrix(rep(fit$time[fit_time_restrict_boot], times=nobd), byrow = FALSE, ncol = nobd)
    Indt_big     <- (data_big <= fit_time_big)
    barNt        <- apply(Indt_big,1, sum)/sum(nn)
    wt_dbarNt    <- diff(c(0, barNt))
    
    barNt_big      <- matrix(rep(barNt, time = nboot), byrow = TRUE, nrow = nboot)
    wt_dbarNt_boot <- t(apply(cbind(0, barNt_big), 1, diff))
    
    wt_dt      <- diff(c(fit$time[fit_time_restrict_boot], Inf))
    t_big      <- matrix(rep(fit$time[fit_time_restrict_boot], time = nboot), byrow = TRUE, nrow = nboot)
    wt_dt_boot <- t(apply(cbind(t_big, Inf), 1, diff))
    Sig_big    <- matrix(rep(sigma2_hat(fit$time[fit_time_restrict_boot], fit, fit1, fit2), times = nboot), byrow = TRUE, nrow = nboot)
    wt_dF      <- diff(c(0, 1 - fit$surv[fit_time_restrict_boot]))
    F_big      <- matrix(rep(1 - fit$surv[fit_time_restrict_boot], time = nboot), byrow = TRUE, nrow = nboot)
    wt_dF_boot <- t(apply(cbind(0, F_big), 1, diff))
    
    Up2_boot_H1_1sided_times_dF <- Up2_boot_H1_1sided*t(apply(cbind(0, F_big), 1, diff)) #
    int_dF_boot_H1  <- apply(as.matrix(Up2_boot_H1_1sided_times_dF[, lowerbindx_boot:upperbindx_boot]), 1, sum) #
    int_dFEL_SOcrit <- as.vector(quantile(int_dF_boot_H1, 1-alpha)) #
    
    Up2_boot_H1_1sided_times_dbarNt <- Up2_boot_H1_1sided*t(apply(cbind(0, barNt_big), 1, diff)) #
    int_dbarNt_boot_H1  <- apply(as.matrix(Up2_boot_H1_1sided_times_dbarNt[, lowerbindx_boot:upperbindx_boot]), 1, sum) #
    int_dbarNtEL_SOcrit <- as.vector(quantile(int_dbarNt_boot_H1, 1-alpha)) #
    
    Up2_boot_H1_1sided_times_dt <- Up2_boot_H1_1sided*t(apply(cbind(t_big, Inf), 1, diff)) #
    int_dt_boot_H1  <- apply(as.matrix(Up2_boot_H1_1sided_times_dt[, lowerbindx_boot:(upperbindx_boot - 1)]), 1, sum) #
    int_dtEL_SOcrit <- as.vector(quantile(int_dt_boot_H1, 1-alpha)) #
    
    teststat_pre <- 1:length(Td_sort_boot)*0
    for (j in 1:(upperbindx_boot - lowerbindx_boot + 1)) {
        # t <- Td_sort_boot[j]
        # lambda0_hat <- lambda0(t, fit1, fit2, tilde_theta = 1)
        # d1 <- fit1$n.event[fit1$time <= t & fit1$n.event != 0]
        # r1 <- fit1$n.risk [fit1$time <= t & fit1$n.event != 0]
        # d2 <- fit2$n.event[fit2$time <= t & fit2$n.event != 0]
        # r2 <- fit2$n.risk [fit2$time <= t & fit2$n.event != 0]
        # A1 <- r1 - d1
        # A2 <- r2 - d2
        # B1 <- d1/(r1 + lambda0_hat)
        # B2 <- d2/(r2 - lambda0_hat)
        # teststat_pre[j] <- -2*sum(d1*log(B1)) - 2*sum(product(A1, log(1 - B1))) -
        #                     2*sum(d2*log(B2)) - 2*sum(product(A2, log(1 - B2))) +
        #                     2*sum(d1*log(d1/r1)) + 2*sum(d2*log(d2/r2)) +
        #                     2*sum(product(A1, log(1 - d1/r1))) + 2*sum(product(A2, log(1 - d2/r2)))
        #
        # if (sided==1 & lambda0_hat >= 0) teststat_pre[j] <- 0
        
        ########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
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
            EL_CBcrit       <- 0
            teststat_pre[j] <- -2*sum(d1*log(B1)) - 2*sum(product(A1, log(1 - B1))) -
            2*sum(d2*log(B2)) - 2*sum(product(A2, log(1 - B2))) +
            2*sum(d1*log(d1/r1)) + 2*sum(d2*log(d2/r2)) +
            2*sum(product(A1, log(1 - d1/r1))) + 2*sum(product(A2, log(1 - d2/r2))) -
            EL_CBcrit
        } else if (t < T_11 & t >= T_21) {
            teststat_pre[j]=0
        }
        # print(j)
        # print(Sys.time())
        ########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
    }
    
    D_hat_Wei <- S1_hat_Wei - S2_hat_Wei
    ourZt     <- product(D_hat_Wei, 1/sqrt(sigma2_hat_S1mS2(fit$time[fit_time_restrict_boot], fit, fit1, fit2)))
    teststat_pre_endpt <- if (sided == 1) (ourZt>0)*(ourZt^2) else ourZt^2
    
    # if (details.return == TRUE) {
    #   return (list(neg2ELratio_at_ts = teststat_pre, stat_compo = stat_compo, compo_vs_not = compo_vs_not,
    #                wt_dbarNt = wt_dbarNt, wt_dt = wt_dt, wt_db = wt_db, wt_dF = wt_dF,
    #                neg2ELratio_bootstrap_at_ts = Up2_boot_H1_1sided, Up2_boot_compo = Up2_boot_compo,
    #                wt_dbarNt_boot = wt_dbarNt_boot, wt_dt_boot = wt_dt_boot, wt_db_boot = wt_db_boot, wt_dF_boot = wt_dF_boot,
    #                lowerbindx_boot = lowerbindx_boot, upperbindx_boot = upperbindx_boot,
    #                ts = Td_sort_boot, bootstrap_ts = fit$time[fit_time_restrict_boot]))
    # } else return (list(neg2ELratio_at_ts = teststat_pre))
    test_nocross = 0
    HW_CBdistr=apply(abs(sum_DWknGImuw_2_big[,lowerbindx_boot:upperbindx_boot]-sum_DWknGImuw_1_big[,lowerbindx_boot:upperbindx_boot]),1,max)
    HW_CBcrit=as.vector(quantile(HW_CBdistr,1-alpha))
    S1_hat=(((c(1,fit1$surv)[cumsum(c(0,fit$time) %in% c(0,fit1$time))])[-1])[fit_time_restrict_boot])[lowerbindx_boot:upperbindx_boot] #except for survival, other "expansions" can't use this
    S2_hat=(((c(1,fit2$surv)[cumsum(c(0,fit$time) %in% c(0,fit2$time))])[-1])[fit_time_restrict_boot])[lowerbindx_boot:upperbindx_boot]
    diff=log(S1_hat)-log(S2_hat)
    HW_CB_ubs=diff+HW_CBcrit/sqrt(sum(nn))
    HW_CB_lbs=diff-HW_CBcrit/sqrt(sum(nn))
    eq_S_all = sum(HW_CB_lbs <= 0) == (upperbindx_boot - lowerbindx_boot + 1) & sum(HW_CB_ubs >= 0) == (upperbindx_boot - lowerbindx_boot + 1)
    g_S_all = sum(HW_CB_ubs >= 0) == (upperbindx_boot - lowerbindx_boot + 1) & sum(HW_CB_lbs <= 0) != (upperbindx_boot - lowerbindx_boot + 1)
    if (g_S_all | eq_S_all) {
        test_nocross = 1
    }
    
    ########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
    inttest_pre_dF     <- teststat_pre*diff(c(0, 1 - fit$surv[fit_time_restrict_boot]))[lowerbindx_boot:upperbindx_boot]
    inttest_pre_dbarNt <- teststat_pre*diff(c(0, barNt))[lowerbindx_boot:upperbindx_boot]
    inttest_pre_dt     <- teststat_pre[-(upperbindx_boot - lowerbindx_boot + 1)]*diff(Td_sort_boot)
    
    suptest        <- max(teststat_pre)
    inttest_dF     <- sum(inttest_pre_dF)
    inttest_dbarNt <- sum(inttest_pre_dbarNt)
    inttest_dt     <- sum(inttest_pre_dt)
    ########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
    
    if (details.return == TRUE) {
        return (list(
        test_nocross = test_nocross,
        neg2ELratio_at_ts = teststat_pre,
        wt_dbarNt = wt_dbarNt,
        wt_dt = wt_dt,
        wt_dF = wt_dF,
        neg2ELratio_bootstrap_at_ts = Up2_boot_H1_1sided,
        wt_dbarNt_boot = wt_dbarNt_boot,
        wt_dt_boot = wt_dt_boot,
        wt_dF_boot = wt_dF_boot,
        lowerbindx_boot = lowerbindx_boot,
        upperbindx_boot = upperbindx_boot,
        ts = Td_sort_boot,
        bootstrap_ts = fit$time[fit_time_restrict_boot],
        
        ########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
        EL_SOcrit           = EL_SOcrit,
        int_dFEL_SOcrit     = int_dFEL_SOcrit,
        int_dbarNtEL_SOcrit = int_dbarNtEL_SOcrit,
        int_dtEL_SOcrit     = int_dtEL_SOcrit,
        
        suptest        = suptest,
        inttest_dF     = inttest_dF,
        inttest_dbarNt = inttest_dbarNt,
        inttest_dt     = inttest_dt,
        
        p_value_suptest        = mean(sup_boot_H1        > suptest),
        p_value_inttest_dF     = mean(int_dF_boot_H1     > inttest_dF),
        p_value_inttest_dbarNt = mean(int_dbarNt_boot_H1 > inttest_dbarNt),
        p_value_inttest_dt     = mean(int_dt_boot_H1     > inttest_dt)
        ########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
        ))
    } else return (list(neg2ELratio_at_ts = teststat_pre))
}
