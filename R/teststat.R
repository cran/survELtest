#' @importFrom survival Surv
#' @importFrom survival survfit
#' @importFrom stats rnorm
#' @importFrom stats quantile
#' @importFrom stats median
#' @importFrom methods is
#' @importFrom Iso "pava"
#' @importFrom plyr "alply"


teststat = function(data, group_order, t1, t2, sided, nboot, alpha, details.return, seed, nlimit) {
    set.seed(seed)
    # note that the third column should be the group number.
    k = length(unique(data[,3]))
    M_vec = rep(1, k)
    dat = Surv(data[, 1], data[, 2])
    fit = survfit(dat ~ 1)
    fitls <- lapply(1:k, FUN = function(i) {
        survfit(dat[data[, 3] == group_order[i]] ~ 1)
    })
    NBOOT = nboot
    iter = 1 : k
    iter_without_last = 1:(k-1)
    k = length(fitls)
    nn = c()
    for (i in iter){
        nn = c(nn,fitls[[i]]$n)
    }
    ### finding t \in [t_1, t_2] that are ordered observed uncensored times
    # T_i1 is a vector withT_11,T_21,...T_{n}1
    T_i1 = unlist(lapply(iter, FUN= function(iter){
        min(fitls[[iter]]$time[fitls[[iter]]$n.event != 0])
    }))
    T_ms = 1:k * 0
    T_ms_all = 1:k * 0
    T_ms = unlist(lapply(iter, FUN = function(iter){
        max(fitls[[iter]]$time[fitls[[iter]]$n.event != 0])
    } ))
    T_ms_all = unlist(lapply(iter, FUN = function(iter){
        max(fitls[[iter]]$time)
    } ))
    fit_time_restrict_boot = (fit$n.event != 0 & fit$time <= t2 & fit$time >= t1)  # ordered observed uncensored times
    mm = length(fit$time[fit_time_restrict_boot])
    lowerb = max(T_i1)  # data driven t_1, as in Remark 3 of Theorem 1
    upperb = min(T_ms)  # data driven t_2, as in Remark 3 of Theorem 1
    upperb_which = which.min(T_ms)
    lowerbindx_boot = which.min(abs(fit$time[fit_time_restrict_boot] - lowerb))
    upperbindx_boot = which.min(abs(fit$time[fit_time_restrict_boot] - upperb))
    if (T_ms[upperb_which] == T_ms_all[upperb_which]) {
        upperbindx_boot = upperbindx_boot - 1
    }
    Td_sort_boot = (fit$time[fit_time_restrict_boot])[lowerbindx_boot:upperbindx_boot]  # t \in [t_1, t_2] that are ordered observed uncensored times
    ### Kaplan--Meier estimates for each treatment group and related quantities for the first step of our composite procedure:
    #  S_hat is a two-dimension with S1_hat, S2_hat,...,S{n}_hat
    S_hat = lapply(iter, FUN= function(iter){
        (((c(1, fitls[[iter]]$surv)[cumsum(c(0, fit$time) %in% c(0, fitls[[iter]]$time))])[-1])[fit_time_restrict_boot])[lowerbindx_boot:upperbindx_boot]
    })
    #  diff is a two-dimension with diff_12, diff_23,..., diff_{n-1}{n}
    diff = lapply(iter_without_last , FUN = function(iter_without_last){
        M_vec[iter_without_last] * log(S_hat[[iter_without_last]]) - M_vec[iter_without_last+1] * log(S_hat[[iter_without_last+1]])
    })
    ### computing bootstrap related quantities for each treatment group:
    # nobd is a vector with nobd1, nobd2,...,nobd{n}
    nobd = unlist(lapply(iter, FUN=function(iter){
        length(rep(fitls[[iter]]$time, times = fitls[[iter]]$n.event))
    }))
    # data_big is a list of matrix with data1_big, data2_big,...,data{n}_big
    data_big =  lapply(iter, FUN= function(iter){
        matrix(rep(rep(fitls[[iter]]$time, times = fitls[[iter]]$n.event), times = mm), byrow = TRUE, nrow = mm)
    })
    # fit_time_big is a list of matrix with fit_time_1_big,fit_time_2_big,...,fit_time_{n}_big
    fit_time_big = lapply(iter, FUN=function(iter){
        matrix(rep(fit$time[fit_time_restrict_boot], times = nobd[iter]), byrow = FALSE, ncol = nobd[iter])
    })
    # Indt_big is a list of matrix with Ind1t_big, Ind2t_big,...,Ind{n}t_big
    Indt_big = lapply(iter, FUN=function(iter){
        data_big[[iter]] <= fit_time_big[[iter]]
    })
    # Indx is a two-dimension list with Ind1x,Ind2x,...,Ind{n}x
    Indx = lapply(iter, FUN=function(iter){
        rep(fitls[[iter]]$n.risk, times = fitls[[iter]]$n.event)
    })
    # Indx_big is a list of matrix with Ind1x_big,Ind2x_big,...,Ind{n}x_big
    Indx_big = lapply(iter, FUN=function(iter){
        matrix(rep(Indx[[iter]], times = mm), byrow = TRUE, nrow = mm)
    })
    nsplit <- ceiling(mm/nlimit)
    
    if (nsplit > 1){
        sum_DWknGImuw_big = lapply(iter, FUN=function(iter){
            array(0,c(NBOOT,mm))
        })
        for (i in 1:nsplit){
            if(i == nsplit){
                nboot <- NBOOT - floor(NBOOT/nsplit)*(nsplit - 1)
                # Gs is a list of matrix with Gs_1,Gs_2,...,Gs_{n}
                Gs = lapply(iter, FUN = function(iter){
                    matrix(rnorm(nobd[iter] * nboot), nrow = nboot, ncol = nobd[iter])
                })
                # Imuw_BIG is a list of array with Imuw_BIG_1,Imuw_BIG_2,...,Imuw_BIG_{n}
                Imuw_BIG = lapply(iter, FUN=function(iter){
                    array(rep(Indt_big[[iter]] / Indx_big[[iter]], each = nboot), c(nboot, mm, nobd[iter]))
                })
                # Gs_BIG is a list of array with Gs_1_BIG,Gs_2_BIG,...,Gs_{n}_BIG
                Gs_BIG = lapply(iter, FUN = function(iter){
                    array(matrix(rep(t(Gs[[iter]]), each = mm), byrow = TRUE, nrow = nboot), c(nboot, mm, nobd[iter]))
                })
                temp_sum_DWknGImuw_big = lapply(iter, FUN = function(iter){
                    sqrt(sum(nn)) * apply(Imuw_BIG[[iter]] * Gs_BIG[[iter]], c(1, 2), sum)
                })
                # sum_DWknGImuw_big is a list of array with sum_DWknGImuw_1_big,sum_DWknGImuw_2_big,...,sum_DWknGImuw_{n}_big
                for (j in 1:k){
                    sum_DWknGImuw_big[[j]][(NBOOT - nboot + 1):NBOOT,] = temp_sum_DWknGImuw_big[[j]]
                }
            }else{
                nboot <- floor(NBOOT/nsplit)
                # Gs is a list of matrix with Gs_1,Gs_2,...,Gs_{n}
                Gs = lapply(iter, FUN = function(iter){
                    matrix(rnorm(nobd[iter] * nboot), nrow = nboot, ncol = nobd[iter])
                })
                # Imuw_BIG is a list of array with Imuw_BIG_1,Imuw_BIG_2,...,Imuw_BIG_{n}
                Imuw_BIG = lapply(iter, FUN=function(iter){
                    array(rep(Indt_big[[iter]] / Indx_big[[iter]], each = nboot), c(nboot, mm, nobd[iter]))
                })
                # Gs_BIG is a list of array with Gs_1_BIG,Gs_2_BIG,...,Gs_{n}_BIG
                Gs_BIG = lapply(iter, FUN = function(iter){
                    array(matrix(rep(t(Gs[[iter]]), each = mm), byrow = TRUE, nrow = nboot), c(nboot, mm, nobd[iter]))
                })
                # sum_DWknGImuw_big is a list of array with sum_DWknGImuw_1_big,sum_DWknGImuw_2_big,...,sum_DWknGImuw_{n}_big
                temp_sum_DWknGImuw_big = lapply(iter, FUN = function(iter){
                    sqrt(sum(nn)) * apply(Imuw_BIG[[iter]] * Gs_BIG[[iter]], c(1, 2), sum)
                })
                for (j in 1:k){
                    sum_DWknGImuw_big[[j]][((i - 1)*nboot + 1):(i*nboot),] = temp_sum_DWknGImuw_big[[j]]
                }
            }
        }
    }else{
        # Gs is a list of matrix with Gs_1,Gs_2,...,Gs_{n}
        Gs = lapply(iter, FUN = function(iter){
            matrix(rnorm(nobd[iter] * NBOOT), nrow = NBOOT, ncol = nobd[iter])
        })
        # Imuw_BIG is a list of array with Imuw_BIG_1,Imuw_BIG_2,...,Imuw_BIG_{n}
        Imuw_BIG = lapply(iter, FUN=function(iter){
            array(rep(Indt_big[[iter]] / Indx_big[[iter]], each = NBOOT), c(NBOOT, mm, nobd[iter]))
        })
        # Gs_BIG is a list of array with Gs_1_BIG,Gs_2_BIG,...,Gs_{n}_BIG
        Gs_BIG = lapply(iter, FUN = function(iter){
            array(matrix(rep(t(Gs[[iter]]), each = mm), byrow = TRUE, nrow = NBOOT), c(NBOOT, mm, nobd[iter]))
        })
        # sum_DWknGImuw_big is a list of array with sum_DWknGImuw_1_big,sum_DWknGImuw_2_big,...,sum_DWknGImuw_{n}_big
        sum_DWknGImuw_big = lapply(iter, FUN = function(iter){
            sqrt(sum(nn)) * apply(Imuw_BIG[[iter]] * Gs_BIG[[iter]], c(1, 2), sum)
        })
    }
    ### bootstrap components in the limiting distribution in Theorem 1:
    theta_hat_js = theta_hat_j_k(Td_sort_boot, fit, fitls, M_vec)
    sigma2_hat_overpjs = sigma2_hat_overpj_k(Td_sort_boot, fit, fitls, M_vec)
    # Dsqsigmadp_big is a list of matrix with Dsqsigmadp1_big,Dsqsigmadp2_big,...,Dsqsigmadp{n}_big
    Dsqsigmadp_big = lapply(iter, FUN = function(iter){
        matrix(rep(1 / sqrt(sigma2_hat_overpjs[, iter]), times = NBOOT), byrow = TRUE, nrow = NBOOT)
    })
    Ujps = array(0, c(k, NBOOT, upperbindx_boot - lowerbindx_boot + 1))
    for (i in iter){
        Ujps[i, , ] = product_mat(as.matrix(-sum_DWknGImuw_big[[i]][, lowerbindx_boot:upperbindx_boot]), Dsqsigmadp_big[[i]])
    }
    wjs = array(0, c(k, NBOOT, upperbindx_boot - lowerbindx_boot + 1))
    for (j in iter){
        wjs[j, , ] = matrix(rep(division00(1, theta_hat_js[, j] ), times = NBOOT), byrow = TRUE, nrow = NBOOT)
    }
    wjs_denom = 0
    for (j in iter){
        wjs_denom = wjs_denom + wjs[j,,]
    }
    for (j in iter){
        wjs[j, , ] =  wjs[j, , ] / wjs_denom
    }
    pava_time1 = Sys.time()
    pava_result = array(unlist(mapply(pava, alply(division00(Ujps, sqrt(wjs)), c(2, 3)), alply(wjs, c(2, 3)), decreasing = TRUE)), c(k, NBOOT, upperbindx_boot - lowerbindx_boot + 1))
    pava_time2 = Sys.time()
    avg_Ujps = apply(sqrt(wjs) * Ujps, c(2, 3), sum)
    avg_Ujps_karray = array(rep(avg_Ujps, each = k), c(k, NBOOT, upperbindx_boot - lowerbindx_boot + 1))
    Up2_boot_H1_1sided = apply(wjs * (pava_result - avg_Ujps_karray) ^ 2, c(2, 3), sum)  # bootstrap SSB(t) in Theorem 1
    ### bootstrap sup_{t \in [t_1, t_2]} SSB(t) in Theorem 1:
    sup_boot_H1 = apply(as.matrix(Up2_boot_H1_1sided), 1, max)  # bootstrap sup_{t \in [t_1, t_2]} SSB(t) in Theorem 1
    EL_SOcrit = as.vector(quantile(sup_boot_H1, probs = 1 - alpha))  # quantiles based on bootstrapped values of sup_{t \in [t_1, t_2]} SSB(t)
    
    nobd_pooled         <- length(rep(fit$time, times = fit$n.event))
    data_big_pooled     <- matrix(rep(rep(fit$time, times = fit$n.event), times = mm), byrow = TRUE, nrow = mm)
    fit_time_big_pooled <- matrix(rep(fit$time[fit_time_restrict_boot], times=nobd_pooled), byrow = FALSE, ncol = nobd_pooled)
    Indt_big_pooled     <- (data_big_pooled <= fit_time_big_pooled)
    barNt        <- apply(Indt_big_pooled, 1, sum) / sum(nn)
    wt_dbarNt    <- diff(c(0, barNt))
    
    barNt_big      <- matrix(rep(barNt, time = NBOOT), byrow = TRUE, nrow = NBOOT)
    wt_dbarNt_boot <- t(apply(cbind(0, barNt_big), 1, diff))
    
    wt_dt      <- diff(c(fit$time[fit_time_restrict_boot], Inf))  # we don't want the last one
    t_big      <- matrix(rep(fit$time[fit_time_restrict_boot], time = NBOOT), byrow = TRUE, nrow = NBOOT)
    wt_dt_boot <- t(apply(cbind(t_big, Inf), 1, diff))
    
    wt_dF      <- diff(c(0, 1 - fit$surv[fit_time_restrict_boot]))
    F_big      <- matrix(rep(1 - fit$surv[fit_time_restrict_boot], time = NBOOT), byrow = TRUE, nrow = NBOOT)
    wt_dF_boot <- t(apply(cbind(0, F_big), 1, diff))
    Up2_boot_H1_1sided_times_dF = Up2_boot_H1_1sided * t(apply(cbind(0,  F_big), 1, diff))[, lowerbindx_boot:upperbindx_boot]
    int_dF_boot_H1 = apply(as.matrix(Up2_boot_H1_1sided_times_dF), 1, sum)  # bootstrap \int_{t_1}^{t_2} SSB(t) dF_0(t) in Theorem 1
    int_dFEL_SOcrit = as.vector(quantile(int_dF_boot_H1, probs = 1 - alpha))  # quantiles based on bootstrapped values of \int_{t_1}^{t_2} SSB(t) dF_0(t) in Theorem 1
    
    Up2_boot_H1_1sided_times_dbarNt <- Up2_boot_H1_1sided*t(apply(cbind(0, barNt_big), 1, diff))[, lowerbindx_boot:upperbindx_boot] #
    int_dbarNt_boot_H1  <- apply(as.matrix(Up2_boot_H1_1sided_times_dbarNt), 1, sum) #
    int_dbarNtEL_SOcrit <- as.vector(quantile(int_dbarNt_boot_H1, 1 - alpha)) #
    
    Up2_boot_H1_1sided_times_dt <- Up2_boot_H1_1sided[,-(upperbindx_boot-lowerbindx_boot+1)]*t(apply(cbind(t_big, Inf), 1, diff))[, lowerbindx_boot:(upperbindx_boot - 1)] #
    int_dt_boot_H1  <- apply(as.matrix(Up2_boot_H1_1sided_times_dt), 1, sum) #
    int_dtEL_SOcrit <- as.vector(quantile(int_dt_boot_H1, 1 - alpha)) #
    
    
    ### functions related to the nonparametric likelihood ratio statistics:
    neg2logR = function(t, fitls, M_vec, EL_CBcrit = 0) {
        k = length(fitls)
        iter = 1:k
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
        for (i in 1:k){
            D[i] = max(d[[i]] - r[[i]]) / M_vec[i] + 0.0001
        }
        cumsum_neg_D = -cumsum(D[k:2])[(k - 1):1]
        init_lambdas = rep(0, times = k)
        for (i in 1:(k - 1)) {
            init_lambdas[i + 1] = median(c(D[i] + init_lambdas[i], cumsum_neg_D[i]))
        }
        init_lambdas = init_lambdas[2:k]
        num_neg_loglik   = neg_log_likelihood_le_t_eq_h_k  (t, fitls, M_vec, init_lambdas = init_lambdas)
        num_warn         = tryCatch(neg_log_likelihood_le_t_eq_h_k(t, fitls, M_vec, init_lambdas = init_lambdas),
        error = function(e) e, warning = function(w) w)
        if (sided == 1) {
            denom_neg_loglik = neg_log_likelihood_le_t_ineq_h_k(t, fitls, M_vec,init_lambdas = init_lambdas)
            denom_warn       = tryCatch(neg_log_likelihood_le_t_ineq_h_k(t, fitls, M_vec, init_lambdas = init_lambdas), error = function(e) e, warning = function(w) w)
        } else {
            denom_neg_loglik = neg_log_likelihood_le_t_glob_h_k(t, fitls, M_vec,init_lambdas = init_lambdas)
            denom_warn       = tryCatch(neg_log_likelihood_le_t_glob_h_k(t, fitls, M_vec, init_lambdas = init_lambdas), error = function(e) e, warning = function(w) w)
        }
        init_dir1 = 1
        bdry_r    = cumsum_neg_D[1] - 0.0001
        bdry_l    = D[1] + 0.0001
        init_wid1 = 0.1 * (bdry_r - bdry_l)
        while (sum(abs(num_neg_loglik$resultEEs) >= 0.0001)        != 0 |
        sum(denom_neg_loglik$resultEEs    >= 10^-8)         != 0 |
        is(num_warn, "warning") * is(denom_warn, "warning") != 0)
        {
            init_lambdas[1] = init_lambdas[1] + init_wid1 * init_dir1 * (-1)^init_dir1
            if (init_lambdas[1] > bdry_r | init_lambdas[1] < bdry_l) break
            if (init_lambdas[1] > bdry_r | init_lambdas[1] < bdry_l) break
            for (i in 2:(k - 1)) {
                init_lambdas[i] = median(c(D[i] + init_lambdas[i - 1], cumsum_neg_D[i]))
            }
            num_neg_loglik   = neg_log_likelihood_le_t_eq_h_k  (t, fitls, M_vec, init_lambdas = init_lambdas)
            denom_neg_loglik = neg_log_likelihood_le_t_ineq_h_k(t, fitls, M_vec, init_lambdas = init_lambdas)
            num_warn         = tryCatch(neg_log_likelihood_le_t_eq_h_k(t, fitls, M_vec, init_lambdas = init_lambdas),
            error = function(e) e, warning = function(w) w)
            if (sided == 1) {
                denom_neg_loglik = neg_log_likelihood_le_t_ineq_h_k(t, fitls, M_vec,init_lambdas = init_lambdas)
                denom_warn       = tryCatch(neg_log_likelihood_le_t_ineq_h_k(t, fitls, M_vec, init_lambdas = init_lambdas), error = function(e) e, warning = function(w) w)
            } else {
                denom_neg_loglik = neg_log_likelihood_le_t_glob_h_k(t, fitls, M_vec,init_lambdas = init_lambdas)
                denom_warn       = tryCatch(neg_log_likelihood_le_t_glob_h_k(t, fitls, M_vec, init_lambdas = init_lambdas), error = function(e) e, warning = function(w) w)
            }
            init_dir1 = init_dir1 + 1
        }
        still_error = (sum(abs(num_neg_loglik$resultEEs) >= 0.0001)        != 0 |
        sum(denom_neg_loglik$resultEEs    >= 10^-8)         != 0 |
        is(num_warn, "warning") * is(denom_warn, "warning") != 0)
        test1t      = 2 * (num_neg_loglik$objective - denom_neg_loglik$objective)
        return(list(out           = test1t - EL_CBcrit,
        still_error   = still_error,
        num_lambdas   = num_neg_loglik$lambdas,
        denom_lambdas = denom_neg_loglik$lambdas))
    }
    
    ### computing quantities related to the NPLR statistics:
    teststat_pre = 1:length(Td_sort_boot) * 0
    error_vec = 1:length(Td_sort_boot) * 0
    # teststat_pre is a two-dimension list with  teststat_pre_12,teststat_pre_13,...,teststat_pre_1{n},teststat_pre_23,...,teststat_pre_2{n},...
    for (j in 1:(upperbindx_boot - lowerbindx_boot + 1)) {
        t = Td_sort_boot[j]  # the given t >= 0 in localized EL statistic
        neg2logRt = neg2logR(t, fitls, M_vec, EL_CBcrit = 0)
        error_vec[j] = neg2logRt$still_error
        if (error_vec[j] == 1) next
        teststat_pre[j] = neg2logRt$out  # the -2 log NPLR statistic in Section 2.2 of the paper at the given t >=0 specified above
    }
    
   
    inttest_pre_dF = teststat_pre * wt_dF[lowerbindx_boot:upperbindx_boot]  # -2logR(t) d\hat{F}_0(t) in (3.1) of the paper
    inttest_pre_dbarNt <- teststat_pre*wt_dbarNt[lowerbindx_boot:upperbindx_boot]
    inttest_pre_dt     <- teststat_pre[-(upperbindx_boot - lowerbindx_boot + 1)]*diff(Td_sort_boot)
    
    ### computing quantities related to the first step of our composite procedure:
    test_nocross = 0
    # pmax_elmts is a two-dimension list used in next parameter HW_CBdistr with apply(abs(M_vec[i] * sum_DWknGImuw_{i}_big - M_vec[i]
    pmax_elmts = lapply(iter_without_last, FUN=function(iter_without_last){
        apply(abs(M_vec[iter_without_last+1] * sum_DWknGImuw_big[[iter_without_last+1]] - M_vec[iter_without_last] * sum_DWknGImuw_big[[iter_without_last]]), 1, max)
    })
    HW_CBdistr = pmax_elmts[[1]]
    if (k >2){
        for (i in 2:(k-1)){
            HW_CBdistr = pmax(HW_CBdistr, pmax_elmts[[i]])
        }
    }
    HW_CBcrit = as.vector(quantile(HW_CBdistr, 1-alpha))
    # HW_CB_ubs_pairwise is a two-dimension list with HW_CB_ubs_12,HW_CB_ubs_23,...HW_CB_ubs{n-1}{n}
    HW_CB_ubs_pairwise = lapply(iter_without_last, FUN = function(iter_without_last){
        diff[[iter_without_last]] + HW_CBcrit / sqrt(sum(nn))
    })
    # HW_CB_lbs_pairwise is a two-dimension list with HW_CB_lbs_12,HW_CB_lbs_23,...HW_CB_lbs_{n-1}{n}
    HW_CB_lbs_pairwise = lapply(iter_without_last, FUN =function(iter_without_last){
        diff[[iter_without_last]] - HW_CBcrit / sqrt(sum(nn))
    })
    g_S_all = sum(HW_CB_ubs_pairwise[[1]] >= 0) == (upperbindx_boot - lowerbindx_boot + 1) & sum(HW_CB_lbs_pairwise[[1]] <= 0) != (upperbindx_boot - lowerbindx_boot + 1)
    if (k >2){
        for (i in 2:(k-1)){
            g_S_all = g_S_all & sum(HW_CB_ubs_pairwise[[i]] >= 0) == (upperbindx_boot - lowerbindx_boot + 1) & sum(HW_CB_lbs_pairwise[[i]] <= 0) != (upperbindx_boot - lowerbindx_boot + 1)
        }
    }
    eq_S_all = sum(HW_CB_lbs_pairwise[[1]] <= 0) == (upperbindx_boot - lowerbindx_boot + 1) & sum(HW_CB_ubs_pairwise[[1]] >= 0) == (upperbindx_boot - lowerbindx_boot + 1)
    if (k > 2){
        for (i in 2: (k-1)){
            eq_S_all = eq_S_all & sum(HW_CB_lbs_pairwise[[i]] <= 0) == (upperbindx_boot - lowerbindx_boot + 1) & sum(HW_CB_ubs_pairwise[[i]] >= 0) == (upperbindx_boot - lowerbindx_boot + 1)
        }
    }
    if (g_S_all | eq_S_all) {
        test_nocross = 1
    }
    
    
    ########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
    inttest_pre_dF     <- teststat_pre*diff(c(0, 1 - fit$surv[fit_time_restrict_boot]))[lowerbindx_boot:upperbindx_boot]
    inttest_pre_dbarNt <- teststat_pre*diff(c(0, barNt))[lowerbindx_boot:upperbindx_boot]
    inttest_pre_dt     <- teststat_pre[-(upperbindx_boot - lowerbindx_boot + 1)]*diff(Td_sort_boot)
    
    suptest        <- max(teststat_pre) # same as test = max(teststat_pre) in original teststat
    inttest_dF     <- sum(inttest_pre_dF)
    inttest_dbarNt <- sum(inttest_pre_dbarNt)
    inttest_dt     <- sum(inttest_pre_dt)
    ########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
    
    if (details.return == TRUE) {
        return (list(
        test_nocross = test_nocross, # output that do not exist in neg2ELratio
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
        
        p_value_suptest        = mean(sup_boot_H1        > suptest), # same as out_sup_pval = mean(sup_boot_H1 >= max(teststat_pre)) in original teststat
        p_value_inttest_dF     = mean(int_dF_boot_H1     > inttest_dF), # same as out_dF_pval = mean(int_dF_boot_H1 >= sum(inttest_pre_dF)) in original teststat
        p_value_inttest_dbarNt = mean(int_dbarNt_boot_H1 > inttest_dbarNt),
        p_value_inttest_dt     = mean(int_dt_boot_H1     > inttest_dt)
        ########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
        ))
    } else return (list(neg2ELratio_at_ts = teststat_pre))
}

