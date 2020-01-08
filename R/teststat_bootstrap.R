#' @importFrom stats rnorm
teststat_bootstrap <- function(nlimit, k, iter, fit, fitls, NBOOT, mm, nn, fit_time_restrict_boot){
  # function for bootstrap
  ### computing bootstrap related quantities for each treatment group:
  nobd = unlist(lapply(iter, FUN=function(iter){
    length(rep(fitls[[iter]]$time, times = fitls[[iter]]$n.event))
  }))
  data_big =  lapply(iter, FUN= function(iter){
    matrix(rep(rep(fitls[[iter]]$time, times = fitls[[iter]]$n.event), times = mm), byrow = TRUE, nrow = mm)
  })
  fit_time_big = lapply(iter, FUN=function(iter){
    matrix(rep(fit$time[fit_time_restrict_boot], times = nobd[iter]), byrow = FALSE, ncol = nobd[iter])
  })
  Indt_big = lapply(iter, FUN=function(iter){
    data_big[[iter]] <= fit_time_big[[iter]]
  })
  Indx = lapply(iter, FUN=function(iter){
    rep(fitls[[iter]]$n.risk, times = fitls[[iter]]$n.event)
  })
  Indx_big = lapply(iter, FUN=function(iter){
    matrix(rep(Indx[[iter]], times = mm), byrow = TRUE, nrow = mm)
  })
  
  nsplit = ceiling(mm/nlimit)
  if (nsplit > 1){
    sum_DWknGImuw_big = lapply(iter, FUN=function(iter){
      array(0,c(NBOOT,mm))
    })
    for (i in 1:nsplit){
      if(i == nsplit){
        nboot = NBOOT - floor(NBOOT/nsplit)*(nsplit - 1)
        Gs = lapply(iter, FUN = function(iter){
          matrix(rnorm(nobd[iter] * nboot), nrow = nboot, ncol = nobd[iter])
        })
        Imuw_BIG = lapply(iter, FUN=function(iter){
          array(rep(Indt_big[[iter]] / Indx_big[[iter]], each = nboot), c(nboot, mm, nobd[iter]))
        })
        Gs_BIG = lapply(iter, FUN = function(iter){
          array(matrix(rep(t(Gs[[iter]]), each = mm), byrow = TRUE, nrow = nboot), c(nboot, mm, nobd[iter]))
        })
        temp_sum_DWknGImuw_big = lapply(iter, FUN = function(iter){
          sqrt(sum(nn)) * apply(Imuw_BIG[[iter]] * Gs_BIG[[iter]], c(1, 2), sum)
        })
        for (j in 1:k){
          sum_DWknGImuw_big[[j]][(NBOOT - nboot + 1):NBOOT,] = temp_sum_DWknGImuw_big[[j]]
        }
      }else{
        nboot = floor(NBOOT/nsplit)
        Gs = lapply(iter, FUN = function(iter){
          matrix(rnorm(nobd[iter] * nboot), nrow = nboot, ncol = nobd[iter])
        })
        Imuw_BIG = lapply(iter, FUN=function(iter){
          array(rep(Indt_big[[iter]] / Indx_big[[iter]], each = nboot), c(nboot, mm, nobd[iter]))
        })
        Gs_BIG = lapply(iter, FUN = function(iter){
          array(matrix(rep(t(Gs[[iter]]), each = mm), byrow = TRUE, nrow = nboot), c(nboot, mm, nobd[iter]))
        })
        temp_sum_DWknGImuw_big = lapply(iter, FUN = function(iter){
          sqrt(sum(nn)) * apply(Imuw_BIG[[iter]] * Gs_BIG[[iter]], c(1, 2), sum)
        })
        for (j in 1:k){
          sum_DWknGImuw_big[[j]][((i - 1)*nboot + 1):(i*nboot),] = temp_sum_DWknGImuw_big[[j]]
        }
      }
    }
  }else{
    Gs = lapply(iter, FUN = function(iter){
      matrix(rnorm(nobd[iter] * NBOOT), nrow = NBOOT, ncol = nobd[iter])
    })
    Imuw_BIG = lapply(iter, FUN=function(iter){
      array(rep(Indt_big[[iter]] / Indx_big[[iter]], each = NBOOT), c(NBOOT, mm, nobd[iter]))
    })
    Gs_BIG = lapply(iter, FUN = function(iter){
      array(matrix(rep(t(Gs[[iter]]), each = mm), byrow = TRUE, nrow = NBOOT), c(NBOOT, mm, nobd[iter]))
    })
    sum_DWknGImuw_big = lapply(iter, FUN = function(iter){
      sqrt(sum(nn)) * apply(Imuw_BIG[[iter]] * Gs_BIG[[iter]], c(1, 2), sum)
    })
  }
  return(list(sum_DWknGImuw_big = sum_DWknGImuw_big))
}
