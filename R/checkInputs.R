
checkInputs = function(data, group_order, t1, t2, sided, nboot, wt, alpha, seed, nlimit){
  #fool-proof checks
  if (!all(sort(unique(data$censor)) == c(0,1))){
    stop("Censoring indicators can only be 0 or 1, and censoring indicators must include at least one 0 and one 1.")
  }
  if (!all(sort(group_order) == sort(unique(data$group)))){
    stop("Parameter \"group_order\" doesn't match the actual grouping numbers in asssigned dataset.")
  }
  if (!(t1 <= max(data$time))){
    stop("Parameter \"t1\" can not be greater than maximum time in assigned dataset.")
  }
  if (!(t2 >= min(data$time))){
    stop("Parameter \"t2\" can not be smaller than minimum time in assigned dataset.")
  }
  if (!(t1 <= t2)){
    stop("Parameter \"t1\" should be smaller or equal than parameter \"t2\".")
  }
  if (!(sided %in% c(1, 2))){
    stop("Parameter \"sided\" can only be 1 or 2.")
  }
  if (!is.numeric(nboot)){
    stop("Parameter \"nboot\" should be a number.")
  }
  if(!is.null(wt)){
    if (!(wt %in% c("p.event","dF", "dt"))){
    stop("Parameter \"wt\" can only be \"p.event\", \"dF\" or \"dt\".")
    }
  }
  if (!(0 <= alpha & alpha<= 1)){
    stop("Parameter \"alpha\" should be between 0 and 1.")
  }
  if (!is.numeric(seed)){
    stop("Parameter \"seed\" should be a number.")
  }
  if (!is.numeric(nlimit)){
    stop("Parameter \"nlimit\" should be a number.")
  }
}
