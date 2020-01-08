#' @importFrom stats terms
checkInput_formula_data = function(fun, formula, data){
  #check if the set of formula and data is correctly expressed
  if(is.null(data) && !is.na(match("data", all.names(formula)))){
    stop("Do not include any vector called data in formula, as it may conflict with the second data argument. \n  If the data frame you mention in formula is called data, please assign it to the second data argument \n  (i.e., data = data).")
  }
  if(length(labels(terms(formula))) != 1){
    stop("Only one variable is allowed to determine group.")
  }
  if((sum(charToRaw(deparse(formula, width.cutoff = 500L)) == charToRaw('$')) == 0) && (is.null(data))){
  }else if(sum(charToRaw(deparse(formula, width.cutoff = 500L)) == charToRaw('$')) >= 3){
    if( (sum(charToRaw(as.character(formula[[2]])[2]) == charToRaw('$')) < 1) || (sum(charToRaw(as.character(formula[[2]])[3]) == charToRaw('$')) < 1) || (sum(charToRaw(labels(terms(formula))) == charToRaw('$')) < 1) ){
      stop(paste("Check the input formula and data. Three accepted sets of formula and data are: \n  1. ", fun, "(Surv(your_time, your_censor) ~ your_group) \n  2. ", fun, "(Surv(your_data$your_time, your_data$your_censor) ~ your_data$your_group) \n  3. ", fun, "(Surv(your_data$your_time, your_data$your_censor) ~ your_data$your_group, your_data)", sep = ""))
    }
  }else{
    stop(paste("Check the input formula and data. Three accepted sets of formula and data are: \n  1. ", fun, "(Surv(your_time, your_censor) ~ your_group) \n  2. ", fun, "(Surv(your_data$your_time, your_data$your_censor) ~ your_data$your_group) \n  3. ", fun, "(Surv(your_data$your_time, your_data$your_censor) ~ your_data$your_group, your_data)", sep = ""))
  }
}
