#' @importFrom stats as.formula
processInput_formula = function(formula, data){
  #reformulate formula to a complete form
  if((sum(charToRaw(deparse(formula, width.cutoff = 500L)) == charToRaw('$')) == 0) && (is.null(data))){
    formula = as.formula(paste("Surv(data$", all.vars(formula)[1],",data$", all.vars(formula)[2], ")~data$", all.vars(formula)[3], sep = ""))
  }
  return(formula)
}
