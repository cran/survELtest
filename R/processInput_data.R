#' @importFrom stats terms
processInput_data = function(formula, data){
  #create data if it is null
  if((sum(charToRaw(deparse(formula, width.cutoff = 500L)) == charToRaw('$')) == 0) && (is.null(data))){
    data = data.frame(time = get(all.vars(formula)[1]), censor = get(all.vars(formula)[2]), group = get(all.vars(formula)[3]))
  }else if((sum(charToRaw(deparse(formula, width.cutoff = 500L)) == charToRaw('$')) >= 3) && (is.null(data))){
    #find time, censor and group automatically
    time = eval(parse(text = as.character(formula[[2]])[2]))
    group = eval(parse(text = labels(terms(formula))))
    censor = eval(parse(text = as.character(formula[[2]])[3]))
    data = data.frame(time = time, censor = censor, group = group)
  }else if((sum(charToRaw(deparse(formula, width.cutoff = 500L)) == charToRaw('$')) >= 3) && (!is.null(data))){
    time = eval(parse(text = as.character(formula[[2]])[2]))
    group = eval(parse(text = labels(terms(formula))))
    censor = eval(parse(text = as.character(formula[[2]])[3]))
    data = data.frame(time = time, censor = censor, group = group)
  }
  return(data)
}
