#' @importFrom stats terms
processInput_group_order = function(formula, group_order){
  #set group_order if it is null
  if(is.null(group_order)){
    if(sum(charToRaw(labels(terms(formula))) == charToRaw('$')) == 0){
      group_order = sort(unique(eval(get(all.vars(formula)[3]))))
    }else if(sum(charToRaw(labels(terms(formula))) == charToRaw('$')) >= 1){
      group_order = sort(unique(eval(parse(text = labels(terms(formula))))))
    }
  }
  return(group_order)
}
