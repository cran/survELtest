
fnc_objective_h = function(h, k, iter, iter_without_last, d, r, A, D, M_vec) {
  lgth_d = unlist(lapply(iter, FUN =function(iter){
    length(d[[iter]])
  }))
  hh = lapply(iter, FUN =function(iter){
    h[ (1 + cumsum(lgth_d[1:k])[iter] - lgth_d[iter])  : cumsum(lgth_d[1:k])[iter] ]
  })
  sum_elmts = unlist(lapply(iter, FUN =function(iter){
    c(sum(d[[iter]]*log(hh[[iter]])),sum(product(A[[iter]], log(1-hh[[iter]]))))
  }))
  out = sum(sum_elmts)
  return (-out)
}
