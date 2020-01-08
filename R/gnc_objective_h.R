
gnc_objective_h = function(h, k, iter, iter_without_last, d, r, A, D, M_vec) {
  lgth_d = unlist(lapply(iter, FUN =function(iter){
    length(d[[iter]])
  }))
  hh = lapply(iter, FUN =function(iter){
    h[ (1 + cumsum(lgth_d[1:k])[iter] - lgth_d[iter])  : cumsum(lgth_d[1:k])[iter] ]
  })
  out = unlist(lapply(iter, FUN =function(iter){
    division00(d[[iter]], hh[[iter]]) - division00(r[[iter]] - d[[iter]], 1 - hh[[iter]])
  }))
  return(-out)
}
