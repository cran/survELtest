
constraints_h = function(h, k, iter, iter_without_last, d, r, A, D, M_vec) {
  lgth_d = unlist(lapply(iter, FUN =function(iter){
    length(d[[iter]])
  }))
  hh = lapply(iter, FUN =function(iter){
    h[ (1 + cumsum(lgth_d[1:k])[iter] - lgth_d[iter])  : cumsum(lgth_d[1:k])[iter] ]
  })
  out = unlist(lapply(iter_without_last, FUN =function(iter_without_last){
    log(division00(prod((1 - hh[[iter_without_last+1]]) ^ M_vec[iter_without_last+1]), prod((1 - hh[[iter_without_last]]) ^ M_vec[iter_without_last])))
  }))
  return(out)
}
