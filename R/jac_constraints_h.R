
jac_constraints_h = function(h, k, iter, iter_without_last, d, r, A, D, M_vec) {
  lgth_d = unlist(lapply(iter, FUN =function(iter){
    length(d[[iter]])
  }))
  hh = lapply(iter, FUN =function(iter){
    h[ (1 + cumsum(lgth_d[1:k])[iter] - lgth_d[iter])  : cumsum(lgth_d[1:k])[iter] ]
  })
  out2D = lapply(iter_without_last,FUN =function(iter_without_last){
    out = c()
    if (iter_without_last != 1){
      out = c(out, rep(0, length = cumsum(lgth_d[1:k])[iter_without_last-1]))
    }
    out = c(out,division00(M_vec[iter_without_last], 1 - hh[[iter_without_last]]), division00(-M_vec[iter_without_last+1], 1 - hh[[iter_without_last+1]]))
    if (iter_without_last != k-1){
      out = c(out,rep(0, length = cumsum(lgth_d[1:k])[k] - cumsum(lgth_d[1:k])[iter_without_last+1]))
    }
    return(out)
  })
  rslt = c()
  for (i in iter_without_last){
    rslt = rbind(rslt, out2D[[i]])
  }
  return(rslt)
}
