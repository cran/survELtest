
division00 <- function(x, y) {
  out <- x/y
  out[as.logical((x == 0)*(x == y))] <- 1
  return (out)  
}