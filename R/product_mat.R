
product_mat <- function(x, y) {
    out <- x*y
    for (i in 1:dim(x)[1]) {
        out[i, (x[i, ] == 0 | y[i, ] == 0)] <- 0
    }
    return (out)
}
