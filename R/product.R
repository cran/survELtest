
product <- function(x, y) {
    out <- x*y
    out[(x == 0 | y == 0)] <- 0
    return (out)
}
