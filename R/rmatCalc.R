rmatCalc <- function(fvect, n, m1, m2) {
    .Call(C_rmat, fvect = as.double(fvect),
               n = as.integer(n),
               m1 = as.integer(m1),
               m2 = as.integer(m2));
}
