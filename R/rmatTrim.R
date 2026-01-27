rmatTrim <- function(rvect, n, m, idx, m2) {
    .Call(C_rmat2, rvect = as.double(rvect),
               n = as.integer(n),
               m = as.integer(m),
               idx = as.integer(idx),
               m2 = as.integer(m2));
}
