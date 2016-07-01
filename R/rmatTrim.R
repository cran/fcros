rmatTrim <- function(rvect, n, m, idx, m2, rvect2) {
    .C("rmat2", rvectC = as.double(rvect),
               nC = as.integer(n),
               mC = as.integer(m),
               idxC = as.integer(idx),
               m2C = as.integer(m2),
               rvect2C = as.double(rvect2));
}
