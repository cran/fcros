moyStdCalc <- function(rvect, n, m) {
    .Call(C_moyStd, rvectC = as.double(rvect),
               nC = as.integer(n),
               mC = as.integer(m));
}
