moyStdCalc <- function(rvect, n, m) {
    .Call(C_moyStd, rvect = as.double(rvect),
               n = as.integer(n),
               m = as.integer(m));
}
