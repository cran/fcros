moyStdCalc <- function(rvect, n, m, moy, std) {
    .C("moyStd", rvectC = as.double(rvect),
               nC = as.integer(n),
               mC = as.integer(m),
               moyC = as.double(moy),
               stdC = as.double(std));
}
