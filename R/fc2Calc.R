fc2Calc <- function(rvect, n, m, idx, m2, fc2) {
    .C("fc2", rvectC = as.double(rvect),
              nC = as.integer(n),
              mC = as.integer(m),
              idxC = as.integer(idx),
              m2C = as.integer(m2),
              fc2C = as.double(fc2));
}
