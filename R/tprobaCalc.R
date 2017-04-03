tprobaCalc <- function(moy, std, n, dl, em) {
    .Call(C_tproba, moyC = as.double(moy),
                 stdC = as.double(std),
                 nC = as.integer(n),
                 dlC = as.integer(dl),
                 emC = as.double(em));
}
