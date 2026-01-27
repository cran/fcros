tprobaCalc <- function(moy, std, n, dl, em) {
    .Call(C_tproba, moy = as.double(moy),
                 stdC = as.double(std),
                 n = as.integer(n),
                 dl = as.integer(dl),
                 em = as.double(em));
}
