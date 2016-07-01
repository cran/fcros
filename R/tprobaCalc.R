tprobaCalc <- function(moy, std, n, dl, em, proba) {
    .C("tproba", moyC = as.double(moy),
                 stdC = as.double(std),
                 nC = as.integer(n),
                 dlC = as.integer(dl),
                 emC = as.double(em),
                 probaC = as.double(proba));
}
