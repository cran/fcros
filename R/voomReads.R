voomReads <- function(x, Rm=1e+06) {
    n <- nrow(x);
    m <- ncol(x);
    x2 <- matrix(c(rep(0, m*n)), ncol = m);
    colnames(x2) <- colnames(x);
    Rj <- apply(x, 2, sum);
    for (j in 1:m) {
        x2[,j] <- log2(Rm*((x[,j] + 0.5)/(Rj[j] + 1)));
    }
    return(x2);
}
