tcnReads <- function(x, maxVal=0) {
    n <- nrow(x);
    m <- ncol(x);
    x2 <- matrix(c(rep(0, m*n)), ncol = m);
    colnames(x2) <- colnames(x);
    sampReads <- apply(x, 2, sum);
    if (maxVal <= 1000) maxVal <- median(sampReads)
    for (j in 1:m) x2[,j] <- (x[,j]/sampReads[j])*maxVal
    x2 <- round(x2,0)
    return(x2);
}
