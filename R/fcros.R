fcros <- function(xdata, cont, test, log2.opt=0, trim.opt=0.25) {
    n <- nrow(xdata);
    idnames <- xdata[,1];   # first column is unique ID for genes
    xcol <- colnames(xdata);
    n.xcol <- length(xcol);
    idx1 <- xcol %in% cont;
    m1 <- sum(idx1 == TRUE);
    idx2 <- xcol %in% test;
    m2 <- sum(idx2 == TRUE);
    m <- m1+m2;
    m1m2 <- m1*m2;

    # form fcros data matrix
    fmat <- matrix(c(rep(0,n*m)), ncol = m);
    x1 <- matrix(c(rep(0,n*m1)), ncol = m1);
    x2 <- matrix(c(rep(0,n*m2)), ncol = m2);
    if (log2.opt) {
         x1 <- log2(xdata[, idx1 == TRUE]);
    } else {
         x1 <- xdata[, idx1 == TRUE];
    }
    fmat[,1:m1] <- as.matrix(x1);
    if (log2.opt) {
         x2 <- log2(xdata[, idx2 == TRUE]);
    } else {
         x2 <- xdata[, idx2 == TRUE];
    }
    fmat[,(m1+1):m] <- as.matrix(x2);

    # compute matrix containing pairwise fold changes
    rvect <- c(rep(0,n*m1m2));
    FC <- matrix(c(rep(0,n)), ncol = 1);
    fvect <- c(fmat[, 1:m])

    rmat.val <- rmatCalc(fvect, n, m1, m2);
    rmat <- matrix(rmat.val$rvectC, ncol = m1m2);
    FC <- rmat.val$FCC;
    rvectC <- rmat.val$rvectC;

    # compute the standard ranks matrix
    rmat.s <- apply(rmat, 2, rank, ties.method = "average")/n;

   # if (trim.opt), reduce the standard rank matrix
    if ((trim.opt > 0) & (trim.opt < 0.5)) {
       deb <- round(trim.opt * m1m2) + 1;
       fin <- m1m2 - deb + 1;
       idx <- deb:fin;
       m2 <- length(idx);
       rvect <- c(rmat.s[,1:m1m2]);
       rvect2 <- rmatTrim(rvect, n, m1m2, idx, m2);
       rmat.sr <- matrix(rvect2, ncol = m2);
       rvect <- rvect2;
       rmat.val <- moyStdCalc(rvect, n, m2);
       moyV <- rmat.val$moyC;
       stdV <- rmat.val$stdC;
       FC2 <- fc2Calc(rvectC, n, m1m2, idx, m2)
       rm(rvect);
       rm(rmat.val);
       rm(rvectC);
    }
    else {
         rmat.sr <- rmat.s;
         m2 <- m1m2;
         idx <- 1:m2;
         rvect <- c(rmat.sr[,1:m2]);
         rmat.val <- moyStdCalc(rvect, n, m2);
         moyV <- rmat.val$moyC;
         stdV <- rmat.val$stdC;
         FC2 <- fc2Calc(rvectC, n, m1m2, idx, m2)
         rm(rvect);
         rm(rmat.val);
         rm(rvectC);
    }

    # compute averages ranks
    ri <- apply(rmat.sr, 1, mean);
    ris <- sort(ri);

    # compute parameters
    lb <- n*ris[1];
    ub <- n*ris[n];
    delta <- (n-1)*mean(ris[-1]-ris[-n]);

    # compute f-value
    moy <- mean(ri);
    std <- sqrt((n-1)/n)*sd(ri);
    f.value <- pnorm(ri, mean = moy, sd = std);

    # perform the Student one sample test to get p-values
    em <- 0.5;
    p.value <- tprobaCalc(moyV, stdV, n, m2-1, em);

    moy_t <- (lb+ub)/(2*n);
    delta_t <- (ub-lb)/(n-1);
    std_t <- delta_t/sqrt(12);
    bounds <- c(lb,ub);
    params <- c(delta,moy,std);
    params_t <- c(delta_t,moy_t,std_t);

    list(idnames=idnames, FC=FC, FC2=FC2, ri=ri, f.value=f.value,
    p.value=p.value, bounds=bounds, params=params, params_t=params_t);
}
