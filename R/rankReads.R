rankReads <- function(xdata, cont, test, nrun=10, err=0.01, trim.opt=0.25, rseed=57) {
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

    # adjust uniform value level if necessary
    if ((err < 0) || (err > 1)) err <- 0.01;

    # form data matrix
    fmat <- matrix(c(rep(0,n*m)), ncol = m);
    x1 <- xdata[, idx1 == TRUE];
    fmat[,1:m1] <- as.matrix(x1);
    x2 <- xdata[, idx2 == TRUE];
    fmat[,(m1+1):m] <- as.matrix(x2);

    set.seed(rseed);
    # compute rank statistics
    if (nrun <= 0) {
        fmat3 <- log2(fmat + matrix(c(runif(n*m, 0, err)), ncol=m))
        fvect <-  c(fmat3)
        rvect <- rmatCalc(fvect, n, m1, m2);
        rmat <- matrix(rvect$rvectC, ncol=m1*m2);
        rmat2 <- apply(rmat, 2, rank, ties.method = "average")/n;
        ri <- apply(rmat2, 1, mean, trim.opt);
        FC  <- rvect$FCC;
        FC2 <- apply(2^rmat, 1, mean, trim.opt);
        score <- rep(1, n);
        deb <- round(trim.opt*m1m2) + 1;
        fin <- m1m2 - deb + 1;
        idx <- deb:fin;
        rmat3 <- t(apply(rmat2, 1, sort));
        rmat4 <- rmat3[, idx];
        ri.std <- apply(rmat4, 1, sd);
    } else {
        for (r in 1:nrun) {
            fmat3 <- log2(fmat + matrix(c(runif(n*m, 0, err)), ncol=m))
            fvect <-  c(fmat3)
            rvect <- rmatCalc(fvect, n, m1, m2);
            rmat <- matrix(rvect$rvectC, ncol=m1*m2);
            rmat2 <- apply(rmat, 2, rank, ties.method = "average")/n;
            if (r==1) {
               riMat <- apply(rmat2, 1, mean, trim.opt);
               fcMat <- rvect$FCC;
               fcMat2 <- apply(2^rmat, 1, mean, trim.opt)
            } else {
               ri <- apply(rmat2, 1, mean, trim.opt);
               fc2 <- apply(2^rmat, 1, mean, trim.opt);
               riMat <- cbind(riMat, ri);
               fcMat <- cbind(fcMat, rvect$FCC);
               fcMat2 <- cbind(fcMat2, fc2);
            }
        }
        ri <- apply(riMat, 1, mean);
        ri.std <- apply(riMat, 1, sd);
        mt <- m1+m2-2*round(trim.opt*(m1+m2));
        yi <- 12*(mt)*ri.std^2;
        score <- pf(yi, mt-1, mt-1);
        FC <- apply(fcMat, 1, mean);
        FC2 <- apply(fcMat2, 1, mean);
    }
    # sort ranks average values
    ris <- sort(ri);

    # compute parameters
    lb <- n*ris[1];
    ub <- n*ris[n];
    delta <- (n-1)*mean(ris[-1]-ris[-n]);

    # compute f-value
    moy <- mean(ri);
    std <- sqrt((n-1)/n)*sd(ri);
    f.value <- pnorm(ri, moy, std);

    # perform the Student one sample test to get p-values
    mt <- round((1-2*trim.opt)*m1m2);
    em <- 0.5;
    p.value <- tprobaCalc(ri, ri.std, n, mt-1, em);

    moy_t <- (lb+ub)/(2*n);
    delta_t <- (ub-lb)/(n-1);
    std_t <- delta_t/sqrt(12);
    bounds <- c(lb, ub);
    params <- c(delta, moy, std);
    params_t <- c(delta_t, moy_t, std_t);

    list(idnames=idnames, FC=FC, FC2=FC2, ri=ri, f.value=f.value,
    p.value=p.value, bounds=bounds, params=params, params_t=params_t, 
    score=score);
}
