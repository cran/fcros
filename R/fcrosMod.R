fcrosMod <- function(fcMat, samp, log2.opt = 0, trim.opt = 0.25) {
    # compute matrix of sorted FC ranks
    n <- nrow(fcMat);
    idnames <- fcMat[,1];      # first column is ID names
    xcol <- colnames(fcMat);
    n.xcol <- length(xcol);
    idx <- xcol %in% samp;
    m <- sum(idx)
    if (log2.opt) {
       fc <- log2(fcMat[, idx == TRUE]);
    } else {
       fc <- fcMat[, idx == TRUE];
    }
    fc.mat <- matrix(c(rep(0,n*m)), ncol = m);
    fc.mat[,1:m] <- as.matrix(fc);
    rvectC <- c(fc.mat[,1:m]);

    rmat <- apply(fc.mat, 2, rank, ties.method = "average")/n

    # if trim.opt > 0, reduce sorted rank matrix
    if ((trim.opt > 0) & (trim.opt < 0.5)) {
       deb <- round(trim.opt * m) + 1;
       fin <- m - deb + 1;
       idx <- deb:fin;
       m2 <- length(idx);
       rvect <- c(rmat[,1:m]);
       rvect2 <- rmatTrim(rvect, n, m, idx, m2);
       rmat.sr <- matrix(rvect2, ncol = m2);
       rvect <- rvect2;
       rmat.val <- moyStdCalc(rvect, n, m2);
       moyV <- rmat.val$moyC;
       stdV <- rmat.val$stdC;
       FC2 <- fc2Calc(rvectC, n, m, idx, m2)
       rm(rvect);
       rm(rmat.val);
       rm(rvectC);
    }
    else {
         rmat.sr <- rmat;
         m2 <- m;
         idx <- 1:m2;
         rvect <- c(rmat.sr[,1:m2]);
         rmat.val <- moyStdCalc(rvect, n, m2);
         moyV <- rmat.val$moyC;
         stdV <- rmat.val$stdC;
         FC2 <- fc2Calc(rvectC, n, m, idx, m2)
         rm(rvect);
         rm(rmat.val);
         rm(rvectC);
    }

    # compute vectors of f-values and p-values
    ri <- apply(rmat.sr, 1, mean);
    ris <- sort(ri);

    # compute parameters
    lb <- n*ris[1];
    ub <- n*ris[n];
    delta <- (n-1)*mean(ris[-1]-ris[-n]);

    moy <- mean(ri);
    std <- sqrt((n-1)/n)*sd(ri);
    f.value <- pnorm(ri, moy, std);

    # perform the Student one sample test
    em <- 0.5;
    p.value <- tprobaCalc(moyV, stdV, n, m2-1, em);

    moy_t <- (lb+ub)/(2*n);
    delta_t <- (ub-lb)/(n-1);
    std_t <- delta_t/sqrt(12);
    bounds <- c(lb,ub);
    params <- c(delta,moy,std);
    params_t <- c(delta_t,moy_t,std_t);
    if (log2.opt) { FC2 <- apply(fc.mat, 1, mean, trim = trim.opt); }
    else {FC2 <- apply(2^fc.mat, 1, mean, trim = trim.opt);}

    list(idnames = idnames, FC2 = FC2, ri = ri, p.value = p.value, f.value = f.value,
    bounds = bounds, params = params, params_t = params_t)
}
