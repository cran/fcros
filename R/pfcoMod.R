pfcoMod <- function(fcMat, samp, log2.opt=0, trim.opt=0.25) {
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

    # Compute sorted ranks matrix
    rmat.s <- (apply(fc.mat, 2, rank, ties.method = "average"))/n;

    # if trim.opt > 0, reduce sorted rank matrix
    if ((trim.opt > 0) & (trim.opt < 0.5)) {
       deb <- round(trim.opt * m) + 1;
       fin <- m - deb + 1;
       idx <- deb:fin;
       m2 <- length(idx);
       rvect <- c(rmat.s[,1:m]);
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
         rmat.sr <- rmat.s;
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

    # compute the symmetric matrix with sorted rank values and its eigen values
    smat <- (t(rmat.sr) %*% rmat.sr)/n;
    smat.eig <- eigen(smat);
    v <- smat.eig$vectors[,1:2];
    if (v[1,1] < 0) v <- -v;
    u1 <- mean(v[,1])* rmat.sr %*% v[,1];
    u2 <- mean(v[,2])* rmat.sr %*% v[,2];
    u1b <- u1 + u2;

    # compute probabilitie for u1 values using normal distribution
    moy <- mean(u1b)
    std <- sqrt((n-1)/n)*sd(u1b)
    f.value <- pnorm(u1b, mean = moy, sd = std)

    # perform the Student one sample test
    em <- 0.5;
    p.value <- tprobaCalc(moyV, stdV, n, m2-1, em);

    # decomposition parameters
    comp <- sqrt(smat.eig$values);
    comp.w <- comp / sum(comp);
    comp.wcum <- cumsum(comp.w);

    list(idnames=idnames, u1=u1b, FC2=FC2, f.value=f.value, p.value=p.value,
                          comp=comp, comp.w=comp.w, comp.wcum=comp.wcum);
}
