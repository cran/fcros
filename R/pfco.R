pfco <- function(xdata, cont, test, log2.opt=0, trim.opt=0.25) {
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

    # form data matrix
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
    rvect <- c(rep(0, n*m1m2));
    FC <- matrix(c(rep(0,n)), ncol = 1);
    fvect <- c(fmat[, 1:m])

    rmat.val <- rmatCalc(fvect, n, m1, m2, rvect, FC);
    rmat <- matrix(rmat.val$rvectC, ncol = m1m2);
    FC <- rmat.val$FCC;
    rvectFC <- rmat.val$rvectC;

    # compute the standard ranks matrix
    rmat.s <- (apply(rmat, 2, rank, ties.method = "average"))/n;

    # if (trim.opt), reduce the standard rank matrix
    if ((trim.opt > 0) & (trim.opt < 0.5)) {
       deb <- round(trim.opt * m1m2) + 1;
       fin <- m1m2 - deb + 1;
       idx <- deb:fin;
       m2 <- length(idx);
       rvect <- c(rmat.s[,1:m1m2]);
       rvect2 <- c(rep(0, n*m2));
       rmat.val <- rmatTrim(rvect, n, m1m2, idx, m2, rvect2);
       rmat.sr <- matrix(rmat.val$rvect2, ncol = m2);
       rvect <- rmat.val$rvect2;
       moyV <- c(rep(0,n));
       stdV <- c(rep(0,n));
       rmat.val <- moyStdCalc(rvect, n, m2, moyV, stdV);
       moyV <- rmat.val$moyC;
       stdV <- rmat.val$stdC;
       FC2 <- c(rep(0,n));
       rmat.val <- fc2Calc(rvectFC, n, m1m2, idx, m2, FC2)
       FC2 <- rmat.val$fc2C;
       rm(rvect);
       rm(rmat.val);
       rm(rvectFC);
    }
    else {
         rmat.sr <- rmat.s;
         m2 <- m1m2;
         idx <- 1:m2;
         rvect <- c(rmat.sr[,1:m2]);
         moyV <- c(rep(0,n));
         stdV <- c(rep(0,n));
         rmat.val <- moyStdCalc(rvect, n, m2, moyV, stdV);
         moyV <- rmat.val$moyC;
         stdV <- rmat.val$stdC;
         FC2 <- c(rep(0,n));
         rmat.val <- fc2Calc(rvectFC, n, m1m2, idx, m2, FC2)
         FC2 <- rmat.val$fc2C;
         rm(rvect);
         rm(rmat.val);
         rm(rvectFC);
    }

    # compute symmetric matrix from rank values and its eigen values
    smat <- (t(rmat.sr) %*% rmat.sr)/n;
    smat.eig <- eigen(smat);
    v1 <- smat.eig$vectors[,1];
    if (v1[1] < 0) v1 <- -v1;
    u1 <- rmat.sr %*% v1;

    # compute probabilities for u1 components using normal distribution
    tt <- c(1.0, 0.9, 0.945, 0.97, 0.97, 0.97)
    at <- floor(10*trim.opt)
    moy <- tt[at+1] * mean(u1)
    std <- sd(u1)  
    f.value <- pnorm(u1, mean = moy, sd = std)

    # perform the Student one sample test
    em <- 0.5;
    prob <- c(rep(0,n));
    pval <- tprobaCalc(moyV, stdV, n, m2-1, em, prob);
    p.value <- pval$probaC;

    # decomposition parameters
    comp <- sqrt(smat.eig$values);
    comp.w <- comp / sum(comp);
    comp.wcum <- cumsum(comp.w);

    list(idnames=idnames, FC=FC, FC2=FC2, u1=u1, f.value=f.value,
    p.value=p.value, comp=comp, comp.w=comp.w, comp.wcum=comp.wcum);
}
