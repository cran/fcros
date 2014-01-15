fcros <-
function(xdata,cont,test,log2.opt=0,trim.opt=0.3) {
    n <- nrow(xdata);

    xcol <- colnames(xdata);
    n.xcol <- length(xcol);
    idx1 <- xcol %in% cont;
    m1 <- sum(idx1==TRUE);
    idx2 <- xcol %in% test;
    m2 <- sum(idx2==TRUE);
    m <- m1+m2;
    m1m2 <- m1*m2;

    # form fcros data matrix
    fmat <- matrix(c(rep(0,n*m)),ncol=m);
    x1 <- matrix(c(rep(0,n*m1)),ncol=m1);
    x2 <- matrix(c(rep(0,n*m2)),ncol=m2);
    if (log2.opt) {
         x1 <- log2(xdata[,idx1==TRUE]);
    } else {
         x1 <- xdata[,idx1==TRUE];
    }
    fmat[,1:m1] <- as.matrix(x1);
    if (log2.opt) {
         x2 <- log2(xdata[,idx2==TRUE]);
    } else {
         x2 <- xdata[,idx2==TRUE];
    }
    fmat[,(m1+1):m] <- as.matrix(x2);

    # compute matrix containing pairwise fold changes
    rmat <- matrix(c(rep(0,n*m1m2)),ncol=m1m2);
    k <- 1;
    for (i in 1:m1) {
        for (j in 1:m2) {
            rmat[,k] <- fmat[,m1+j]-fmat[,i];
            k <- k+1;
        }
    }
    # compute the fold changes
    FC <- matrix(c(rep(0,n)),ncol=1);
    for (i in 1:n) {
        x1 <- fmat[i,1:m1];
        x2 <- fmat[i,(m1+1):m];
        FC[i] <- mean(2^x2)/mean(2^x1);
    }
    FC2 = apply(2^rmat[1:n,],1,mean,trim=trim.opt);

    # compute sorted ranks matrix
    rmat.s <- matrix(c(rep(0,n*m1m2)),ncol=m1m2);
    for (k in 1:m1m2) {
        rmat.s[,k] <- rank(rmat[,k],ties.method="first");
    }
    # compute averages ranks
    ri <- matrix(c(rep(0,n)),ncol=1);
    ri <- apply(rmat.s,1,mean,trim=trim.opt)/n;
    ris <- sort(ri);

    # compute parameters
    lb <- n*ris[1];
    ub <- n*ris[n];
    delta <- (n-1)*mean(ris[-1]-ris[-n]);

    # compute p-values
    moy <- mean(ri);
    std <- sd(ri);
    f.value = pnorm(ri,mean=moy,sd=std);
    moy_t <- (lb+ub)/(2*n);
    delta_t <- (ub-lb)/(n-1);
    std_t <- delta_t/sqrt(12);
    bounds <- c(lb,ub);
    params <- c(delta,moy,std);
    params_t <- c(delta_t,moy_t,std_t);

    list(FC=FC, FC2=FC2, ri=ri, f.value=f.value, bounds=bounds, params=params, params_t=params_t);
}
