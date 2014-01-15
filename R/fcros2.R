fcros2 <-
function(xdata1,xdata2,cont,test,log2.opt=0,trim.opt=0.3) {

    n1 <- nrow(xdata1);
    n2 <- nrow(xdata2);
    if (n1 != n2) stop('xdata1 and xdata2 should have the same number of rows');

    n <- n1;
    # compute rank matrix from dataset 1
    r1 <- fcrosRMat(xdata1,cont,test,log2.opt,trim.opt);

    # compute rank matrix from dataset 2
    r2 <- fcrosRMat(xdata2,cont,test,log2.opt,trim.opt);

    m1 <- ncol(r1$rmat);
    m2 <- ncol(r2$rmat);
    m <- m1+m2;
    rmat <- matrix(c(r1$rmat,r2$rmat),ncol=m);

    # compute the fold changes
    FC <- 0.5*(r1$FC+r2$FC);
    FC2 <- 0.5*(r1$FC2+r2$FC2);

    # compute sorted ranks matrix
    rmat.s <- matrix(c(rep(0,n*m)),ncol=m);
    for (k in 1:m) {
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
