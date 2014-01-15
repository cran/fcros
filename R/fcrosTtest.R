fcrosTtest <-
function(xdata, cont, test, log2.opt = 0) {
    n <- nrow(xdata);

    xcol <- colnames(xdata);
    n.xcol <- length(xcol);
    idx1 <- xcol %in% cont;
    m1 <- sum(idx1==TRUE);
    idx2 <- xcol %in% test;
    m2 <- sum(idx2==TRUE);
    m <- m1+m2;
    m1m2 <- m1*m2;

    # form ttest data matrix
    tmat <- matrix(c(rep(0,n*m)),ncol=m);
    if (log2.opt) {
       x1 <- log2(xdata[,idx1==TRUE]);
    } else {
       x1 <- xdata[,idx1==TRUE];
    }
    tmat[,1:m1] <- as.matrix(x1);
    if (log2.opt) {
       x2 <- log2(xdata[,idx2==TRUE]);
    } else {
       x2 <- xdata[,idx2==TRUE];
    }
    tmat[,(m1+1):m] <- as.matrix(x2);

    FC <- matrix(c(rep(0,n)),ncol=1);
    p.value <- matrix(c(rep(0,n)),ncol=1);
    stat <- matrix(c(rep(0,n)),ncol=1);
    # compute fold changes and p-values
    for (i in 1:n) {
        x1 <- tmat[i,1:m1];
        x2 <- tmat[i,(m1+1):m];
        tt <- t.test(x1,x2,var.equal=TRUE);
        p.value[i] <- tt$p.value;
        stat[i] <- tt$statistic;
        FC[i] <- mean(2^x2)/mean(2^x1);
    }
     list(FC = FC, p.value = p.value, stat = stat);
}
