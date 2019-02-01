rankReads <- function(xdata, cont, test, meth=0, Ttimes=10, err=0.1,
                             trim.opt=0.25, rseed=60) {
    n <- nrow(xdata);
    idnames <- rownames(xdata);
    xcol <- colnames(xdata);
    n.xcol <- length(xcol);
    idx1 <- xcol %in% cont;
    m1 <- sum(idx1);
    idx2 <- xcol %in% test;
    m2 <- sum(idx2);
    m <- m1+m2;

    # adjust uniform value level if necessary
    if ((err < 0) || (err > 1)) err <- 0.1;

    # form data matrix
    fmat <- matrix(c(rep(0,n*m)), ncol = m);
    x1 <- xdata[, idx1];
    fmat[,1:m1] <- as.matrix(x1);
    x2 <- xdata[, idx2];
    fmat[,(m1+1):m] <- as.matrix(x2);
    colnames(fmat) <- c(cont, test);
    rownames(fmat) <- rownames(xdata)

    set.seed(rseed);
    if (Ttimes > 0) { # perform Ttimes runs
       rmat <- matrix(0, nrow=n, ncol=Ttimes)
       for (i in 1:Ttimes) {
           fmat2 <- fmat + matrix(c(runif(n*m, 0, err)), ncol=m);
           if (meth == 0) {
              rmat.tmp <- (apply(fmat2, 2, rank, ties.method = "average"))/n;
              rmat[,i] <- apply(rmat.tmp, 1, mean, trim.opt)
           }
           else {
              tmp <- log2(fmat2)
              af <- pfco(tmp, cont, test, trim.opt=trim.opt);
              rmat[,i] <- af$ri;
           }
       }
       moy <- apply(rmat, 1, mean);
       std <- apply(rmat, 1, sd);
       if (meth) {
          stat <- 12*Ttimes*std^2;
          pval <- pf(stat, Ttimes-1, Ttimes);
          moyT <- mean(moy)
          stdT <- sqrt((n-1)/n)*sd(moy)
          f.value <- pnorm(moy, mean = moyT, sd = stdT)
          FC <- af$FC; FC2 <- af$FC2; p.value <- af$p.value; ri <- af$ri;

          list(idnames=idnames, FC=FC, FC2=FC2, ri=ri, f.value=f.value,
                                p.value=p.value, score=pval);
       } else {
          score <- std/moy;

          list(idnames=idnames, moy=moy, score=score);
       }
    } else { # allow to perform one run of the fcros method only
       fmat2 <- log2(fmat + matrix(c(runif(n*m, 0, err)), ncol=m));
       af <- pfco(fmat2, cont, test, trim.opt=trim.opt);

       list(idnames=af$idnames, FC=af$FC, FC2=af$FC2, ri=af$ri, 
            f.value=af$f.value, p.value=af$p.value, bounds=af$bounds, 
            params=af$params, params_t=af$params_t);
    }
}
