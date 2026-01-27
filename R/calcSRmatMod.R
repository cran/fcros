calcSRmatMod <- function(xdata, samp, log2.opt=0, trim.opt=0.25) {
   n <- nrow(xdata);
   xcol <- colnames(xdata);
   n.xcol <- length(xcol);
   idx <- xcol %in% samp;
   m <- sum(idx);

   # form data matrix
   fmat <- matrix(c(rep(0,n*m)), ncol=m);
   if (log2.opt) {
      x1 <- log2(xdata[,idx]);
   } else {
      x1 <- xdata[,idx];
   }
   fmat <- as.matrix(x1);

   # compute matrix containing pairwise fold changes
   rvectC <- c(fmat[, 1:m]);
   rmat <- apply(fmat, 2, rank, ties.method="average")/n;

   # if (trim.opt), reduce the standard rank matrix
   if ((trim.opt > 0) & (trim.opt < 0.5)) {
      deb <- round(trim.opt * m) + 1;
      fin <- m - deb + 1;
      idx <- deb:fin;
      m2c <- length(idx);
      rvect <- c(rmat[,1:m]);
      rvect2 <- rmatTrim(rvect, n, m, idx, m2c);
      rmat.sr <- matrix(rvect2, ncol = m2c);
      rvect <- rvect2;
      rmat.val <- moyStdCalc(rvect, n, m2c);
      moyV <- rmat.val$moyC;
      stdV <- rmat.val$stdC;
      FC2 <- fc2Calc(rvectC, n, m, idx, m2c);
      rm(rvect);
      rm(rmat.val);
      rm(rvectC);
   }
   else {
      rmat.sr <- rmat;
      m2c <- m;
      idx <- 1:m2c;
      rvect <- c(rmat.sr[,1:m2c]);
      rmat.val <- moyStdCalc(rvect, n, m2c);
      moyV <- rmat.val$moyC;
      stdV <- rmat.val$stdC;
      FC2 <- fc2Calc(rvectC, n, m, idx, m2c);
      rm(rvect);
      rm(rmat.val);
      rm(rvectC);
   }
   list(rmat.sr=rmat.sr, moyV=moyV, stdV=stdV, FC2=FC2, m2c=m2c);
}
