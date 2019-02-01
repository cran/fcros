pfco <- function(xdata, cont, test, log2.opt=0, trim.opt=0.25) {
   n <- nrow(xdata);
   idnames <- rownames(xdata);
   xcol <- colnames(xdata);
   n.xcol <- length(xcol);
   idx1 <- xcol %in% cont;
   m1 <- sum(idx1);
   idx2 <- xcol %in% test;
   m2 <- sum(idx2);

   # compute matrix of sorted FC ranks
   rankmat <- calcSRmat(xdata, cont, test, log2.opt, trim.opt);
   rmat.sr <- rankmat$rmat.sr;
   FC <- rankmat$FC;
   FC2 <- rankmat$FC2;
   moyV <- rankmat$moyV;
   stdV <- rankmat$stdV;
   m2c <- rankmat$m2c;

   # compute symmetric matrix from rank values and its eigen values
   smat <- (t(rmat.sr) %*% rmat.sr)/n;
   smat.eig <- eigen(smat);
   v <- smat.eig$vectors[,1:2];
   if (v[1,1] < 0) v <- -v;
   u1 <- mean(v[,1])*rmat.sr %*% v[,1];
   u2 <- mean(v[,2])*rmat.sr %*% v[,2];
   u1b <- u1 +  u2;

   # compute probabilities for u1 components using normal distribution
   moy <- mean(u1b)
   std <- sqrt((n-1)/n)*sd(u1b)
   f.value <- pnorm(u1b, mean = moy, sd = std)

   # perform the Student one sample test
   em <- 0.5;
   p.value <- tprobaCalc(moyV, stdV, n, m2c-1, em);

   # decomposition parameters
   comp <- sqrt(smat.eig$values);
   comp.w <- comp / sum(comp);
   comp.wcum <- cumsum(comp.w);

   list(idnames=idnames, FC=FC, FC2=FC2, ri=u1b, f.value=f.value,
    p.value=p.value, comp=comp, comp.w=comp.w, comp.wcum=comp.wcum);
}
