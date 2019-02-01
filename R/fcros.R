fcros <- function(xdata, cont, test, log2.opt=0, trim.opt=0.25) {
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

   # compute averages ranks
   ri <- apply(rmat.sr, 1, mean);
   ris <- sort(ri);

   # compute parameters
   lb <- n*ris[1];
   ub <- n*ris[n];
   delta <- (n-1)*mean(ris[-1]-ris[-n]);

   # compute f-value
   moy <- mean(ri);
   std <- sqrt((n-1)/n)*sd(ri);
   f.value <- pnorm(ri, mean = moy, sd = std);

   # perform the Student one sample test to get p-values
   em <- 0.5;
   p.value <- tprobaCalc(moyV, stdV, n, m2c-1, em);

   moy_t <- (lb+ub)/(2*n);
   delta_t <- (ub-lb)/(n-1);
   std_t <- delta_t/sqrt(12);
   bounds <- c(lb,ub);
   params <- c(delta, moy, std);
   params_t <- c(delta_t,moy_t,std_t);

   list(idnames=idnames, FC=FC, FC2=FC2, ri=ri, f.value=f.value,
    p.value=p.value, bounds=bounds, params=params, params_t=params_t);
}
