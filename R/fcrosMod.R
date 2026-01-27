fcrosMod <- function(fcMat, samp, log2.opt=0, trim.opt=0.25) {
    n <- nrow(fcMat);
    idnames <- rownames(fcMat);

    # compute matrix of sorted FC ranks
    fmod <- calcSRmatMod(fcMat, samp, log2.opt, trim.opt);

    rmat.sr <- fmod$rmat.sr;
    moyV <- fmod$moyV;
    stdV <- fmod$stdV;
    FC2 <- fmod$FC2;
    m2 <- fmod$m2;

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

    list(idnames=idnames, FC2=FC2, ri=ri, p.value=p.value, f.value=f.value,
    bounds=bounds, params=params, params_t=params_t)
}
