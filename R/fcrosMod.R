fcrosMod <- function(fcMat, log2.opt = 0, trim.opt = 0.25) {
    # compute matrix of sorted FC ranks
    n <- nrow(fcMat);
    rmat <- apply(fcMat, 2, rank, ties.method = "average")/n

    # compute vectors of f-values and p-values
    ri <- apply(rmat, 1, mean, trim = trim.opt);
    ris <- sort(ri);

    # compute parameters
    lb <- n*ris[1];
    ub <- n*ris[n];
    delta <- (n-1)*mean(ris[-1]-ris[-n]);

    moy <- mean(ri)
    std <- sd(ri)
    f.value <- pnorm(ri, moy, std)
    p.value <- 2*f.value;
    idx <- (p.value > 1);
    p.value[idx] <- (2.0-p.value[idx]);
    moy_t <- (lb+ub)/(2*n);
    delta_t <- (ub-lb)/(n-1);
    std_t <- delta_t/sqrt(12);
    bounds <- c(lb,ub);
    params <- c(delta,moy,std);
    params_t <- c(delta_t,moy_t,std_t);
    if (log2.opt) { FC2 <- apply(fcMat, 1, mean, trim = trim.opt); }
    else {FC2 <- apply(2^fcMat, 1, mean, trim = trim.opt);}

    list(FC2=FC2, ri=ri, p.value=p.value, f.value=f.value, bounds=bounds, params=params, params_t=params_t)
}
