pfcoMod <- function(fcMat, samp, log2.opt=0, trim.opt=0.25) {
    n <- nrow(fcMat);
    idnames <- rownames(fcMat);

    # compute matrix of sorted FC ranks
    fmod <- calcSRmatMod(fcMat, samp, log2.opt, trim.opt);

    rmat.sr <- fmod$rmat.sr;
    moyV <- fmod$moyV;
    stdV <- fmod$stdV;
    FC2 <- fmod$FC2;
    m2 <- fmod$m2;

    # compute the symmetric matrix with sorted rank values and its eigen values
    smat <- (t(rmat.sr) %*% rmat.sr)/n;
    smat.eig <- eigen(smat);
    v <- smat.eig$vectors[,1:2];
    if (v[1,1] < 0) v <- -v;
    u1 <- mean(v[,1])* rmat.sr %*% v[,1];
    u2 <- mean(v[,2])* rmat.sr %*% v[,2];
    u1b <- u1 + u2;

    # compute probabilitie for u1 values using normal distribution
    moy <- mean(u1b)
    std <- sqrt((n-1)/n)*sd(u1b)
    f.value <- pnorm(u1b, mean = moy, sd = std)

    # perform the Student one sample test
    em <- 0.5;
    p.value <- tprobaCalc(moyV, stdV, n, m2-1, em);

    # decomposition parameters
    comp <- sqrt(smat.eig$values);
    comp.w <- comp / sum(comp);
    comp.wcum <- cumsum(comp.w);

    list(idnames=idnames, ri=u1b, FC2=FC2, f.value=f.value, p.value=p.value,
                          comp=comp, comp.w=comp.w, comp.wcum=comp.wcum);
}
