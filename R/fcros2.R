fcros2 <-
function(xdata1, xdata2, cont, test, log2.opt = 0, trim.opt = 0.25) {
    n1 <- nrow(xdata1);
    n2 <- nrow(xdata2);
    if (n1 != n2) stop('xdata1 and xdata2 should have the same number of rows');

    n <- n1;
    # compute FC matrix from dataset 1
    fc1 <- fcrosFCmat(xdata1, cont, test, log2.opt, trim.opt);

    # compute FC matrix from dataset 2
    fc2 <- fcrosFCmat(xdata2, cont, test, log2.opt, trim.opt);

    m1 <- ncol(fc1$fcMat);
    m2 <- ncol(fc2$fcMat);
    m <- m1+m2;
    fc <- matrix(c(fc1$fcMat, fc2$fcMat), ncol = m);

    # compute the fold changes
    FC <-  0.5*(fc1$FC + fc2$FC);
    FC2 <- 0.5*(fc1$FC2+ fc2$FC2);

    # perform analysis with fold changes matrix
    af <- fcrosMod(fc, log2.opt, trim.opt);
    ri <- af$ri;
    p.value <- af$p.value;
    f.value <- af$f.value;
    bounds <- af$bounds;
    params <- af$params;
    params_t <- af$params_t;

    list(FC=FC, FC2=FC2, ri=ri, p.value=p.value, f.value=f.value, bounds=bounds, params=params, params_t=params_t);
}
