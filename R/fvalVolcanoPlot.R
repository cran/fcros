fvalVolcanoPlot <- function(af, thr=0.05) {
    ratio <- af$FC2;
    pval <- af$f.value;
    idx <- (pval > 0.5);
    pval[idx] <- (1.0 - pval[idx]);
    de.idx <- (pval <= thr);
    x <- log2(ratio);
    y <- -log10(pval);
    plot(x, y, ylab="-log10(p-Value)", xlab="log2(FC)", col="blue", main="Volcano plot");
    lines(x[de.idx], y[de.idx], type="p", col="red");
}
