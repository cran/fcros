#############################################################################
#DD:20190114
#############################################################################
# dscore <- -log10(sort(res$score)) where res is output of rankReads()
# "deb" and "fin" are bounds of the region containing the inflection point
#############################################################################
scoreThr <- function (dscore, deb, fin) {
    n <- length(dscore)
    dds <- dscore[-n] - dscore[-1]
    idx <- rep(0, n-1)
    dmx <- max(dds[deb:fin])
    idx <- dds == dmx
    tt <- which(idx == TRUE)
    pp <- dscore[tt]
    list(pos=tt, thr=pp)
}
