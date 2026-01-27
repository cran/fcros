chrMerge <- function(nbSeg, idStart, idEnd, lBound, uBound, segVal,
                            segProba, fcall, L2R, nd, dm, sigma) {
    .Call(C_merge, nbSeg = as.integer(nbSeg),
                segIdS = as.integer(idStart),
                segIdE = as.integer(idEnd),
                segLB = as.double(lBound),
                segUB = as.double(uBound),
                segVal = as.double(segVal),
                segProba = as.double(segProba),
                fcall = as.integer(fcall),
                L2R = as.double(L2R),
                nd = as.integer(nd),
                dm = as.double(dm),
                sigma = as.double(sigma));
}
