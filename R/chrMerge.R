chrMerge <- function(nbSeg, idStart, idEnd, lBound, uBound, segVal,
                            segProba, fcall, L2R, nd, dm, sigma) {
    .Call(C_merge, nbSegC = as.integer(nbSeg),
                segIdSC = as.integer(idStart),
                segIdEC = as.integer(idEnd),
                segLBC = as.double(lBound),
                segUBC = as.double(uBound),
                segValC = as.double(segVal),
                segProbaC = as.double(segProba),
                fcallC = as.integer(fcall),
                L2RC = as.double(L2R),
                ndC = as.integer(nd),
                dmC = as.double(dm),
                sigmaC = as.double(sigma));
}
