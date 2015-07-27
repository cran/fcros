fcrosChrSummary <- function(af, chrInfo, chromosomes = c(1:22,"X","Y"), alpha = 0.1) {
   xinfo <- chrInfo
   n <- nrow(xinfo)
   xinfo$Position <- 0.5*(as.numeric(xinfo$Start) + as.numeric(xinfo$End))
   xinfo$f.L2R <- log2(af$FC2)
   xinfo$f.value <- af$f.value
   a1 <- 0.5*alpha
   a2 <- 1-0.5*alpha
   seg <- c(rep(0, n))
   seg[af$f.value <= a1] <- -1
   seg[af$f.value >= a2] <- 1
   xinfo$f.call <- seg
   chrSumm <- matrix(c(rep(0,3*length(chromosomes))), ncol = 3)
   colnames(chrSumm) <- c("Chr", "nLoss", "nGain")
   chr <- paste("chr", chromosomes[1], sep = "")
   idx <- which(xinfo$Chromosome == chr)
   xda <- xinfo[idx,]
   xdr <- order(xda$Position)
   xinfo2 <- xda[xdr,]
   chrSumm[1,1] <- chromosomes[1]
   chrSumm[1,2] <- sum(xda$f.call == -1)
   chrSumm[1,3] <- sum(xda$f.call == 1)
   for (i in 2:length(chromosomes)) {
       chr <- paste("chr", chromosomes[i], sep = "")
       idx <- which(xinfo$Chromosome == chr)
       xda <- xinfo[idx,]
       xdr <- order(xda$Position)
       xinfo2 <- rbind(xinfo2, xda[xdr,])
       chrSumm[i,1] <- chromosomes[i]
       chrSumm[i,2] <- sum(xda$f.call == -1)
       chrSumm[i,3] <- sum(xda$f.call == 1)
   }
   list(chrData = xinfo2, chrSumm = chrSumm)
}
