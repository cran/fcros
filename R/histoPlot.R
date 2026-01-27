histoPlot <- function(af, nbins=50) {
   hist(af$ri, nclass=nbins, xlab="", main="FCROS/PFCO statistics", xlim=c(0, 1));
}
