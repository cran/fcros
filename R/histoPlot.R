histoPlot <- function(af, nbins = 50) {
    if (length(af$ri) > 0) {
       hist(af$ri, nclass = nbins, xlab = "", main = "FCROS statistics", xlim = c(0, 1));
    }
    else {
       hist(af$u1, nclass = nbins, xlab = "", main = "PFCO statistics");
    }
}
