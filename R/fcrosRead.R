fcrosRead <-
function(file) {
    xdata <- read.table(file, sep="\t", header=TRUE);
    list(xdata = xdata);
}
