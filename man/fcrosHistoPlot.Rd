\name{fcrosHistoPlot}
\alias{fcrosHistoPlot}

\title{fcrosHistoPlot(): Histogram plot function of the fcros results}

\description{This function allows to have a histogram plot. It uses the 
statistics ri obtained using one of the following functions: fcros(), 
fcros2() and fcrosMod().}

\usage{fcrosHistoPlot(af, nbins = 50)}

\arguments{
  \item{af}{ This is an object obtained using the function fcros(),
        fcros2() or fcrosMod(): \code{af = fcrosMod(xdata, samp)}}
  \item{nbins}{ This parameter is used for the number of bins in
                the histogram. Default setting is 50: \code{nbins = 50}}
}

\value{ This function plots a histogram on the screen.}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D and Kastner P, Fold change rank ordering statistics:
            a new method for detecting differentially expressed genes,
            BMC bioinformatics, 2014, 15:14}

\examples{
   data(fdata);

   cont <- c("cont01", "cont07", "cont03", "cont04", "cont08");
   test <- c("test01", "test02", "test08", "test09", "test05");
   log2.opt <- 0;

   # perform fcros
   af <- fcros(fdata, cont, test, log2.opt);
   
   # Volcano plot
   fcrosHistoPlot(af);
}
