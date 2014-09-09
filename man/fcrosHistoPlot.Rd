\name{fcrosHistoPlot}
\alias{fcrosHistoPlot}

\title{fcrosHistoPlot: Histogram plot of the FCROS statistics}

\description{This function allows to have a histogram plot. It uses the 
ri statistics obtained using function fcros(), fcros2() or fcrosMod().}

\usage{fcrosHistoPlot(af)}

\arguments{
  \item{af}{ This is an object obtained using function fcros(),
        fcros2() or fcrosMod(): \code{af}}
}

\value{ This function plots a histogram on the screen.}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D and Kastner P, Fold change rank ordering statistics:
            a new method for detecting differentially expressed genes,
            BMC Bioinformatics 2014, 15:14}

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
