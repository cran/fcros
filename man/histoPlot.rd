\name{histoPlot}
\alias{histoPlot}

\title{Histogram plot function of the fcros package results}

\description{This function allows to have a histogram plot. It uses the 
statistics "ri" or "u1" obtained using one of the following functions: fcros(),
fcros2(), fcrosMod(), pfco() or pfcoMod().}

\usage{histoPlot(af, nbins = 50)}

\arguments{
  \item{af}{ This is an object obtained using the function fcros(),
        fcros2(), fcrosMod(), pfco() or pfcoMod()}
  \item{nbins}{ This parameter is used for the number of bins in
                the histogram. Default setting is 50: \code{nbins = 50}}
}

\value{ This function plots a histogram on the screen.}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D and Kastner P, Fold change rank ordering statistics:
            a new method for detecting differentially expressed genes,
            BMC Bioinformatics, 2014, 15:14\cr

            Dembele D and Kastner P, Comment on: Fold change rank ordering statistics:
                    a new method for detecting differentially expressed
                    genes, BMC Bioinformatics, 2016, 17:462}

\examples{
   data(fdata);
   rownames(fdata) <- fdata[,1];

   cont <- c("cont01", "cont07", "cont03", "cont04", "cont08");
   test <- c("test01", "test02", "test08", "test09", "test05");
   log2.opt <- 0;

   # perform fcros() and pfco()
   af <- fcros(fdata, cont, test, log2.opt);
   af2 <- pfco(fdata, cont, test, log2.opt);

   # Histogram plots
   op <- par(mfrow = c(1,2))
      histoPlot(af);
      histoPlot(af2);
   par(op);
}
