\name{fdata}
\alias{fdata}

\docType{data}

\title{Example of dataset for the FCROS method.}

\description{This is a synthetic dataset generated using package "madsim".
It is a two biological conditions data matrix using 10000 genes and 20 samples
(10 control and 10 test samples). 108 and 102 genes of this dataset are down-
and up regulated, respectively.}

\usage{data(fdata)}
\format{
  A data frame with 10000 observations on the following 22 variables.
  \describe{
    \item{\code{index}}{a text for the genes, unique ID}
    \item{\code{cont01}}{a numeric vector with log2 intensities for control sample 01}
    \item{\code{cont02}}{a numeric vector with log2 intensities for control sample 02}
    \item{\code{cont03}}{a numeric vector with log2 intensities for control sample 03}
    \item{\code{cont04}}{a numeric vector with log2 intensities for control sample 04}
    \item{\code{cont05}}{a numeric vector with log2 intensities for control sample 05}
    \item{\code{cont06}}{a numeric vector with log2 intensities for control sample 06}
    \item{\code{cont07}}{a numeric vector with log2 intensities for control sample 07}
    \item{\code{cont08}}{a numeric vector with log2 intensities for control sample 08}
    \item{\code{cont09}}{a numeric vector with log2 intensities for control sample 09}
    \item{\code{cont10}}{a numeric vector with log2 intensities for control sample 10}
    \item{\code{test01}}{a numeric vector with log2 intensities for test sample 01}
    \item{\code{test02}}{a numeric vector with log2 intensities for test sample 02}
    \item{\code{test03}}{a numeric vector with log2 intensities for test sample 03}
    \item{\code{test04}}{a numeric vector with log2 intensities for test sample 04}
    \item{\code{test05}}{a numeric vector with log2 intensities for test sample 05}
    \item{\code{test06}}{a numeric vector with log2 intensities for test sample 06}
    \item{\code{test07}}{a numeric vector with log2 intensities for test sample 07}
    \item{\code{test08}}{a numeric vector with log2 intensities for test sample 08}
    \item{\code{test09}}{a numeric vector with log2 intensities for test sample 09}
    \item{\code{test10}}{a numeric vector with log2 intensities for test sample 10}
    \item{\code{DE_status}}{a numeric vector with values -1, 0 and 1. Value 0 is used for
    no change genes, while -1 and 1 are used for down- and up-regulated genes, respectively.}
  }
}

\details{This dataset is obtained using the microarray data simulation model implemented
in package "madsim". A real microarray data, "madsim_test", was used
as seed and the number of control and test samples were set to 10. Parameter "sdn" was set to
0.3 and all other parameters were set to their default settings.}

\references{Dembele D, A flexible microarray data simulation model. 
Microarrays 2013, v.2, n.2, pp.115-130}

\examples{
   data(fdata)
   
   op = par(mfrow=c(2,1));
   
   #  perform MA plot using samples "cont01" and "cont05"
   A1 = 0.5*(fdata$cont01 + fdata$cont05);
   M1 = fdata$cont05 - fdata$cont01;

   plot(A1, M1, col="red", xlim=c(2,18), ylim=c(-5,5), 
            main="MA plot with two control samples");
   lines(x = c(2,18), y = c(0,0), col = "blue")
   lines(x = c(2,18), y = c(-1,-1), col = "blue")
   lines(x = c(2,18), y = c(1,1), col = "blue")

   #  perform MA plot using samples "cont01" and "test05"
   A2 = 0.5*(fdata$cont01 + fdata$test05);
   M2 = fdata$test05 - fdata$cont01;

   plot(A2, M2, col="red", xlim=c(2,18), ylim=c(-5,5), 
            main="MA plot with one control and one test samples");
   lines(x = c(2,18), y = c(0,0), col = "blue")
   lines(x = c(2,18), y = c(-1,-1), col = "blue")
   lines(x = c(2,18), y = c(1,1), col = "blue")
   par(op)
}
\keyword{datasets}
