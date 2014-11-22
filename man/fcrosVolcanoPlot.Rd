\name{fcrosVolcanoPlot}
\alias{fcrosVolcanoPlot}

\title{fcrosVolcanoPlot(): Performs a volcano plot of the fcros statistics}

\description{This function allows to have a volcano like plot using the output 
results the functions fcros(), fcros2() or fcrosMod(): p-values versus robust
fold changes (FC2). The p-value are transformed using -log10(), while FC2 are
transformed using log2().}

\usage{fcrosVolcanoPlot(af, thr = 0.05)}

\arguments{
  \item{af}{ This is an object obtained using the functions fcros(),
        fcros2() or fcrosMod(): \code{af}}
  \item{thr}{ The threshold to obtain the DE genes in the dataset (red plots):
              \code{thr = 0.05}}
}

\value{ This function displays on the screen a volcano like plot using the p-values
        and the robust fold changes (FC2).}

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
   fcrosVolcanoPlot(af, thr = 0.01);

}
