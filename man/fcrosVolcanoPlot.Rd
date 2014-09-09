\name{fcrosVolcanoPlot}
\alias{fcrosVolcanoPlot}

\title{fcrosVolcanoPlot: Performs Volcano plot of dataset}

\description{This function allows to have a Volcano like plot of dataset.
It uses the p-values obtained using fcros(), fcros2() or fcrosMod().}

\usage{fcrosVolcanoPlot(af, thr = 0.05)}

\arguments{
  \item{af}{ This is an object obtained using functions fcros(),
        fcros2() or fcrosMod(): \code{af}}
  \item{thr}{ Threshold to get the DE genes in the dataset (red plots):
              \code{thr = 0.05}}
}

\value{ This function performs a volcano plot on the screen using p-values
        and fold changes computed from the dataset.}

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
