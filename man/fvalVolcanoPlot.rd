\name{fvalVolcanoPlot}
\alias{fvalVolcanoPlot}

\title{Performs a volcano plot of the FCROS/PFCO statistics}

\description{This function allows to have a volcano like plot using the output 
results of the function fcros(), fcros2(), fcrosMod(), pfco() or pfcoMod():
p-values versus robust fold changes (FC2). 
The p-value are transformed using -log10(), while FC2 are transformed using log2().}

\usage{fvalVolcanoPlot(af, thr = 0.05)}

\arguments{
  \item{af}{ This is an object obtained using the functions fcros(),
        fcros2(), fcrosMod(), pfco() or pfcoMod():
        \code{af = fcros(xdata, cont, test)}}
  \item{thr}{ The threshold to obtain the DE genes in the dataset (red plots):
              \code{thr = 0.05}}
}

\value{ This function displays on the screen a volcano like plot using the f-values
        and the robust fold changes (FC2).}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D and Kastner P, Fold change rank ordering statistics: 
            a new method for detecting differentially expressed genes,
            BMC bioinformatics, 2014, 15:14\cr

            Dembele D and Kastner P, Comment on: Fold change rank ordering statistics:
                    a new method for detecting differentially expressed
                    genes, BMC Bioinformatics, 2016, 17:462}

\examples{
   data(fdata);
   rownames(fdata) <- fdata[,1];

   cont <- c("cont01", "cont07", "cont03", "cont04", "cont08");
   test <- c("test01", "test02", "test08", "test09", "test05");
   log2.opt <- 0;

   # perform fcros()
   af <- fcros(fdata, cont, test, log2.opt);
   
   # Volcano plot
   fvalVolcanoPlot(af, thr = 0.01);
}
