\name{fcrosRMat}
\alias{fcrosRMat}

\title{fcrosRMat: This function is used by fcros2()}

\description{This function allows to compute a matrix of
pairwise fold changes from a two biological conditions data in a single
table or matrix matrix.}

\usage{fcrosRMat(xdata, cont, test, log2.opt=0, trim.opt=0.3)}

\arguments{
  \item{xdata}{ A matrix or table containing a two biological conditions dataset to
                process for detecting differentially expressed genes: \code{xdata}}
  \item{cont}{ A vector containing the names of the control samples:
       \code{cont} = c("cont01","cont02",...)}
  \item{test}{ A vector containing the names of the test samples:
       \code{test} = c("test01", "test02", "test03", ...)}
  \item{log2.opt}{ A scalar equals to 0 or 1. Value 0 (default) means that
              data in "xdata" are expressed in a log2 scale: \code{log2.opt} = 0}
  \item{trim.opt}{ A scalar between 0 and 0.5. Value 0.3 (default) means that
              30% of lower and upper ranks of gene are not used for computing "ri",
              i.e. the average rank, see below: \code{trim.opt} = 0.3}
}

\details{Names appearing in the parameters "cont" and "test" should match
some names of the columns in the matrix or table "xdata".}

\value{ This function returns a data frame containing 3 components
    \item{rmat }{A matrix of fold change ranks associated to genes in the data 
               matrix "xdata". These rank will be used to obtain f-value for genes.}
    \item{FC }{Fold changes for genes in the data matrix "xdata". These fold changes are calculated as
             a ratio of averages from the test and control samples.}
    \item{FC2 }{Fold changes for genes in the data matrix "xdata". These fold changes are calculated as
              a robust average of the fold changes obtained using several pairs of samples.}
}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D and Kastner P, Fold change rank ordering statistics: a new method for detecting
                    differentially expressed genes, BMC Bioinformatics 2014, 15:14}

\examples{
   data(fdata);

   cont <- c("cont01", "cont07", "cont03", "cont04", "cont08");
   test <- c("test01", "test02", "test08", "test09", "test05");
   log2.opt <- 0;

   # perform fcrosRMat()
   rmat <- fcrosRMat(fdata, cont, test, log2.opt);

   # plot histogram of the log2 fold change in dataset
   hist(log2(rmat$FC),nclass=100);
}
