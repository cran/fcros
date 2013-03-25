\name{fcrosRMat}
\alias{fcrosRMat}

\title{fcrosRMat: This function is used by fcros2()}

\description{This function allows to compute a matrix of
pairwise fold changes from a two biological conditions dataset.}

\usage{fcrosRMat(xdata,cont,test,log2.opt)}

\arguments{
  \item{xdata}{ Two biological conditions dataset to process for detecting
                differentially expressed genes: \code{xdata}}
  \item{cont}{ Vector containing the names of the control samples:
       \code{cont} = c("cont01","cont02",...)}
  \item{test}{ Vector containing the names of the test samples:
       \code{test} = c("test01", "test02", "test03", ...)}
  \item{log2.opt}{ A scalar equals to 0 or 1. Value 0 (default) means that
              data in "xdata" are expressed in a log2 scale: \code{log2.opt} = 0}
}

\details{Names appearing in the parameters "cont" and "test" should match
some names of the columns in the "xdata" datatset}

\value{ This function return a data frame containg 3 components
    \item{rmat }{The matrix of fold change ranks associated to genes in the dataset.
              These rank will be used to obtain f-value for genes}
    \item{FC }{Fold change of genes in the dataset. It is calculated as a ratio of
             averages from the test and control samples.}
    \item{FC2 }{Fold change of genes in the dataset. It is calculated as an average
              of the fold changes obtained using several pairs of samples.}
}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele et al. (2013), Fold change ranks ordering statistics for detecting
                    differentially expressed genes, Manuscript in preparation}

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
