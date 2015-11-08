\name{fcrosFCmat}
\alias{fcrosFCmat}

\title{fcrosFCmat(): Computation of a matrix of change values}

\description{This function allows to compute a matrix of
pairwise fold changes from a two biological conditions dataset and
is used internally by fcros2()}

\usage{fcrosFCmat(xdata, cont, test, log2.opt=0, trim.opt=0.25)}

\arguments{
  \item{xdata}{ A matrix or a table containing a two biological conditions
                dataset to process for obtaining a matrix of fold changes
                that results from pairwise comparisons of samples:
                \code{xdata}. The first column of the matrix "xdata"
                should contain the gene IDs or symbols.}
  \item{cont}{ A vector containing label names of the control samples:
       \code{cont} = c("cont01", "cont02", ...)}
  \item{test}{ A vector containing label names of the test samples:
       \code{test} = c("test01", "test02", "test03", ...)}
  \item{log2.opt}{ A scalar equals to 0 or 1. The value 0 (default) means that
              the data values in matrix "xdata" are expressed in a log2 scale:
              \code{log2.opt} = 0}
  \item{trim.opt}{ A scalar between 0 and 0.5. The value 0.25 (default) means
                   that 25\% of the lower and the upper rank values of each gene 
                   are not used for computing the statistic "ri", i.e. the 
                   inter-quartile rank values are averaged: 
                   \code{trim.opt} = 0.25}
}

\details{The label names appearing in the parameters "cont" and "test" should match
some label names of the columns in the matrix "xdata". The matrix xdata can
contain other label names not used.}

\value{ This function returns a data frame containing 3 components
    \item{idnames}{ A vector containing the list of IDs or symbols associated with genes}
    \item{fcMat }{A matrix of fold changes associated with genes. This matrix is
                obtained using pairwise comparisons of samples
                in the dataset "xdata".}
    \item{FC }{The fold changes for genes in the matrix "xdata". These
             fold changes are calculated as a ratio of averages from the
             test and the control samples. Non log scale values are used
             in this calculation.}
    \item{FC2 }{The robust fold changes for genes in the matrix "xdata". These
              fold changes are calculated as a trimmed mean of the fold
              changes obtained using pairs of samples. Non log scale values
              are used in this calculation.}
}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D and Kastner P, Fold change rank ordering statistics: 
                    a new method for detecting differentially expressed 
                    genes, BMC Bioinformatics, 2014, 15:14}

\examples{
   data(fdata);

   cont <- c("cont01", "cont07", "cont03", "cont04", "cont08");
   test <- c("test01", "test02", "test08", "test09", "test05");
   log2.opt <- 0;

   # perform fcrosRMat()
   fc <- fcrosFCmat(fdata, cont, test, log2.opt);

   # plot histogram of the fold change (log2 scale) in dataset
   hist(log2(fc$FC),nclass = 50);
}
