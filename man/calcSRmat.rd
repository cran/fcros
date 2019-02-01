\name{calcSRmat}
\alias{calcSRmat}

\title{Calculation of the sorted rank matrix from the dataset}

\description{This is an internal function used to calculate the sorted
rank matrix. It is used in the functions: fcros() and pfco().}

\usage{calcSRmat(xdata, cont, test, log2.opt=0, trim.opt=0.25)}

\arguments{
  \item{xdata}{ A matrix or a table containing two biological conditions
                dataset to process for detecting differentially expressed
                genes: \code{xdata}.}
  \item{cont}{ A vector containing the label names of the control samples:
               \code{cont} = c("cont01", "cont02", ...).}
  \item{test}{ A vector containing the label names of the test samples:
               \code{test} = c("test01", "test02", "test03", ...).}
  \item{log2.opt}{ A scalar equals to 0 or 1. The value 0 (default) means that
              data in the matrix "xdata" are expressed in a log2 scale:
              \code{log2.opt} = 0}
  \item{trim.opt}{ A scalar between 0 and 0.5. The value 0.25 (default) means
              that 25\% of the lower and the upper rank values of each gene are not
              used for computing its statistics "ri", i.e. the interquartile range
              rank values are averaged: \code{trim.opt} = 0.25}
}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D and Kaster P, Fold change rank ordering statistics:
                    a new method for detecting differentially expressed 
                    genes, BMC bioinformatics, 2014, 15:14}

\examples{
#    data(fdata);
}
