\name{calcSRmatMod}
\alias{calcSRmatMod}

\title{Calculation of the sorted rank matrix from the dataset}

\description{This is an internal function used to calculate the sorted
rank matrix. It is used in the functions: fcrosMod() and pfcoMod().}

\usage{calcSRmatMod(xdata, samp, log2.opt=0, trim.opt=0.25)}

\arguments{
  \item{xdata}{ A matrix containing fold changes or ratios from a biological
                dataset to process for searching differentially expressed
                genes or for detecting recurrent copy number aberrations 
                regions: \code{fcMat}.}
  \item{samp}{ A vector of sample label names which should appear in the columns
               of the matrix fcMat: \code{samp}.}
  \item{log2.opt}{ A scalar equals to 0 or 1. The value 0 (default) means that
                values in the matrix "fcMat" are expressed in a log2 scale:
                \code{log2.opt} = 0}
  \item{trim.opt}{ A scalar between 0 and 0.5. The value 0.25 (default) means
                that 25\% of the lower and the upper rank values for each gene 
                are not used for computing the statistic "ri", i.e. the 
                interquartile range rank values are averaged:
                \code{trim.opt} = 0.25}
}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D and Kaster P, Fold change rank ordering statistics:
                    a new method for detecting differentially expressed 
                    genes, BMC bioinformatics, 2014, 15:14}

\examples{
#    data(fdata);
}
