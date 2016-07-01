\name{rmatCalc}
\alias{rmatCalc}

\title{Calculation the matrix of FC using pairwise comparisons}

\description{This is an internal function for using a C code to
calculate the matrix of fold changes using pairwise comparison of data samples.}

\usage{rmatCalc(fvect, n, m1, m2, rvect, FC)}

\arguments{
  \item{fvect}{ Two biological conditions dataset matrix}
  \item{n}{ Number of genes or probes in the dataset}
  \item{m1}{ Number of samples in the first biological condition}
  \item{m2}{ Number of samples in the second biological condition}
  \item{rvect}{ Vector with the matrix elements that result from the
                pairwise comparison of samples}
  \item{FC}{ Fold Change values associated with genes or probes}
}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D and Kaster P, Fold change rank ordering statistics:
                    a new method for detecting differentially expressed 
                    genes, BMC bioinformatics, 2014, 15:14}

\examples{
#    data(fdata);
}
