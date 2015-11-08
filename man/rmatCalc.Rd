\name{rmatCalc}
\alias{rmatCalc}

\title{Using a C code for calculatin the pairwise matrix of FC}

\description{This is an internal function for using the C code while
calculating the pairwise comparison of data samples.}

\usage{rmatCalc(fvect, n, m1, m2, rvect, FC)}

\arguments{
  \item{fvect}{ Two biological conditions data matrix}
  \item{n}{ Number of genes/probes in the dataset}
  \item{m1}{ Number of samples in the first biological condition}
  \item{m2}{ Number of samples in the second biological condition}
  \item{rvect}{ Matrix with pairwise comparison fold change}
  \item{FC}{ Fold Change values associated with genes/probes}
}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D and Kaster P, Fold change rank ordering statistics:
                    a new method for detecting differentially expressed 
                    genes, BMC bioinformatics, 2014, 15:14}

\examples{
#    data(fdata);
}
