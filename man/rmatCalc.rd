\name{rmatCalc}
\alias{rmatCalc}

\title{Calculation of the FC matrix using pairwise comparisons}

\description{This is an internal function for using a C code to calculate a
vector form of the matrix of fold changes using pairwise comparison of data samples.}

\usage{rmatCalc(fvect, n, m1, m2)}

\arguments{
  \item{fvect}{ Two biological conditions dataset matrix}
  \item{n}{ Number of genes or probes in the dataset}
  \item{m1}{ Number of samples in the first biological condition}
  \item{m2}{ Number of samples in the second biological condition}
}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D and Kaster P, Fold change rank ordering statistics:
                    a new method for detecting differentially expressed 
                    genes, BMC bioinformatics, 2014, 15:14\cr

            Dembele D and Kastner P, Comment on: Fold change rank ordering statistics:
                    a new method for detecting differentially expressed
                    genes, BMC Bioinformatics, 2016, 17:462}

\examples{
#    data(fdata);
}
