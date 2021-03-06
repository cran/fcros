\name{moyStdCalc}
\alias{moyStdCalc}

\title{Calculation of the mean and the standard deviation of the full
                   or reduced matrix of sorted ranks}

\description{This is an internal function for using a C code in the
calculation of the mean and the standard deviation of the full or reduced 
matrix with sorted rank values. The calculations are performed for each row.}

\usage{moyStdCalc(rvect, n, m)}

\arguments{
  \item{rvect}{ Vector containing the full or reduced matrix with sorted rank values}
  \item{n}{ Number of genes or probes in the dataset}
  \item{m}{ Number of columns of the full or reduced matrix of sorted rank values}
}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D, Analysis of high biological data using their rank
values, Stat Methods Med Res, accepted for publication, 2018}

\examples{
#    data(fdata);
}
