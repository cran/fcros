\name{rmatTrim}
\alias{rmatTrim}

\title{Calculation of the reduced matrix containing sorted rank values}

\description{This is an internal function for using a C code for
calculating the reduced matrix of sorted rank values.}

\usage{rmatTrim(rvect, n, m, idx, m2, rvect2)}

\arguments{
  \item{rvect}{ Vector containing the full matrix with sorted rank values}
  \item{n}{ Number of genes or probes in the dataset}
  \item{m}{ Number of columns of the full matrix of sorted rank values}
  \item{idx}{ Indexes of the columns to keep in the reduced matrix}
  \item{m2}{ Number of columns of the reduced matrix}
  \item{rvect2}{ Vector containing reduced matrix with sorted rank values}
}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D and al, submitted}

\examples{
#    data(fdata);
}
