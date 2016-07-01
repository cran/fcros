\name{fc2Calc}
\alias{fc2Calc}

\title{Calculation of fold change using pairwise comparison values}

\description{This is an internal function for using a C code to
calculate fold changes using pairwise comparison of samples.}

\usage{fc2Calc(rvect, n, m, idx, m2, fc2)}

\arguments{
  \item{rvect}{ Vector containing the full or reduced matrix with the pairwise
                comparison of samples results}
  \item{n}{ Number of genes or probes in the dataset}
  \item{m}{ Number of columns of the full or reduced matrix of pairwise 
            comparison of samples results}
  \item{idx}{ Indexes of the columns to keep in the pairwise comparison
              of samples}
  \item{m2}{ Number of columns in the full or reduced matrix of 
             comparison of samples}
  \item{fc2}{ Vector containing fold change values for rows}
}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D and al, Manuscript submitted}

\examples{
#    data(fdata);
}
