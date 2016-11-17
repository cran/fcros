\name{tprobaCalc}
\alias{tprobaCalc}

\title{Calculation of the Student one sample test probabilities}

\description{This is an internal function for using a C code for
performing a Student one sample test for each row of the full or 
reduced matrix with sorted rank values.}

\usage{tprobaCalc(moy, std, n, dl, em, proba)}

\arguments{
  \item{moy}{ Vector containing average of rank values for rows}
  \item{std}{ Vector containing standard deviation of rank values for rows}
  \item{n}{ Number of genes or probes in the dataset}
  \item{dl}{ Degree of freedom in the test. This is equal to the number of 
             the columns in the full or reduced matrix of sorted rank values 
             minus one}
  \item{em}{ Expected average rank values for each row}
  \item{proba}{ Vector containing Student t-test probability for rows}
}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D, Use of the Perron-Frobenius theorem in the analysis
of high throughput biological data, Manuscript submitted}

\examples{
#    data(fdata);
}
