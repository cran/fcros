\name{scoreThr}
\alias{scoreThr}

\title{Determine numerically a threshold for the ranking score}

\description{Given a log 10 transformed of the sorted ranking score values
(dscore), this function determines numerically a threshold (inflection
point) for significance of sequencing value level. User should specify the 
lower and the upper bounds of the slope region containing the inflection point.}

\usage{scoreThr(dscore, deb, fin)}

\arguments{
  \item{dscore}{ transformed score values obtained using the function
                 rankReads()}
  \item{deb}{ integer for lower bound of the slope region containing the
              inflection point. User should perfom a plot to determine this
              value}
  \item{fin}{ integer for upper bound of the slope region containing the
              inflection point. User should perfom a plot to determine this
              value}
}

\value{ This function returns two values: ``pos'' and ``thr'' which are the
        index or the number of genes with significant values and the score
        (-log10(score)) threshold
}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D, Analysis of high-throughput biological data using their 
                    rank values, Stat Meth Med Res, 2019, 28(8)2276-2291}

\examples{
#    data(fdata);
}
