\name{tcnReads}
\alias{tcnReads}

\title{Performs a total count normalization of reads}

\description{Given a data table x with count reads, one column for each
sample, this function adjust values in such a way they become comparable
between samples. User can specify the maximum value for total count reads
to be used for each sample.}

\usage{tcnReads(x, maxVal=0)}

\arguments{
  \item{x}{ Data table with count reads, one column for one sample}
  \item{maxVal}{ Target total number reads for each sample. The value 0
                 (default) means that the median total of samples is
                 used: \code{maxVal} = 0}
}

\value{ This function returns a data table of the same size as input
    \item{x2 }{Data table of the same size as x}
}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D, Analysis of high biological data using their rank
                    values, Stat Methods Med Res, accepted for publication, 2018\cr

        Dembele D, Analysis of high-throughput biological data using their rank
                    values, Stat Meth Med Res, 2019, 28(8)2276-2291}

\examples{
#    data(fdata);
}
