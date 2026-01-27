\name{varBeta}
\alias{varBeta}

\title{Compute variance of a beta distribution from data}

\description{This is an internal function. 
Given a vector with values between 0 and 1 and assumed to
come from a beta distribution, this function return the variance of 
the distribution. A trim parameter allow to have a robust value}

\usage{varBeta(x, trim.opt)}

\arguments{
  \item{x}{ vector with components between 0 and 1}
  \item{trim.opt}{ value between 0 and 0.5 used as trim parameter}
}

\value{ This function returns one value equalt to the variance of the best
        beta distribution of entry x
}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D, Analysis of high-throughput biological data using their
                    rank values, Stat Meth Med Res, 2019, 28(8)2276-2291}

\examples{
#    data(fdata);
}
