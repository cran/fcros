\name{voomReads}
\alias{voomReads}

\title{Transformation of read count values}

\description{This function allows to transform the count values
associated with Sequencing reads. The purpose of this transformation
is to allow applying normal-based microarray-like statistical methods
to RNA-seq read counts. Log2 values are returned.}

\usage{voomReads(x, Rm=1e+06)}

\arguments{
  \item{x}{ This is a read counts matrix}
  \item{Rm}{ A constant used in the transformation}
}

\value{ This function returns a data matrix of the same size as input 
        matrix x. The values of this matrix are expressed in log2 scale.
}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Charity W Law, Yunshun Chen, Wei Shi and Gordon K Smyth,
            voom: precision weights unlock linear model analysis tools
            for RNA-seq read counts,Genome Biology, 2014, 15R29}

\examples{
#   data(fdata);
}
