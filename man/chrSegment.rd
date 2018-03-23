\name{chrSegment}
\alias{chrSegment}

\title{Segmentation of a chromosome data}

\description{This function allows to segment a chromosome data}

\usage{chrSegment(chrData, nd = 10)}

\arguments{
  \item{chrData}{ A chromosome data obtained from an output of the function
             dataSummary(): \code{xinfo2} = dataSummary(af, xinfo, chromosomes, alpha)\cr
                                \code{idx} = which(xinfo2$xinfo.s$Chromosome == "chr1")\cr
                                \code{chrData} = xinfo2$xinfo.s[idx, ]}
  \item{nd}{ The acceptable number of non-detected probes which can separate two
             significant probes in a segment. Default setting value is 10: \code{nd} = 10}
}

\value{ This function returns a data frame containing 6 information for each segment
    \item{idStart}{ The start position indexes associated with segments}
    \item{idEnd}{ The End position indexes associated with segments}
    \item{lBounds}{ The lower bound positions associated with segments}
    \item{uBounds}{ The upper bound positions associated with segments}
    \item{segL2R}{ The change values associated with segments}
    \item{segProba}{ The probabilities associated with segments}
}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D, Analysis of high biological data using their rank
values, Stat Methods Med Res, 2018}

\examples{
    # load CGH data and info files
    data(cghData)
    data(cghInfo)
    noms = colnames(cghData)
    m = length(noms)
    samp  <- noms[2:m]

    # associate statistics with probes in the dataset
    af <- pfcoMod(cghData, samp, log2.opt = 0, trim.opt = 0.25)

    chromosomes = c(7:9)
    alpha = 0.05

    # summarize results for each chromosome
    xinfo2 = dataSummary(af, cghInfo, chromosomes, alpha)

    # focused on chromosome 7 data
    idx = which(xinfo2$xinfo.s$Chromosome == "7")
    chrData = xinfo2$xinfo.s[idx, ]

    # segment chromosome 7 data
    chrSeg = chrSegment(chrData, nd = 15)

    # show first 10 segment results
    chrSeg[1:10,]
}
