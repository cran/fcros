\name{chrMerge}
\alias{chrMerge}

\title{Using a C code for merging chromosome segments}

\description{This is an internal function for using a C code while merging
chromosome segments in the segmentation step.}

\usage{chrMerge(nbSeg, idStart, idEnd, lBound, uBound, segVal, segProba,
                         fcall, L2R, nd, dm, sigma)}

\arguments{
  \item{nbSeg}{ Number of current segments}
  \item{idStart}{ Position indexes of the first probes for segments}
  \item{idEnd}{ Positions indexes of the last probes for segments}
  \item{lBound}{ Lower bound position for segments}
  \item{uBound}{ Upper position for segments}
  \item{segVal}{ Change values associated with segments}
  \item{segProba}{ Probabilities associated with segments}
  \item{fcall}{ Detection status associated with probes}
  \item{L2R}{ Change values associated with probes}
  \item{nd}{ Number of acceptable non-detection between two significant of a segment}
  \item{dm}{ Average distance between two consecutive probes of the chromosome}
  \item{sigma}{ Standard deviation of the residual observations, see reference}
}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D, Use of the Perron-Frobenius theorem in the analysis
of high throughput biological data, Manuscript submitted}

\examples{
#    data(fdata);
}
