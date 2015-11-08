\name{fcrosChrMerge}
\alias{fcrosChrMerge}

\title{Using a C code for merging}

\description{This is an internal function for using the C code while merging
chromosome segments in the segmentation step.}

\usage{fcrosChrMerge(nbSeg, idStart, idEnd, lBound, uBound, segVal, segProba,
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

\references{Dembele D and al, Detection of recurrent chromosomal copy number
                    aberrations from high-throughput data,
                    Manuscript submitted}

\examples{
#    data(fdata);
#    cont <- c("cont01", "cont07", "cont03", "cont04", "cont08");
#    test <- c("test01", "test02", "test08", "test09", "test05");
#    log2.opt <- 0;
#    trim.opt <- 0.25;
#    af <- fcros(fdata, cont, test, log2.opt, trim.opt);
#
#    fcrosWrite(af, file = "test2delete_values.txt");
#    fcrosWrite(af, file = "test2delete_params.txt", values = FALSE);
}
