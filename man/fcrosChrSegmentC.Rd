\name{fcrosChrSegmentC}
\alias{fcrosChrSegmentC}

\title{Segmentation of a chromosome data}

\description{This function allows to segment a chromosome data}

\usage{fcrosChrSegmentC(chrData, nd = 10)}

\arguments{
  \item{chrData}{ A chromosome data obtained from an output of the function
             fcrosChrSummary(): \code{xinfo2} = chrSummary(af, xinfo, chromosomes, alpha)\cr
                                \code{idx} = which(xinfo2$chrData$Chromosome == "chr1")\cr
                                \code{chrData} = xinfo2$chrData[idx, ]}
  \item{nd}{ The acceptable number of non-detected probes which can separate two
             significant probes in a segment. Default setting value is 5: \code{nd} = 5}
}

\value{ This function returns a data frame containing 6 information for each segment
    \item{idStart}{ The start position index of the segment}
    \item{idEnd}{ The End position index of the segment}
    \item{lBound}{ The lower bound position index of the segment}
    \item{uBound}{ The upper bound position of the segment}
    \item{segL2R}{ The change value associated with the segment}
    \item{segProba}{ The probability associated with the segment}
}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D, A method for detection of recurrent chromosomal copy
                    number aberrations from high-throughput biological data,
                    Manuscript submitted}

\examples{
#    data(fdata);
#    cont <- c("cont01", "cont07", "cont03", "cont04", "cont08");
#    test <- c("test01", "test02", "test08", "test09", "test05");
#    log2.opt <- 0;
#    trim.opt <- 0.25;
#    af <- fcros(fdata, cont, test, log2.opt, trim.opt);

#    fcrosWrite(af, file = "test2delete_values.txt");
#    fcrosWrite(af, file = "test2delete_params.txt", values = FALSE);
}
