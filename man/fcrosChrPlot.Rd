\name{fcrosChrPlot}
\alias{fcrosChrPlot}

\title{Plot a chromosome data}

\description{This function generates a picture using the chromosome data.}

\usage{fcrosChrPlot(chrData, thr = 0.05, deb = 100, fin = 1e7)}

\arguments{
  \item{chrData}{ A chromosome data obtained from an output of the function
             fcrosChrSummary(): \code{xinfo2} = chrSummary(af, xinfo, chromosomes, alpha)\cr
                                \code{idx} = which(xinfo2$chrData$Chromosome == "chr1")\cr
                                \code{chrData} = xinfo2$chrData[idx, ]}
  \item{thr}{ The probability threshold leading to the selection of the significant probes:
              \code{thr} = 0.05}
  \item{deb}{ This parameter allows to specify the start position of the chromosome region
              for plotting. It can be used for zooming. Negative value will lead to the plot
              of all chromosome data. \code{deb} = 100}
  \item{fin}{ This parameter allows to specify the end position of the chromosome region for
              plot. It can be used for zooming. Negative value will lead to the plot of all
              chromosome data. \code{thr} = 1e7}
}

\value{ This function generates a picture on the screen}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D, Detection of recurrent chromosomal copy number
                    aberrations from high-throughput biological data,
                    Manuscript in preparation}

\examples{
    # load CGH data and info files
    data(cghData)
    data(cghInfo)
    noms = colnames(cghData)
    m = length(noms)
    samp  <- noms[2:m]

    # associate statistics with probes
    af <- fcrosMod(cghData, samp, log2.opt = 0, trim.opt = 0.25)

    chromosomes = c(7:9)
    alpha = 0.05
    
    # summarize results for each chromosome
    xinfo2 = fcrosChrSummary(af, cghInfo, chromosomes, alpha)

    # focused on chromosome 7 data
    idx = which(xinfo2$chrData$Chromosome == "chr7")
    chrData = xinfo2$chrData[idx, ]

    # Plot chromosome 7 data
    fcrosChrPlot(chrData, thr = alpha)
}
