\name{dataSummary}
\alias{dataSummary}

\title{Summarization of the detection results for a list of chromosomes}

\description{ From an outpout object of the function fcrosMod() or pfcoMod(),
the chromosomes information object, the list of chromosomes and a threshold, 
this function creates two objects containing ordered chromosome data 
and summary results.}

\usage{dataSummary(af, xinfo, chromosomes = c(1:22,"X","Y"), alpha = 0.05)}

\arguments{
  \item{af}{ An output object of the function fcrosMod() or pfcoMod():\cr
             \code{af} = fcrosMod(xdata, samp, log2.opt, trim.opt)\cr
             \code{af} = pfcoMod(xdata, samp, log2.opt, trim.opt)}
  \item{xinfo}{ A data frame containing chromosomes information (probe name,
                  gene symbol, chromosome index, start position, end position
                  and the cytoband). These information should appear with the labels
                  ProbeName, GeneSymbol, Chromosome, Start, End, Cytoband. Additional
                  information may be used. Only labels Chromosome, Start and End are 
                  mandotory.}
  \item{chromosomes}{ A list of chromosomes. Default setting is a list with all
                      chromosomes:\cr 
                      \code{chromosomes} = (1:22,"X","Y")}
  \item{alpha}{ A threshold allowing to select significant probes based on probabilities.
                Default setting is 0.05 (5\% of error) \code{thr} = 0.05}
}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D, Use of the Perron-Frobenius theorem in the analysis
of high throughput biological data, Manuscript submitted}

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

    # display the number of significant probes for each chromosome
    xinfo2$chrSumm
}
