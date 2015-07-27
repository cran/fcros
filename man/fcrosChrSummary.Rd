\name{fcrosChrSummary}
\alias{fcrosChrSummary}

\title{Summarization of the detection results for a list of chromosomes}

\description{ From an outpout object of the function fcrosMod(), the chromosomes
information object, the list of chromosomes and a threshold, this function 
creates two objects containing ordered chromosome data and summary results.}

\usage{fcrosChrSummary(af, chrInfo, chromosomes = c(1:22,"X","Y"), alpha = 0.1)}

\arguments{
  \item{af}{ An output object of the function fcrosMod():
             \code{af} = fcrosMod(xdata, samp, log2.opt, trim.opt)}
  \item{chrInfo}{ A data frame containing chromosomes information (probe name,
                  gene symbol, chromosome index, start position, end position
                  and the cytoband). These information should appear with the labels
                  ProbeName, GeneSymbol, Chromosome, Start, End, Cytoband. Additional
                  information may be used. Only labels Chromosome, Start and End are 
                  mandotory.}
  \item{chromosomes}{ A list of chromosomes. Default setting is a list with all
                      chromosomes: \code{chromosomes} = (1:22,"X","Y")}
  \item{alpha}{ A threshold allowing to select significant probes based on probabilities.
                Default setting is 0.1 (10\% of error) \code{thr} = 0.1}
}

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

    # associate statistics with probes in the dataset
    af <- fcrosMod(cghData, samp, log2.opt = 0, trim.opt = 0.25)

    chromosomes = c(7:9)
    alpha = 0.05
    
    # summarize results for each chromosome
    xinfo2 = fcrosChrSummary(af, cghInfo, chromosomes, alpha)

    # display the number of significant probes for each chromosome
    xinfo2$chr$Summ
}
