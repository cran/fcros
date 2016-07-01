\name{cghInfo}
\alias{cghInfo}

\docType{data}

\title{Example of aCGH "info file" for the fcros package.}

\description{This dataset is part of the Agilent 244K design array
probes description file which is used in the publicatin of Sircoulomb et al. 
2010, BMC Cancer, 10:539.
For our illustration, we used only 3 chromosomes (7, 8 and 9) data. The 
complete data are available in the file "GPL9158-3352.txt", see the Gene 
Expression Omnibus website and accession number GSE17907.}

\usage{data(cghInfo)}
\format{
  A data frame with 33,3613 rows for 7 columns with information on probes.
  \describe{
    \item{\code{Index}:}{a numeric used for the probe}
    \item{\code{ProbeName}:}{a text for the unique probe ID}
    \item{\code{GeneSymbol}:}{a text with the gene symbol associated with the probe}
    \item{\code{Chromosome}:}{a text with the chromosome index associated with the probe}
    \item{\code{Start}:}{a start position value for the sequence associated with the probe}
    \item{\code{End}:}{an end position value for the sequence associated with the probe}
    \item{\code{Cytoband}:}{a text for the cytoband associated with the probe}
  }
}

\details{These information data are obtained from the Agilent 244K design array
probes description file, see reference below for details.}

\references{F Sircoulomb, I Bekhouche, P Finetti, J Adelaide, AB Hamida,
J Bonansea, S Raynaud, C Innocenti, E Charafe-Jauffret, C Tarpin, FB Ayed,
P Viens, J Jacquemier, F Bertucci, D Birnbaum and M Chaffanet; Genome 
profiling of ERBB2-amplifief breast cancers. BMC Cancer, 2010, 10:539}

\examples{
   data(cghInfo)

   summary(cghInfo)
}
\keyword{datasets}
