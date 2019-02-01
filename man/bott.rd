\name{bott}
\alias{bott}

\docType{data}

\title{Example of sequencing data to test the rankReads function.}

\description{This is a subset of data taken from the Bottomly dataset
see http://bowtie-bio.sourceforge.net/recount/. The complete dataset has
36,536 rows or ENSEMBL identifiers (genes) and 21 columns (samples). For the
"bott" data, the first 5,000 rows, first 3 and last 3 samples were used.}

\usage{data(bott)}
\format{
  A data frame with 5,000 rows and 7 columns.
  \describe{
    \item{\code{gene}:}{ENSEMBL ID}
    \item{\code{SRX033480}:}{sequencing values for the first B6 mouse}
    \item{\code{SRX033488}:}{sequencing values for the second B6 mouse}
    \item{\code{SRX033481}:}{sequencing values for the third B6 mouse}
    \item{\code{SRX033493}:}{sequencing values for the first D2 mouse}
    \item{\code{SRX033486}:}{sequencing values for the second D2 mouse}
    \item{\code{SRX033494}:}{sequencing values for the third B6 mouse}
  }
}

\details{"bott" is a subset of a complet dataset obtained using a single 
RNA-seq reads from C57BL/6J (B6) and DBA/J2 (D2) mice.}

\references{Bottomly et al. Evaluating Gene Expression in C57BL/6J and DBA/J2
Mouse Striatum Using RNA-seq and Microarrays, PLoS One, 6(3)e17820, 2011}

\examples{
   data(bott)
   summary(bott)
}
\keyword{datasets}
