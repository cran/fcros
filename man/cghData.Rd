\name{cghData}
\alias{cghData}

\docType{data}

\title{Example of aCGH "data file" for the fcros package.}

\description{This dataset is part of the study performed and published 
by Sircoulomb et al. 2010, BMC Cancer, 10:539.
For our illustration, we used the first 10 patients'
data and for only 3 chromosomes: 7, 8 and 9. The complete data are
available from the Gene Expression Omnibus website under the accession 
number GSE17907.}

\usage{data(cghData)}
\format{
  A data frame with 33,3613 observations for 10 samples.
  \describe{
    \item{\code{Probes}:}{a text for the unique probe ID}
    \item{\code{GSM447252}:}{a numeric vector with log2 change values}
    \item{\code{GSM447253}:}{a numeric vector with log2 change values}
    \item{\code{GSM447254}:}{a numeric vector with log2 change values}
    \item{\code{GSM447255}:}{a numeric vector with log2 change values}
    \item{\code{GSM447256}:}{a numeric vector with log2 change values}
    \item{\code{GSM447257}:}{a numeric vector with log2 change values}
    \item{\code{GSM447258}:}{a numeric vector with log2 change values}
    \item{\code{GSM447259}:}{a numeric vector with log2 change values}
    \item{\code{GSM447260}:}{a numeric vector with log2 change values}
    \item{\code{GSM447261}:}{a numeric vector with log2 change values}
  }
}

\details{This dataset is obtained using array Comparative Genomic
Hybridization technology. Agilent 244K design arrays have been used
to monitor breast cancer patients, see reference below for details.}

\references{F Sircoulomb, I Bekhouche, P Finetti, J Adelaide, AB Hamida, 
J Bonansea, S Raynaud, C Innocenti, E Charafe-Jauffret, C Tarpin, FB Ayed,
P Viens, J Jacquemier, F Bertucci, D Birnbaum and M Chaffanet; Genome 
profiling of ERBB2-amplifief breast cancers. BMC Cancer, 2010, 10:539}

\examples{
   data(cghData)

   summary(cghData)

   #  perform boxplot of data sample values
   boxplot(cghData[,2:11])
}
\keyword{datasets}
