\name{rankReads}
\alias{rankReads}

\title{Use of rank ordering statistics for count reads data}

\description{Implementation of a method based on fold change rank ordering
             statistics for detecting differentially expressed genes in a
             count reads data. This function is used with two biological 
             conditions RNA-seq dataset.
             Using pairwise combinations of samples from the
             two biological conditions, fold changes (FC) are calculated.
             For each combination, the FC obtained are sorted in increasing
             order and corresponding rank values are associated with genes.
             Then, a statistic is assigned to the robust average ordered rank
             values for each gene/probe. Uniform small values are added to reads
             to overcome ties. A score is computed to localyze genes with almost
             zero count reads for all samples.}

\usage{rankReads(xdata, cont, test, nrun=10, err=0.01, trim.opt=0.25,
                        rseed=57)}

\arguments{
  \item{xdata}{ A matrix or a table containing two biological conditions
                dataset (count reads) to process for detecting differentially 
                expressed genes: \code{xdata}. The first column of the matrix "xdata"
                should contain the gene IDs or their names.}
  \item{cont}{ A vector containing the label names of the control samples:
               \code{cont} = c("cont01", "cont02", ...).}
  \item{test}{ A vector containing the label names of the test samples:
               \code{test} = c("test01", "test02", "test03", ...).}
  \item{nrun}{ An integer used to associate a score with each gene, that allows
               to localize genes with almost zero count value for all samples.
               The value 20 (default) means that the dataset is used 20 times
               and small uniform values are added at each time:
              \code{nrun} = 20}
  \item{err}{ This is the amount of the uniform values to add to count values.
               The value 0.01 (default) is used:
              \code{err} = 0.01}
  \item{trim.opt}{ A scalar between 0 and 0.5. The value 0.25 (default) means
              that 25\% of the lower and the upper rank values of each gene are not
              used for computing its statistics "ri", i.e. the inter-quartile range
              rank values are averaged: \code{trim.opt} = 0.25}
  \item{rseed}{ This value allow to set the computer random generation seed
               value in order to be able to have the same results for runs
               performed at different times: \code{rseed} = 57}
}

\details{Label names appearing in the parameters "cont" and "test" should match
with some label names in the columns of the data matrix "xdata". It is not 
necessary to use all label names appearing in the columns of the dataset matrix.}

\value{ This function returns a data frame containing 10 components
    \item{idnames}{ A vector containing the list of IDs or symbols associated with genes}
    \item{ri }{The average of rank values associated with genes.
             These values are rank values statistics leading to f-values
             and p-values.}
    \item{FC }{The fold changes for genes in the dataset. These fold changes are
             calculated as a ratio of averages from the test and the control
             samples. Non log scale values are used in the calculation.}
    \item{FC2 }{The robust fold changes for genes. These fold changes are
             calculated as a trimmed mean of the fold changes or ratios
             obtained from the dataset samples.  Non log scale values are used 
             in the calculation.}
    \item{f.value }{The f-values are probabilities associated with genes using
             the "mean" and the "standard deviation" ("sd") of the statistics "ri".
             The "mean" and "sd" are used as a normal distribution parameters.}
    \item{p.value }{The p-values associated with genes. These values are obtained
             from the fold change rank values and one sample t-test.}
    \item{bounds }{Two values, which are the lower and the upper bounds or
             the minimum and the maximum values of the non standardized "ri".}
    \item{params }{Three values, which are the estimates for the parameters "delta"
            (average difference between consecutive ordered average of rank values)
            "mean" (mean value of "ri") and the standard deviation ("sd")
            of "ri".}
    \item{params_t }{Three values which are theoretical levels for 
                   parameters "delta", "mean" and "sd".}
    \item{score }{probability of the Fisher test on variances. Higher values
            are associated with genes with almost zero count reads for all
            samples.}
}


\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D, Manuscript in Preparation}

\examples{
#    data(fdata);
}
