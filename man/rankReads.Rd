\name{rankReads}
\alias{rankReads}

\title{This function computes a score to assess the significance of sequencing values}

\description{Implementation of two methods based (1) on the coefficient of variation
             or (2) on the fold change rank ordering statistics for detecting genes
             with significant sequencing values (gwssv). A score is obtained for each
             gene and a threshold allows to select the number of gwssv.}

\usage{rankReads(xdata, cont, test, meth=0, Ttimes=10, err=0.1, trim.opt=0,
                        rseed=60)}

\arguments{
  \item{xdata}{ A matrix or a table containing sequencing dataset.
                The rownames of xdata is used for the output idnames.}
  \item{cont}{ A vector containing the label names of the control samples:
               \code{cont} = c("cont01", "cont02", ...).}
  \item{test}{ A vector containing the label names of the test samples:
               \code{test} = c("test01", "test02", "test03", ...).}
  \item{meth}{ This parameter allows to specify the approach to use. The value
               0 (defaul) means the coefficient of variation is used. When
               non-zero value is given, the fcros method used: \code{meth} = 0}
  \item{Ttimes}{ The number of perturbed data to use.
               The value 10 (default) means that the dataset is used 20 times
               and small uniform values are added at each time:
              \code{Ttimes} = 10}
  \item{err}{ This is the amount of the uniform values to add to count values.
               The value 0.1 (default) is used:
              \code{err} = 0.1}
  \item{trim.opt}{ A scalar between 0 and 0.5. The value 0.25 (default) means
              that 25\% of the lower and the upper rank values of each gene are not
              used for computing its statistics "ri", i.e. the inter-quartile range
              rank values are averaged: \code{trim.opt} = 0.25}
  \item{rseed}{ This value allow to set the computer random generation seed
               value in order to be able to have the same results for runs
               performed at different times: \code{rseed} = 58}
}

\details{Label names appearing in the parameters "cont" and "test" should match
with some label names in the columns of the data matrix "xdata". It is not
necessary to use all label names appearing in the columns of the dataset matrix.
For a general purpose dataset, one of these parametere can be empty.}

\value{ This function returns a data frame containing 10 components when meth=1
        and 3 components when meth=0
    \item{idnames}{ A vector containing the list of IDs or symbols associated with genes}
    \item{score }{coefficient of variation (meth=0) or Fisher-Snedecor test p-value
            (meth=1). Smaller (higher) values are associated with genes with significant
            (non significant) sequencing values.}
    \item{moy }{trimmed means associated with genes (when meth=0).}
    \item{ri }{The average of rank values associated with genes when meth=1.
             These values are rank values statistics leading to f-values
             and p-values (when meth=1).}
    \item{FC }{The fold changes for genes in the dataset. These fold changes are
             calculated as a ratio of averages from the test and the control
             samples. Non log scale values are used in the calculation (when meth=1).}
    \item{FC2 }{The robust fold changes for genes. These fold changes are
             calculated as a trimmed mean of the fold changes or ratios
             obtained from the dataset samples.  Non log scale values are used
             in the calculation (when meth=1).}
    \item{f.value }{The f-values are probabilities associated with genes using
             the "mean" and the "standard deviation" ("sd") of the statistics "ri".
             The "mean" and "sd" are used as a normal distribution parameters 
             (when meth=1).}
    \item{p.value }{The p-values associated with genes. These values are obtained
             from the fold change rank values and one sample t-test (when meth=1).}
}


\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D, Analysis of high-throughput biological data using their
                    rank values, Stat Meth Med Res, 2019, 28(8)2276-2291}

\examples{
   data(bott);
   cont <- c("SRX033480", "SRX033488", "SRX033481");
   test <- c("SRX033493", "SRX033486", "SRX033494");
   n <- nrow(bott);

   x2 <- tcnReads(bott[,c(cont,test)])
   idx.ok <- (apply(x2, 1, sum) != 0)
   xdata <- x2[,c(cont,test)]
   rownames(xdata) <- bott[,1]
   idx.ok <- (apply(x2, 1, sum) != 0)
   tt2 <- sum(idx.ok)

   raf10.cv <- rankReads(xdata, cont, test, meth=0)
   raf10.pv <- rankReads(xdata, cont, test, meth=1)
   score.cv <- -log10(sort(raf10.cv$score))
   score.pv <- -log10(sort(raf10.pv$score))
   tmp <- scoreThr(score.cv, 2500, 3500)
   tmp

   tmp <- scoreThr(score.pv, 2500, 3500)
   tmp

   op <- par(mfrow = c(1,2))
   plot(score.cv, xlab = "index of genes",
      ylab = "-log10(sorted(score)", main = "rs.cv", type = "l",
      col = "blue", panel.first = grid())
   plot(score.pv, xlab = "index of genes",
      ylab = "-log10(sorted(score)", main = "rs.pv", type = "l",
      col = "blue", panel.first = grid())
   par(op)
}
