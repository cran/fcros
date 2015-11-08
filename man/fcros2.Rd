\name{fcros2}
\alias{fcros2}

\title{fcros2() to search for differentially expressed genes/probes}

\description{Implementation of a method based on fold change or ratio rank
             ordering statistics for detecting differentially expressed 
             genes. This function should be used with dataset in two separate
             tables and from two biological conditions datasets (microarray,
             RNA-seq, ...). Pairwise combinations of samples from
             the two biological conditions are used to obtain a matrix of
             fold changes. For each combination, the FCs obtained are
             sorted in an increasing order and the corresponding rank values are
             associated with genes/probes. 
             Then, a statistic is associated with each gene/probe.}

\usage{fcros2(xdata1, xdata2, cont, test, log2.opt = 0, trim.opt = 0.25)}

\arguments{
    \item{xdata1}{ A matrix or a table containing two biological conditions
                   dataset to process for detecting differentially expressed
                   genes: \code{xdata1}. The first column of the matrix "xdata1" 
                   should contain the gene IDs or symbols.}
    \item{xdata2}{ A matrix or a table containing two biological conditions
                   dataset to process fordetecting differentially expressed
                   genes: \code{xdata2}. The first column of the matrix "xdata2" 
                   should contain the gene IDs or symbols.}
  \item{cont}{ A vector containing the label names of the control samples:
       \code{cont} = c("cont01", "cont02", ...)}
  \item{test}{ A vector containing the label names of the test samples:
       \code{test} = c("test01", "test02", "test03", ...)}
  \item{log2.opt}{ A scalar equals to 0 or 1. The value 0 (default) means that
              data in the tables "xdata1" and "xdata2" are expressed in a log2 scale:
              \code{log2.opt} = 0}
  \item{trim.opt}{ A scalar between 0 and 0.5. The value 0.25 (default) means
                   that 25\% of the lower and  the upper rank values for a gene
                   are not used for computing  the statistics "ri", i.e. the 
                   inter-quartile range
                   rank values are averaged: \code{trim.opt} = 0.25}
}

\details{The label names appearing in the parameters "cont" and "test" should 
match some label names in the columns of the data matrices "xdata1" and "xdata2". 
It is not necessary to use all column label names appearing in matrices "xdata1" 
and "xdata2". However, it is assumed that the same genes (same IDs or symbol) are 
used in xdata1 and xdata2.}

\value{ This function returns a data frame containing 9 components
    \item{idnames}{ A vector containing the list of the IDs or symbols associated with genes}
    \item{ri }{The average of rank values associated with genes in the
             datasets. These values are rank statistics leading to f-values
             and p-values.}
    \item{FC }{The fold changes for genes. These fold changes are
             calculated as a ratio of averages from the test and control
             samples. Non log scale values are used in the calculation.}
    \item{FC2 }{The robust fold changes for genes. These fold changes are
             calculated as a trimmed mean of the fold changes or ratios
             obtained from the dataset samples using the pairwise comparisons.
              Non log scale values are used in the calculation.}
    \item{f.value }{The f-values are probabilities associated with genes. These 
              values are obtained using the "mean" and the
              "standard deviation" ("sd") of the statistics "ri". The "mean" and
              "sd" are used as a normal distribution parameters.}
    \item{p.value }{The p-values associated with genes. These values are obtained
              from the f-values.}
    \item{bounds }{Two values which are the lower and the upper bound or the
              minimum and the maximum values of non standardized "ri".}
    \item{params }{Three values, which are the estimates for the parameters "delta"
             (average difference between consecutive ordered average of rank) 
             "mean" (mean value of variable "ri") and the standard deviation ("sd")
             of variable "ri".}
    \item{params_t }{Three values, which are the theoretical levels for
                   parameters "delta", "mean" and "sd".}
}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D and Kastner P, Fold change rank ordering statistics: 
             a new method for detecting differentially expressed genes, 
             BMC Bioinformatics, 2014, 15:14}

\examples{
   data(fdata);

   cont <- c("cont01", "cont02", "cont03", "cont04", "cont05",
             "cont06", "cont07", "cont08", "cont09", "cont10");
   test <- c("test01", "test02", "test03", "test04", "test05",
             "test06", "test07", "test08", "test09", "test10");
   log2.opt <- 0;
   trim.opt <- 0.25;

   # perform fcros2()
   xdata1 <- fdata[,c(2:5, 12:17)];
   xdata2 <- fdata[,c(6:11, 18:21)];
   af2 <- fcros2(xdata1, xdata2, cont, test, log2.opt, trim.opt);

   # now select top 20 down and/or up regulated genes
   top20 <- fcrosTopN(af2$f.value, 20);
   alpha1 <- top20$alpha[1];
   alpha2 <- top20$alpha[2];
   id.down  <- matrix(0,1);
   id.up <- matrix(0,1);
   n <- length(af2$FC);
   f.value <- af2$f.value;

   idown <- 1;
   iup <- 1;
   for (i in 1:n) {
       if (f.value[i] <= alpha1) { id.down[idown] <- i; idown <- idown + 1; }
       if (f.value[i] >= alpha2) { id.up[iup] <- i; iup <- iup + 1; }
   }

   data.down <- fdata[id.down[1:(idown-1)], ];
   ndown <- nrow(data.down);
   data.up <- fdata[id.up[1:(iup-1)], ];
   nup <- nrow(data.up);

   # now plot down regulated genes
   t <- 1:20;
   op = par(mfrow = c(2,1));
   plot(t, data.down[1,2:21], type = "l", col = "blue", xlim = c(1,20), 
        ylim = c(0,18), main = "Top down-regulated genes");
   for (i in 2:ndown) {
       lines(t,data.down[i,2:21], type = "l", col = "blue")
   }

   # now plot down and up regulated genes
   plot(t, data.up[1,2:21], type = "l", col = "red", xlim = c(1,20), 
       ylim = c(0,18), main = "Top up-regulated genes");
   for (i in 2:nup) {
       lines(t, data.up[i,2:21], type = "l", col = "red")
   }
   par(op)
}
