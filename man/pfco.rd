\name{pfco}
\alias{pfco}

\title{Searching for differentially expressed genes/probes using an
                 approach based on the Perron-Frobenius theorem}

\description{Implementation of a method based on fold change and the Perron
             theorem for detecting differentially expressed genes in a
             dataset.
             This function should be used with two biological conditions
             dataset (microarray or RNA-seq, ...).
             Using pairwise combinations of samples from the
             two biological conditions, fold changes (FC) are calculated.
             For each combination, the FC obtained are sorted in increasing
             order and corresponding rank values are associated with genes.
             Then, a statistic is assigned to the robust average ordered rank
             values for each gene/probe.}

\usage{pfco(xdata, cont, test, log2.opt = 0, trim.opt = 0.25)}

\arguments{
  \item{xdata}{ A matrix or a table containing two biological conditions
                dataset to process for detecting differentially expressed
                genes. The rownames of xdata are used for the output idnames.}
  \item{cont}{ A vector containing the label names of the control samples:
               \code{cont} = c("cont01", "cont02", ...).}
  \item{test}{ A vector containing the label names of the test samples:
               \code{test} = c("test01", "test02", "test03", ...).}
  \item{log2.opt}{ A scalar equals to 0 or 1. The value 0 (default) means that
              data in the matrix "xdata" are expressed in a log2 scale:
              \code{log2.opt} = 0}
  \item{trim.opt}{ A scalar between 0 and 0.5. The value 0.25 (default) means
              that 25\% of the lower and the upper rank values of each gene are not
              used for computing its statistics "ri", i.e. the interquartile range
              rank values are averaged: \code{trim.opt} = 0.25}
}

\details{Label names appearing in the parameter "samp" should match
with some label names in the columns of the data matrix "xdata". It is not 
necessary to use all label names appearing in the columns of the dataset matrix.}

\value{ This function returns a data frame containing 9 components
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
             using a one sample Student t-test on the fold change rank values.}
    \item{comp }{Singular values.}
    \item{comp.w }{Singular values weights.}
    \item{comp.wcum }{Cumulative sum of the singular values weights.}
}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D, Analysis of high biological data using their rank
values, Stat Methods Med Res, accepted for publication, 2018}

\examples{
   data(fdata);
   rownames(fdata) <- fdata[,1];

   cont <- c("cont01", "cont07", "cont03", "cont04", "cont08");
   test <- c("test01", "test02", "test08", "test09", "test05");
   log2.opt <- 0;
   trim.opt <- 0.25;

   # perform pfco()
   af <- pfco(fdata, cont, test, log2.opt, trim.opt);

   # now select top 20 down and/or up regulated genes
   top20 <- fcrosTopN(af, 20);
   alpha1 <- top20$alpha[1];
   alpha2 <- top20$alpha[2];
   id.down  <- matrix(0, 1);
   id.up <- matrix(0, 1);
   n <- length(af$FC);
   f.value <- af$f.value;

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
