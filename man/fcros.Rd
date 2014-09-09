\name{fcros}
\alias{fcros}

\title{FCROS for detecting differentially expressed genes}

\description{Implementation of a method based on fold change ranks ordering 
             statistics for detecting differentially expressed genes in a
             dataset in a single table.
             This function should be used with two biological conditions
             dataset (microarray or RNA-seq, ...).
             Using several combinations of pairs of samples from the
             two biological conditions, fold changes (FC) are calculated.
             For each combination, the FC obtained are sorted in increasing
             order and corresponding ranks are associated to genes. Then, a
             statistic is associated to the robust average ordered ranks for
             each gene.}

\usage{fcros(xdata, cont, test, log2.opt = 0, trim.opt = 0.25)}

\arguments{
  \item{xdata}{ A matrix or table containing two biological conditions 
                dataset to process for detecting differentially expressed
                genes: \code{xdata}}
  \item{cont}{ A vector containing the names of the control samples:
               \code{cont} = c("cont01", "cont02", ...)}
  \item{test}{ A vector containing the names of the test samples:
               \code{test} = c("test01", "test02", "test03", ...)}
  \item{log2.opt}{ A scalar equals to 0 or 1. Value 0 (default) means that
              data in "xdata" are expressed in a log2 scale:
              \code{log2.opt} = 0}
  \item{trim.opt}{ A scalar between 0 and 0.5. Value 0.25 (default) means 
              that 25\% of lower and upper ranks of gene are not used for
              computing "ri", i.e. inter quartile range ranks are averaged,
              (see output value of this function below): 
              \code{trim.opt} = 0.25}
}

\details{Names appearing in the parameters "cont" and "test" should match
some names of the columns in the data matrix "xdata".}

\value{ This function returns a data frame containing 8 components
    \item{ri }{The average of ranks associated to genes.
             These values are rank statistics leading to f-values
             and p-values.}
    \item{FC }{Fold changes for genes in the dataset. These fold changes are
             calculated as a ratio of averages from the test and the control
             samples.}
    \item{FC2 }{Robust fold changes for genes. These fold changes are
             calculated as a robust average of the fold changes or ratios 
             obtained from dataset samples.}
    \item{f.value }{F-values are probabilities associated to genes using 
             the "mean" and the
             "standard deviation" ("sd") of values in "ri". The "mean" and
             "sd" are used as normal distribution parameters.}
    \item{p.value }{P-values associated to genes. These values are obtained
             from the f-values.}
    \item{bounds }{Two values which are the lower and the upper bounds or 
             the minimum and the maximum values of "ri".}
    \item{params }{Three values which are estimates for parameters "delta"
            (average difference between consecutive ordered average of rank)
            "mean" (mean value of "ri") and the standard deviation ("sd") 
            of "ri".}
    \item{params_t }{Three values which are theoretical levels for 
                   parameters "delta", "mean" and "sd".}
}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D and Kastner P, Fold change rank ordering statistics: 
                    a new method for detecting differentially expressed 
                    genes, BMC Bioinformatics 2014, 15:14}

\examples{
   data(fdata);

   cont <- c("cont01", "cont07", "cont03", "cont04", "cont08");
   test <- c("test01", "test02", "test08", "test09", "test05");
   log2.opt <- 0;
   trim.opt <- 0.25;

   # perform fcros
   af <- fcros(fdata, cont, test, log2.opt, trim.opt);

   # now select top 20 down and/or up regulated genes
   top20 <- fcrosTopN(af$f.value,20);
   alpha1 <- top20$alpha[1];
   alpha2 <- top20$alpha[2];
   id.down  <- matrix(0,1);
   id.up <- matrix(0,1);
   n <- length(af$FC);
   f.value <- af$f.value;

   idown <- 1;
   iup <- 1;
   for (i in 1:n) {
       if (f.value[i] <= alpha1) { id.down[idown] <- i; idown <- idown+1; }
       if (f.value[i] >= alpha2) { id.up[iup] <- i; iup <- iup+1; }
   }

   data.down <- fdata[id.down[1:(idown-1)],];
   ndown <- nrow(data.down);
   data.up <- fdata[id.up[1:(iup-1)],];
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
