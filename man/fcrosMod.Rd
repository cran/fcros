\name{fcrosMod}
\alias{fcrosMod}

\title{Search for differentially expressed genes or to detect
                  recurrent copy number aberration probes}

\description{Implementation of a method based on fold change rank ordering
             statistics to search for differentially expressed genes or to 
             detect chromosomal recurrent copy number aberration probes. This
             function should be used with a matrix of fold changes or ratios
             from biological dataset (microarray, RNA-seq, ...). The function
             fcrosMod() is an extention of the function fcros() to a
             dataset which does not contain replicate samples or to a
             dataset with one biological condition dataset. Statistics
             are associated with genes/probes to characterize their change levels.}

\usage{fcrosMod(fcMat, samp, log2.opt = 0, trim.opt = 0.25)}

\arguments{
  \item{fcMat}{ A matrix containing fold changes or ratios from a biological
                dataset to process for searching differentially expressed
                genes or for detecting recurrent copy number aberrations 
                regions: \code{fcMat}. The first column of the matrix "fcMat" 
                should contain the gene IDs or symbols.}
  \item{samp}{ A vector of sample label names which should appear in the columns
               of the matrix fcMat: \code{samp}.}
  \item{log2.opt}{ A scalar equals to 0 or 1. The value 0 (default) means that
                values in the matrix "fcMat" are expressed in a log2 scale:
                \code{log2.opt} = 0}
  \item{trim.opt}{ A scalar between 0 and 0.5. The value 0.25 (default) means
                that 25\% of the lower and the upper rank values for each gene 
                are not used for computing the statistic "ri", i.e. the 
                inter-quartile range rank values are averaged:
                \code{trim.opt} = 0.25}
}

\details{The label names appearing in the parameter "samp" should
match some label names of the columns in the data matrix "xdata". It is not 
necessary to use all label names appearing in the columns of the dataset matrix.}

\value{ This function returns a data frame containing 8 components
    \item{idnames}{ A vector containing the list of IDs or symbols associated with genes}
    \item{ri }{The average of ordered rank values associated with genes in the
             dataset. These values are rank statistics leading to the f-values
             and the p-values.}
    \item{FC2 }{The robust fold changes for genes in matrix "fcMat". These fold
             changes are calculated as a trimed mean of the values in
             "fcMat". Non log scale values are used in this calculation.}
    \item{f.value }{The f-values are probabilities associated with genes using
             the "mean" and the "standard deviation" ("sd") of values in "ri". 
             The "mean" and "sd" are used as a normal distribution parameters.}
    \item{p.value }{The p-values associated with genes. The p-values are obtained
             using a one sample t-test on the fold change rank values.}
    \item{bounds }{Two values, which are the lower and the upper bounds or
             the minimum and the maximum values of the non standardized "ri".}
    \item{params }{Three values, which are the estimates for the parameters "delta"
             (average difference between consecutive ordered average of 
             rank) "mean" (mean value of the "ri") and the standard deviation
             ("sd") of the "ri".}
    \item{params_t }{Three values, which are theoretical levels for
             the parameters "delta", "mean" and "sd".}
}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D, Analysis of high biological data using their rank
values, Stat Methods Med Res, 2018}

\examples{
   data(fdata);

   cont <- c("cont01", "cont07", "cont03", "cont04", "cont08");
   test <- c("test01", "test02", "test08", "test09", "test05");
   log2.opt <- 0;
   trim.opt <- 0.25;

   # perform fcrosMod()
   fc <- fcrosFCmat(fdata, cont, test, log2.opt, trim.opt);
   m <- ncol(fc$fcMat)
   samp <- paste("Col",as.character(1:m), sep = "");
   fc.val <- cbind(fdata[,1], data.frame(fc$fcMat))
   colnames(fc.val) <- c("ID", samp)

   af <- fcrosMod(fc.val, samp, log2.opt, trim.opt);

   # now select top 20 down and/or up regulated genes
   top20 <- fcrosTopN(af, 20);
   alpha1 <- top20$alpha[1];
   alpha2 <- top20$alpha[2];
   id.down  <- matrix(0,1);
   id.up <- matrix(0,1);
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
   plot(t, data.down[1, 2:21], type = "l", col = "blue", xlim = c(1,20),
        ylim = c(0,18), main = "Top down-regulated genes");
   for (i in 2:ndown) {
       lines(t,data.down[i, 2:21], type = "l", col = "blue")
   }

   # now plot down and up regulated genes
   plot(t, data.up[1,2:21], type = "l", col = "red", xlim = c(1,20), 
       ylim = c(0,18), main = "Top up-regulated genes");
   for (i in 2:nup) {
       lines(t, data.up[i,2:21], type = "l", col = "red")
   }
   par(op)
}
