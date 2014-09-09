\name{fcrosMod}
\alias{fcrosMod}

\title{FCROS for detecting differentially expressed genes in a table with 
             ratios or fold changes}

\description{Implementation of a method based on fold change ranks ordering 
             statistics for detecting differentially expressed genes. This
             function should be used with a matrix of fold changes or ratios
             from biological dataset (microarray, RNA-seq, ...). This 
             function is an extention of the basic fcros() function to a
             dataset which does not contain replicate samples or to a
             dataset with one biological condition dataset.
             This allows its application to experiment over time or to array 
             Comparative Genomic Hybridization (aCGH) datasets. Statistics
             are associated to genes to characterize their change level.}

\usage{fcrosMod(fcMat, log2.opt = 0, trim.opt = 0.25)}

\arguments{
  \item{fcMat}{ A matrix containing fold changes or ratios from a biological
                dataset to process for detecting differentially expressed 
                genes: \code{fcMat}}
  \item{log2.opt}{ A scalar equals to 0 or 1. Value 0 (default) means that
                values in "fcMat" are expressed in a log2 scale:
                \code{log2.opt} = 0}
  \item{trim.opt}{ A scalar between 0 and 0.5. Value 0.25 (default) means
                that 25\% of lower and upper ranks of gene are not used for
                computing "ri", i.e. inter quartile range ranks are averaged
                (see output value of this function below):
                \code{trim.opt} = 0.25}
}

\details{Names appearing in the parameters "cont" and "test" should match
some names of the columns in the data matrix "xdata".}

\value{ This function returns a data frame containing 7 components
    \item{ri }{The average of ordered ranks associated to genes in the 
             dataset. These values are rank statistics leading to f-values 
             and p-values.}
    \item{FC2 }{Robust fold changes for genes in "fcMat". These fold
             changes are calculated as a robust average of the values in
             "fcMat".}
    \item{f.value }{F-values are probabilities associated to genes using 
                  the "mean" and the
             "standard deviation" ("sd") of values in "ri". The "mean" and
             "sd" are used as normal distribution parameters.}
    \item{p.value }{P-values associated to genes. The p-values are obtained
             using the f-values.}
    \item{bounds }{Two values which are the lower and the upper bounds or 
             the minimum and the maximum values of "ri".}
    \item{params }{Three values which are estimates for parameter "delta" 
             (average difference between consecutive ordered average of 
             rank) "mean" (mean value of "ri") and the standard deviation 
             ("sd") of "ri".}
    \item{params_t }{Three values which are theoretical levels for 
                   parameters "delta", "mean" and "sd".}
}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D and Kastner P, Fold change rank ordering statistics: 
              a new method for detecting differentially expressed genes, 
              BMC Bioinformatics 2014, 15:14}

\examples{
   data(fdata);

   cont <- c("cont01", "cont07", "cont03", "cont04", "cont08");
   test <- c("test01", "test02", "test08", "test09", "test05");
   log2.opt <- 0;
   trim.opt <- 0.25;

   # perform fcrosMod
   fc <- fcrosFCmat(fdata, cont, test, log2.opt, trim.opt);
   af <- fcrosMod(fc$fcMat, log2.opt, trim.opt);

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
