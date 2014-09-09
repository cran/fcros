\name{fcros-package}
\alias{fcros-package}

\docType{package}

\title{FCROS for detecting differentially expressed genes}

\description{Implementation of a method based on fold change ranks ordering 
             statistics for detecting differentially expressed genes. FCROS 
             can be used with two biological conditions
             dataset: microarray or RNA-seq, ... Using combinations of
             pairs of samples from the two biological conditions, fold change
             (FC) are calculated. A matrix of fold changes or ratios can also
             be used as input. The FC are sorted in an increasing order and
             corresponding ranks are associated to genes. Then, a statistic
             is computed for the robust
             average ordered ranks of each gene in the dataset.
}

\details{
    \tabular{ll}{
        Package: \tab fcros\cr
        Type: \tab Package\cr
        Version: \tab 1.2\cr
        Date: \tab 2014-09-09\cr
        License: \tab GPL (>=2.0)\cr
    }
       
    Package fcros has the following functions:
    \tabular{ll}{
         fcros(): \tab Using dataset from two biological conditions and in a
                  single table, this function\cr\tab performs pairwise 
                  conparison of samples to obtain a matrix of fold changes.
                  \cr\tab Fold changes of this matrix are sorted
                  and used to associate statistic to genes.\cr
         fcros2(): \tab Using dataset from two biological conditions and in
                  two tables, this function\cr\tab performs pairwise
                  conparison of samples to obtain a matrix of fold changes.
                  \cr\tab Fold changes of this matrix are sorted
                  and used to associate statistic to genes.\cr
         fcrosMod(): \tab This function uses fold changes or ratios matrix to
                  associate statistic to genes.\cr
         fcrosFCmat(): \tab This function allows to compute a matrix with
                  pairwise fold changes from\cr\tab two biological conditions
                  dataset.\cr
         fcrosTtest(): \tab This function allows to use Student t-test to 
                       calculate p-value for genes in a\cr\tab dataset.\cr
         fcrosRead(): \tab Allows to read a tab delimited text file to be 
                  use in fcros() or fcros2().\cr
         fcrosWrite(): \tab Allows to save results from function fcros(),
                  fcros2() or fcrosMod() in a tab\cr\tab delimited text 
                  file.\cr
         fcrosTopN(): \tab Allows to search for the top N down- and/or
                  up-regulated genes from the\cr\tab results of function
                  fcros(), fcros2() or fcrosMod().\cr
         fcrosHistoPlot(): \tab Plots on the screen the histogram of the
                  FCROS statistics obtained using results\cr\tab of
                  function fcros(), fcros2() or fcrosMod()\cr
         fcrosVolcanoPlot(): \tab Perform a volcano plot of the FCROS
                  results from function fcros(), fcros2() or\cr\tab
                  fcrosMod()\cr
    }
}

\author{Doulaye Dembele
Maintainer: Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D and Kastner P, Fold change rank ordering statistics: 
                    a new method for detecting differentially expressed 
                    genes, BMC Bioinformatics 2014, 15:14}
\keyword{ package }

\examples{
   data(fdata);

   cont <- c("cont01", "cont07", "cont03", "cont04", "cont08");
   test <- c("test01", "test02", "test08", "test09", "test05");
   log2.opt <- 0;
   trim.opt <- 0.25;

   # perform fcros
   af <- fcros(fdata, cont, test, log2.opt, trim.opt);

   # perform Volcano plot
   fcrosVolcanoPlot(af, thr = 0.01)

   # save fcros values in a file
   fcrosWrite(af, file = "test2delete_values.txt");
   
   # save fcros parameters in a file
   fcrosWrite(af, file = "test2delete_params.txt",values = FALSE);

   # now select top 20 down and/or up regulated genes
   top20 <- fcrosTopN(af$f.value,20);
   alpha1 <- top20$alpha[1];
   alpha2 <- top20$alpha[2];
   id.down  <- matrix(c(0,11),ncol=1);
   id.up <- matrix(c(rep(0,11)),ncol=1);
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
       lines(t, data.down[i,2:21], type = "l", col = "blue")
   }

   # now plot down and up regulated genes
   plot(t, data.up[1,2:21], type = "l", col = "red", xlim = c(1,20),
       ylim = c(0,18), main = "Top up-regulated genes");
   for (i in 2:nup) {
       lines(t, data.up[i,2:21], type = "l", col = "red")
   }
   par(op)
}
