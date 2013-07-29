\name{fcros-package}
\alias{fcros-package}

\docType{package}

\title{FCROS for detecting differentially expressed genes}

\description{Implementation of a method based on fold change ranks ordering statistics
             for detecting differentially expressed genes. FCROS can be used with
             two biological conditions transcriptomic dataset: microarray or RNA-seq.
             Using several combinations of pairs of samples from the two biological
             conditions, fold change (FC) are calculated. For each combination, the FC
             obtained are sorted in increasing order and corresponding ranks are associated
             to genes. Then, a statistic is associated to the average ordered rank of 
             each gene in the dataset(s).
}

\details{
    \tabular{ll}{
        Package: \tab fcros\cr
        Type: \tab Package\cr
        Version: \tab 1.0\cr
        Date: \tab 2013-07-29\cr
        License: \tab GPL (>=2.0)\cr
    }
       
    Package fcros has the following functions:
    \tabular{ll}{
         fcros(): \tab The main function based on fold change to associated statistics
                  to genes.\cr\tab Data are in a single table.\cr
         fcros2(): \tab The main function based on fold change to associated statistics
                  to genes.\cr\tab Data are in two separate tables.\cr
         fcrosTtest(): \tab Uses Student t-test to calculate p-value for genes in a dataset.\cr
         fcrosRead(): \tab Allows to read a tab delimited text file to be use in fcros() or fcros2().\cr
         fcrosWrite(): \tab Allows to save fcros() results in a tab delimited text file\cr
         fcrosTopN(): \tab Allows to search for the top N down- and/or up-regulated genes using\cr \tab
                      the function fcros() or fcros2() results\cr
    }
}

\author{Doulaye Dembele
Maintainer: Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D and Kastner P, Fold change ranks ordering statistics for detecting
                    differentially expressed genes, 2013, Manuscript submitted}
\keyword{ package }

\examples{
   data(fdata);

   cont <- c("cont01", "cont07", "cont03", "cont04", "cont08");
   test <- c("test01", "test02", "test08", "test09", "test05");
   log2.opt <- 0;
   trim.opt <- 0.2;

   # perform fcros
   af <- fcros(fdata, cont, test, log2.opt, trim.opt);

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
