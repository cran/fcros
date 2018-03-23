\name{fcros-package}
\alias{fcros-package}

\docType{package}

\title{A Method to Search for Differentially Expressed Genes and to Detect Recurrent
             Chromosomal Copy Number Aberrations}

\description{Implementation of a method based on fold change rank ordering
             statistics to search for differentially expressed genes or to
             detect recurrent chromosomal copy number aberrations. This package
             can be used for two biological conditions high-throughput dataset 
             (microarray, RNA-seq, ...), for expression profiling dataset over time
             without replicates or for cytogenetics dataset (aCGH, Sequencing).
}

\details{
    \tabular{ll}{
        Package: \tab fcros\cr
        Type: \tab Package\cr
        Version: \tab 1.5.6\cr
        Date: \tab 2018-03-21\cr
        License: \tab GPL (>= 2)\cr
    }
       
    Package fcros has the following functions:
    \tabular{ll}{
         fcros(): \tab The function to use with a dataset from two biological condition samples. The\cr\tab
                  dataset should be in a single table. The function fcros() performs a pairwise\cr\tab
                  conparison of samples to obtain a matrix of fold changes. The fold changes\cr\tab
                  are sorted and their rank values are used to associate statistic with genes/probes.\cr
         fcros2(): \tab The function to use with datasets from two biological biological conditions. The\cr\tab
                  datasets should be in two separate tables as inputs. The\cr\tab
                  function fcros2() performs a pairwise comparison of samples from each table\cr\tab to
                  obtain fold changes. The fold changes are sorted, their rank values are\cr\tab
                   combined and then used to associate statistic with genes/probes.\cr
         pfco(): \tab The function to use with a dataset from two biological condition samples. The\cr\tab
                  dataset should be in a single table. The function pfco() performs a pairwise\cr\tab
                  conparison of samples to obtain a matrix of fold changes. The fold changes\cr\tab
                  are sorted and their rank values are used to associate statistic with genes/probes\cr
                  \tab using a singular value decomposition.\cr
         fcrosMod(): \tab This function uses fold changes or ratios matrix as input to associate
                  statistic\cr\tab with genes/probes.\cr
         pfcoMod(): \tab This function uses fold changes or ratios matrix as input to associate
                  statistic\cr\tab with genes/probes using a singular value decomposition.\cr
         fcrosFCmat(): \tab This function allows to compute a matrix of fold changes using\cr\tab
                   pairwise comparisons of the two biological condition samples in a dataset.\cr
         fcrosTtest(): \tab This function allows to use the Student t-test to calculate p-values\cr\tab 
                  for the genes in a dataset.\cr
         fcrosRead(): \tab This function allows to read a tab delimited text file to be use as an\cr\tab 
                  input for the function fcros(), fcros2() or fcrosMod().\cr
         fcrosWrite(): \tab This function allows to save the results obtained using the function
                  fcros(),\cr\tab fcros2() or fcrosMod() in a tab delimited text file.\cr
         fcrosTopN(): \tab This function allows to search for the top N down- and/or up-regulated genes
                  \cr\tab from the results obtained using the function fcros(), fcros2(), pfco(),\cr\tab
                  fcrosMod() or pfcoMod().\cr
         fvalTopN(): \tab This function allows to search for the top N down- and/or up-regulated genes
                  \cr\tab from the results obtained using the function fcros(), fcros2(), pfco(),\cr\tab
                  fcrosMod() or pfcoMod().\cr
         pvalTopN(): \tab This function allows to search for the top N down- and/or up-regulated genes
                  \cr\tab from the results obtained using the function fcros(), fcros2(), pfco(),\cr\tab
                  fcrosMod() or pfcoMod().\cr
         histoPlot(): \tab This function plots on the screen the histogram of the FCROS statistics
                  \cr\tab obtained using the results of the function fcros(), fcros2(), pfco(),\cr\tab
                  fcrosMod() or pfcoMod()\cr
         fvalVolcanoPlot(): \tab This function performs a volcano plot of the results obtained\cr\tab
                  using the function fcros(), fcros2(), pfco(), fcrosMod() or pfcomod()\cr
         pvalVolcanoPlot(): \tab This function performs a volcano plot of the results obtained\cr\tab
                  using the function fcros(), fcros2(), pfco(), fcrosMod() or pfcoMod()\cr
         chrSummary(): \tab This function summarizes detection results by chromosome\cr
         chrSegment(): \tab This function segments a chromosome data\cr
         chrPlot(): \tab This function performs a plot of the chromosome probes data\cr
         chrPlot2(): \tab This function performs a plot of the chromosome segmentation results\cr
         voomReads(): \tab This function performs a transformation of the read counts\cr
         tcnReads(): \tab This function performs a total count normalization of reads\cr
         rankReads(): \tab The function to use with a dataset from two biological condition samples. 
                  The\cr\tab dataset should be in a single table. The function rankReads() performs a\cr\tab 
                  pairwise conparison of samples to obtain a matrix of fold changes. Small uniform\cr\tab
                  values are added to read counts. This is repeated nrun time. The fold changes\cr\tab
                  are sorted and their rank values are used to associate statistic with genes/probes.
    }
}

\author{Doulaye Dembele
Maintainer: Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D and Kastner P, Fold change rank ordering statistics: 
                    a new method for detecting differentially expressed 
                    genes, BMC Bioinformatics, 2014, 15:14\cr

            Dembele D, Analysis of high biological data using their rank
                    values, Stat Methods Med Res, 2018}

\keyword{ package }

\examples{
   data(fdata);

   cont <- c("cont01", "cont07", "cont03", "cont04", "cont08");
   test <- c("test01", "test02", "test08", "test09", "test05");
   log2.opt <- 0;
   trim.opt <- 0.25;

   # perform fcros()
   af <- fcros(fdata, cont, test, log2.opt, trim.opt);

   # perform Volcano plot
   fvalVolcanoPlot(af, thr = 0.01)

   # save fcros values in a file
   fcrosWrite(af, file = "test2delete_values.txt");
   
   # save fcros parameters in a file
   fcrosWrite(af, file = "test2delete_params.txt", values = FALSE);

   # now select top 20 down and/or up regulated genes
   top20 <- fcrosTopN(af, 20);
   alpha1 <- top20$alpha[1];
   alpha2 <- top20$alpha[2];
   id.down  <- matrix(c(0,11), ncol = 1);
   id.up <- matrix(c(rep(0,11)), ncol = 1);
   n <- length(af$FC);
   f.value <- af$f.value;

   idown <- 1;
   iup <- 1;
   for (i in 1:n) {
       if (f.value[i] <= alpha1) { id.down[idown] <- i; idown <- idown+1; }
       if (f.value[i] >= alpha2) { id.up[iup] <- i; iup <- iup+1; }
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
