\name{pvalTopN}
\alias{pvalTopN}

\title{Search for the top N changed genes or probes using p-values}

\description{This function allows to seach for the top N 
differentially expressed genes or changed probes.
It uses the p-value. These p-values can be obtained using one of
the following functions fcros(), fcros2(), fcrosMod(), pfco(),
pfcoMod() or another statistical method.}

\usage{pvalTopN(pval, topN)}

\arguments{
  \item{pval}{ This is vector with p-values obtained using the functions fcros(),
        fcros2(), fcrosMod(), pfco(), pfcoMod() or another statistical method:
        \code{pval = af$p.value}}
  \item{topN}{ The expected number of the top DE genes/probes in the dataset used:
               \code{topN}}
}

\value{ This function returns a data frame containing 2 components
    \item{alpha }{Two threshold values for the down- and the up-regulated allowing
                to have the top N DE genes}
    \item{index }{The indexes of the top N DE genes}
}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele D and Kastner P, Fold change rank ordering statistics: 
            a new method for detecting differentially expressed genes,
            BMC Bioinformatics, 2014, 15:14\cr

            Dembele D and Kastner P, Comment on: Fold change rank ordering statistics:
                    a new method for detecting differentially expressed
                    genes, BMC Bioinformatics, 2016, 17:462}

\examples{
   data(fdata);

   rownames(fdata) <- fdata[,1]
   cont <- c("cont01", "cont07", "cont03", "cont04", "cont08");
   test <- c("test01", "test02", "test08", "test09", "test05");
   log2.opt <- 0;

   # perform fcros()
   af <- fcros(fdata, cont, test, log2.opt);
   
   # now select top 10 down and/or up regulated genes
   top10 <- pvalTopN(af$p.value, 12);

   # display thresholds
   top10$alpha
   
   # display index of top10 genes
   fdata[top10$index, 1]
   
   # display fvalue of the top10 genes
   (af$f.value)[top10$index]
}
