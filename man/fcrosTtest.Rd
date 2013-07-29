\name{fcrosTtest}
\alias{fcrosTtest}

\title{Student t-test for detecting differentially expressed genes}

\description{The function uses the basic R t.test() function to perform the
             Student t-test for genes in a data matrix. It should be used
             for two biological conditions transcriptomic dataset (microarray,
             or RNA-seq). Fold change and p-values are returned
             for each gene in the data matrix.}

\usage{fcrosTtest(xdata, cont, test, log2.opt = 0)}

\arguments{
  \item{xdata}{ A matrix or table containing a two biological conditions data to
                process for detecting differentially expressed genes:
       \code{xdata}}
  \item{cont}{ A vector containing the names of control samples:
       \code{cont} = c("cont01", "cont02", ...)}
  \item{test}{ A vector containing the names of test samples:
       \code{test} = c("test01", "test02", "test03", ...)}
  \item{log2.opt}{ A scalar equals to 0 or 1. Value 0 (default) means that
              data in "xdata" are expressed in a log2 scale: \code{log2.opt} = 0}
}

\details{Names appearing in the parameters "cont" and "test" should match
column names of the data matrix "xdata".}

\value{     
    \item{FC }{Fold changes for genes in the data matrix.}
    \item{p.value }{Student t-test p-value associated to genes.}
}

\author{Doulaye Dembele}

\examples{
   data(fdata);

   cont <- c("cont01", "cont07", "cont03", "cont04", "cont08");
   test <- c("test01", "test02", "test08", "test09", "test05");
   log2.opt <- 0;

   # perform fcros
   at <- fcrosTtest(fdata, cont, test, log2.opt);
   
   # now select some differentially expressed genes
   id.de  <- matrix(0,1);
   n <- length(at$FC);
   for (i in 1:n) {
       if ((at$p.value)[i] <= 0.0005) { id.de <- rbind(id.de,i); }
   }

   data.de <- fdata[id.de,];
   nde <- nrow(data.de);

   # now plot the DE genes
   t <- 1:20;
   plot(t, data.de[1,2:21], type="l", col="blue", xlim=c(1,20), 
        ylim=c(0,18), main="Down- and up-regulated genes");
   for (i in 2:nde) {
       lines(t, data.de[i,2:21], type = "l", col = "blue")
   }
}
