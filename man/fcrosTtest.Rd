\name{fcrosTtest}
\alias{fcrosTtest}

\title{Student t-test for detecting differentially expressed genes}

\description{The function uses the basic R t.test() function to perform 
             the Student t-test. It should be
             used for two biological conditions dataset
             (microarray, or RNA-seq). The Fold changes, statistics and p-values
             are returned for each gene in the dataset.}

\usage{fcrosTtest(xdata, cont, test, log2.opt = 0)}

\arguments{
  \item{xdata}{ A table containing a two biological conditions
          dataset to process for detecting differentially expressed genes.
          The rownames of xdata are used for the output idnames.}
  \item{cont}{ A vector containing the label names of the control samples:
       \code{cont} = c("cont01", "cont02", ...)}
  \item{test}{ A vector containing the label names of the test samples:
       \code{test} = c("test01", "test02", "test03", ...)}
  \item{log2.opt}{ A scalar equals to 0 or 1. The value 0 (default) means that
       data in the matrix "xdata" are expressed in a log2 scale: \code{log2.opt} = 0}
}

\details{Label names appearing in the parameters "cont" and "test" should match
column label names of the data matrix "xdata". It is not necessary to use all
column label names of the dataset "xdata".}

\value{     
    \item{idnames}{ A vector containing the list of IDs or symbols associated with genes}
    \item{FC }{The fold changes for the genes in the dataset.}
    \item{stat }{The Student t-test statistics associated with genes.}
    \item{p.value }{The Student t-test p-values associated with genes.}
}

\author{Doulaye Dembele}

\examples{
   data(fdata);

   rownames(fdata) <- fdata[,1];

   cont <- c("cont01", "cont07", "cont03", "cont04", "cont08");
   test <- c("test01", "test02", "test08", "test09", "test05");
   log2.opt <- 0;

   # perform fcrosTtest()
   at <- fcrosTtest(fdata, cont, test, log2.opt);
   
   # now select some differentially expressed genes
   id.de  <- matrix(0, 1);
   n <- length(at$FC);
   for (i in 1:n) {
       if ((at$p.value)[i] <= 0.0005) { id.de <- rbind(id.de, i); }
   }

   data.de <- fdata[id.de, ];
   nde <- nrow(data.de);

   # now plot the DE genes
   t <- 1:20;
   plot(t, data.de[1, 2:21], type = "l", col = "blue", xlim = c(1,20),
        ylim = c(0,18), main = "Down- and up-regulated genes");
   for (i in 2:nde) {
       lines(t, data.de[i,2:21], type = "l", col = "blue")
   }
}
