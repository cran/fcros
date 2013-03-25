\name{fcrosTopN}
\alias{fcrosTopN}

\title{fcrosTopN: Search for the "topN" differentially expressed genes}

\description{This function allows to seach for the "topN" DE genes.
It uses the fcros() f-value results for that purpose.}

\usage{fcrosTopN(fvalue, topN)}

\arguments{
  \item{fvalue}{ This is a f-value vector obtained using fcros() function: \code{fvalue}}
  \item{topN}{ Expected number of DE genes in the dataset:
       \code{topN}}
}

\value{ This function return a data frame containg 2 components
    \item{alpha }{Two thrsehold values allowing to have "topN" DE genes}
    \item{index }{Indexes of the "topN" DE genes}
}

\author{Doulaye Dembele doulaye@igbmc.fr}

\references{Dembele et al. (2013), Fold change ranks ordering statistics for detecting
                    differentially expressed genes, Manuscript in preparation}

\examples{
   data(fdata);

   cont <- c("cont01", "cont07", "cont03", "cont04", "cont08");
   test <- c("test01", "test02", "test08", "test09", "test05");
   log2.opt <- 0;

   # perform fcros
   af <- fcros(fdata, cont, test, log2.opt);
   
   # now select top 10 down and/or up regulated genes
   top10 <- fcrosTopN(af$f.value,10);

   # display thresholds
   top10$alpha
   
   # display index of top10 genes
   fdata[top10$index,1]
   
   # display fvalue of the top10 genes
   (af$f.value)[top10$index]
}
