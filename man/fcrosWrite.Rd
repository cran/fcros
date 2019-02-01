\name{fcrosWrite}
\alias{fcrosWrite}

\title{Writing the fcros() or pfco() results in a tab delimited text file}

\description{This function creates a tab-delimited text file with the
results of the function fcros(), fcros2(), fcrosMod(), pfco() or pfcoMod().
The results are the values associated with the parameters
"idnames", "ri", "FC", "FC2", "f.value" and "p.value".}

\usage{fcrosWrite(af, file = "fcrosResults.txt", thr = 1)}

\arguments{
  \item{af}{ An output object of the functions fcros(), fcros2(), fcrosMod(), 
             pfco()\cr
             or pfcoMod():
             \code{af = fcros(xdata, cont, test, log2.opt, trim.opt)}\cr
             \code{af = pfco(xdata, cont, test, log2.opt, trim.opt)}}
  \item{file}{ The output file name:
             \code{file = "fcrosResults.txt"}}
  \item{thr}{ A threshold allowing to filter data based on p-values. Default
              setting is 1 (no filtering) \code{thr = 1}}
}

\value{ This function creates and saves a tab-delimited text file on the disk.}

\author{Doulaye Dembele doulaye@igbmc.fr}

\examples{
    data(fdata);

    rownames(fdata) <- fdata[,1];
    cont <- c("cont01", "cont07", "cont03", "cont04", "cont08");
    test <- c("test01", "test02", "test08", "test09", "test05");
    log2.opt <- 0;
    trim.opt <- 0.25;
    af <- fcros(fdata, cont, test, log2.opt, trim.opt);
    af2 <- pfco(fdata, cont, test, log2.opt, trim.opt);

    fcrosWrite(af, file = "test2delete_values.txt");
    fcrosWrite(af2, file = "test2delete2_values.txt");
}
