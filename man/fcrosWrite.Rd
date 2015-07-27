\name{fcrosWrite}
\alias{fcrosWrite}

\title{fcrosWrite(): to write the fcros results in file}

\description{This function creates a tab-delimited text file with the
results of the function fcros(), fcros2() or fcrosMod().
Two kind of file can be obtained. In the first, the results of the parameters
"idnames", "ri", "FC", "FC2", "f-value" and "p-value" are saved. In the second,
the results of the parameters "lb", "ub", "delta", "meand" and "sd" are saved.}

\usage{fcrosWrite(af, file = "fcrosResults.txt", values = TRUE, thr = 1)}

\arguments{
  \item{af}{ An output object of the functions fcros(), fcros2() or fcrosMod():
             \code{af = fcros(xdata, cont, test, log2.opt, trim.opt)}}
  \item{file}{ The output file name:
             \code{file = "fcrosResults.txt"}}
  \item{values}{ A flag allowing to create a file with the parameters "ri", "FC",
                 "FC2", "f-value" and "p-value" (default) or with the parameters
                 "lb", "ub", "delta", "meand" and "sd": \code{values} = TRUE}
  \item{thr}{ A threshold allowing to filter data based on p-values. Default
              setting is 1 (no filtering) \code{thr = 1}}
}

\value{ This function creates and saves a tab-delimited text file on the disk.}

\author{Doulaye Dembele doulaye@igbmc.fr}

\examples{
    data(fdata);
    cont <- c("cont01", "cont07", "cont03", "cont04", "cont08");
    test <- c("test01", "test02", "test08", "test09", "test05");
    log2.opt <- 0;
    trim.opt <- 0.25;
    af <- fcros(fdata, cont, test, log2.opt, trim.opt);

    fcrosWrite(af, file = "test2delete_values.txt");
    fcrosWrite(af, file = "test2delete_params.txt", values = FALSE);
}
