\name{fcrosWrite}
\alias{fcrosWrite}

\title{This function allows to write FCROS results in a file}

\description{This function will create a file containing FCROS results
This file contains values ("ri", "FC", "FC2", "f-value" and "p-value") or
parameters ("lb", "ub", "delta", "meand" and "sd").}

\usage{fcrosWrite(af, file = "filename", values = TRUE)}

\arguments{
  \item{af}{ An output object of the functions fcros() or fcros2():
             \code{af} = fcros(x$xdata, cont, test, log2.opt, trim.opt)}
  \item{file}{ The output file name:
             \code{file} = "test2delete_values.txt"}
  \item{values}{ A flag allowing to create a file with values "ri", "FC", 
                 "FC2", "f-value" and "p-value" (default) or with parameters
                 "lb", "ub", "delta", "meand" and "sd": \code{value} = TRUE}
}

\value{ This function creates and saves a tab delimited text file on the 
current disk.}

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
