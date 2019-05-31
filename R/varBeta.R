varBeta <- function(x, trim.opt=0) {
   xbar <- mean(x, trim=trim.opt)
   tmp <- (x - xbar)^2
   xvar <- mean(tmp)
   return(xvar)
}
