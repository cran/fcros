fvalTopN <-
function(fval, topN) {
      n <- length(fval);
      s.fval <- sort(fval, method = "sh", index.return = TRUE);
      i <- 1;
      j <- 1;
      while ((i+j) < topN) {
            if ((1-s.fval$x[i+1]) > (s.fval$x[n-j])) {
               i <- i+1; 
            } else { j <- j+1; }
      }
      alpha <- c(s.fval$x[i], s.fval$x[n-j+1]);
      index <- c(s.fval$ix[1:i],s.fval$ix[(n-j+1):n]);
      list(alpha = alpha, index = index)
}
