pvalTopN <-
function(pval, topN) {
      s.pval <- sort(pval, method = "sh", index.return = TRUE);
      alpha <- s.pval$x[topN];
      index <- c(s.pval$ix[1:topN]);
      list(alpha=alpha, index=index)
}
