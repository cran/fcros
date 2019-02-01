fcrosWrite <-
function(af, file = "fcrosResults.txt", thr = 1) {
   n <- length(af$idnames);
   id <- 1:n
   idx <- id[af$p.value <= thr];
   idnames <- af$idnames[idx];
   FC2 <- af$FC2[idx];
   fVal <- af$f.value[idx];
   pVal <- af$p.value[idx];
   ri <- af$ri[idx];
   if (nrow(summary(af)) > 8) {
      FC <- af$FC[idx];
      results <- matrix(c(as.character(idnames),ri, FC, FC2, fVal, pVal), ncol=6);
      colnames(results) <- c("idnames", "ri", "FC", "FC2", "f.value", "p.value");
   }
   else {
      results <- matrix(c(as.character(idnames), ri, FC2, fVal, pVal), ncol=5);
      colnames(results) <- c("idnames", "ri", "FC2", "f.value", "p.value");
   }
   write.table(results, file, quote = FALSE, sep = "\t", eol = "\n", col.names = TRUE, row.names = FALSE);
}
