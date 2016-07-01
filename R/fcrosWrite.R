fcrosWrite <-
function(af, file = "fcrosResults.txt", values = TRUE, thr = 1) {
   n <- length(af$idnames);
   if (values) {
      id <- 1:n
      idx <- id[af$p.value <= thr];
      idnames <- af$idnames[idx];
      FC2 <- af$FC2[idx];
      fVal <- af$f.value[idx];
      pVal <- af$p.value[idx];
      if (length(af$ri) > 0) {
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
      }
      else {
         u1 <- af$u1[idx];
         if (nrow(summary(af)) > 8) {
            FC <- af$FC[idx];
            results <- matrix(c(as.character(idnames),u1, FC, FC2, fVal, pVal), ncol=6);
            colnames(results) <- c("idnames", "u1", "FC", "FC2", "f.value", "p.value");
         }
         else {
              results <- matrix(c(as.character(idnames), u1, FC2, fVal, pVal), ncol=5);
              colnames(results) <- c("idnames", "u1", "FC2", "f.value", "p.value");
         }
      }
      write.table(results, file, quote = FALSE, sep = "\t", eol = "\n", col.names = TRUE, row.names = FALSE);
   }
   else {
      if (length(af$ri) > 0) {
         results <- matrix(c((af$bounds)[1], (af$bounds)[2], (af$params)[1], (af$params)[2], (af$params)[3]), ncol=5);
         colnames(results) <- c("lb","ub","delta","mean","sd");
      }
      else {
         results <- matrix(c(af$comp, af$comp.w, af$comp.wcum), ncol=3);
         colnames(results) <- c("comp", "comp.w", "comp.wcum");
      }
      write.table(results, file, quote = FALSE, sep="\t", eol = "\n", col.names = TRUE, row.names = FALSE);
   }
}
