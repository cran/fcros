fcrosWrite <-
function(af, file = "filename", values = TRUE) {
   ## Write fcros results in tab delimited files
   
   n <- length(af$ri);
   geneid <- paste("g",as.character(1:n), sep="");
   if (values) {
      if (nrow(summary(af)) > 7) {
         af.values <- matrix(c(geneid, af$ri, af$FC, af$FC2, af$f.value, af$p.value), ncol = 6);
         colnames(af.values) <- c("geneID", "ri-stat", "FC", "FC2", "f-Value", "p-value");
      } else {
         af.values <- matrix(c(geneid, af$ri, af$FC2, af$f.value, af$p.value), ncol = 5);
         colnames(af.values) <- c("geneID", "ri-stat", "FC2", "f-Value", "p-value");
      }

      write.table(af.values, file, quote = FALSE, sep="\t", eol = "\n", col.names = TRUE, row.names = FALSE);

   } else {
      af.params <- matrix(c((af$bounds)[1], (af$bounds)[2], (af$params)[1], (af$params)[2], (af$params)[3]), ncol = 5);
      colnames(af.params) <- c("lb","ub","delta","mean","sd");

      write.table(af.params, file, quote = FALSE, sep="\t", eol = "\n", col.names = TRUE, row.names = FALSE);

   }
}
