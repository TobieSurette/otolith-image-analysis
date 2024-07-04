print.piecewise <- function(p){
   # PRINT.PIECEWISE - Display a 'piecewise' object to the R console.
   
   cat(paste("knots : [", paste(p$knots, collapse = ", "), "]\n", sep = ""))
   cat(paste("coefficients : [", paste(p$coefficients[1,], collapse = ", "), "]\n", sep = ""))
   if (dim(p$coefficients)[1] > 1){
      for (i in 2:dim(p$coefficients)[1]){
         cat(paste("             : [", paste(p$coefficients[i,], collapse = ", "), "]\n", sep = ""))
      }
   }
}
