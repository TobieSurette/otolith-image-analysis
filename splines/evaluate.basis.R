evaluate.basis <- function(b, x){
   # EVALUATE.BASIS - Evaluate a 'basis' object at specified values.

   y <- matrix(NA, nrow = length(x), ncol = length(b))
   for (i in 1:length(b)){
      y[,i] <- evaluate(b[[i]], x)
   }
   y[is.na(y)] <- 0

   return(y)
}
