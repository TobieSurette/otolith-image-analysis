derive.basis <- function(b){
   # DERIVE.BASIS - Returns the derivative of a 'basis' object.

   # Integrate each basis function:
   for (i in 1:length(b)) b[[i]] <- derive(b[[i]])

   return(b)
}