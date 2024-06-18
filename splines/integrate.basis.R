integrate.basis <- function(b){
   # INTEGRATE.BASIS - Integral of a 'basis' object.

   # Integrate each basis function:
   for (i in 1:length(b)) b[[i]] <- integrate(b[[i]])

   return(b)
}
