degree.piecewise <- function(x){
   # DEGREE.PIECEWISE - Returns the degree of a 'piecewise' object.
   
   res <- dim(x$coefficients)[2] - 1

   return(res)
}
