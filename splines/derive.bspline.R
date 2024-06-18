derive.bspline <- function(b){
   # DERIVE.BSPLINE - Returns the derivative of a 'bspline' object.

   b$basis <- derive(b$basis)
   b$degree <- max(0, b$degree - 1)
   
   return(b)
}
