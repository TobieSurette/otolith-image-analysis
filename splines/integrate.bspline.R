integrate.bspline <- function(b){
   # INTEGRATE.BSPLINE - Integrate a B-spline object.
   
   b$basis <- integrate(b$basis)
   b$degree <- b$degree + 1
   
   return(b)
}
