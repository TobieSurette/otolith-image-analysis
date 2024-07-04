support.bspline <- function(x, ...){
   # SUPPORT - Returns the support of a piecewise object.

   return(c(min(x$knots), max(x$knots)))
}
