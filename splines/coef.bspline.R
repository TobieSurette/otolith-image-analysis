coef.bspline <- function(x, fixed = NULL){
   # COEF.BSPLINE - Return the coefficients of a 'bspline' object.
   
   if (is.null(fixed)) return(x$coefficients)
   if (fixed) return(x$coefficients[x$is.fixed])
   if (!fixed) return(x$coefficients[!x$is.fixed])
}
