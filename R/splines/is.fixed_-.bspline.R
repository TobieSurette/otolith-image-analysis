"is.fixed<-.bspline" <- function(x, value){
   # IS.FIXED<-.BSPLINE - Fixed parameter assignement for a 'bspline' object.

   if ((length(value) == length(x$coefficients)) & (all(is.logical(value)))){
      x$is.fixed <- value
   }else{
      stop(paste("'is.fixed' must contain", length(x$coefficients), "logical values."))
   }

   return(x)
}
