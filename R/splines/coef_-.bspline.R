"coef<-.bspline" <- function(b, value){
   # COEF<-.BSPLINE - Coefficient assignement for a 'bspline' object.

   k <- b$dim - sum(b$is.fixed)
   if (length(value) != k) stop("'bspline' coefficient vector has ", k, " free values.")
   b$coefficients[!b$is.fixed] <- value

   return(b)
}
