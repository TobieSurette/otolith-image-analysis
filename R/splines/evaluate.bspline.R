evaluate.bspline <- function(b, x){
   # EVALUATE.BSPLINE - Evaluate a 'bspline' object at specified values.

   y <- evaluate(b$basis, x) %*% as.matrix(b$coefficients)

   return(y)
}
