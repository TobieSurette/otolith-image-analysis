as.function.piecewise <- function(p){
   # AS.FUNCTION.PIECEWISE - Returns a function which evaluates a 'piecewise' object.
   
   f <- function(x){
      y <- evaluate(p, x)
      return(y)
   }
   
   return(f)
}
