"degree<-.piecewise" <- function(p, value){
   # DEGREE<-.PIECEWISE - Degree assignement for a 'piecewise' object.
   
   if (value < degree(p)) stop("Cannot reduce degree of 'piecewise' object.")
   if (value > degree(p)){
      m <- matrix(0, nrow = dim(p$coefficients)[1], ncol = value - degree(p))
      p$coefficients <- cbind(p$coefficients, m)
   }
   
   return(p)
}
