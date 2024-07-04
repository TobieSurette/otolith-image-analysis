mspline <- function(knots, ...){
   # MSPLINE - Evaluate m-spline basis functions.
   
   B <- bspline(knots, ...)
   I <- integrate(B)
   for (i in 1:length(I$basis)){
      k <- I$basis[[i]]$knots
      w <- evaluate(I$basis[[i]], k[length(k)])
      B$basis[[i]] <- (1/w)*B$basis[[i]]
   }
   
   return(B)
}
