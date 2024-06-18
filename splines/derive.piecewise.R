derive.piecewise <- function(p){
   # DERIVE.PIECEWISE - Returns the derivative of a 'piecewise' object.

   # Extract piecewise coefficients:
   beta <- p$coefficients
   
   if (dim(beta)[2] > 1){
       beta <- beta[, 2:dim(beta)[2], drop = FALSE]
       for (i in 1:dim(beta)[1]){
          beta[i, ] <- 1:dim(beta)[2] * beta[i, ]
       }
       p$coefficients <- beta
   }else{
       p$coefficients <- 0 * p$coefficients
   }

   # Remove zero-value non-structural knot intervals:
   if (!is.finite(p$knots[1]) & all(p$coefficients[1,] == 0)){
      p$knots <- p$knots[2:length(p$knots)]
      p$coefficients <- p$coefficients[2:dim(p$coefficients)[1], ]
   }
   if (!is.finite(p$knots[length(p$knots)]) & all(p$coefficients[dim(p$coefficients)[1],] == 0)){
      p$knots <- p$knots[1:(length(p$knots)-1)]
      p$coefficients <- p$coefficients[1:(dim(p$coefficients)[1]-1), ]
   }
   
   return(p)
}
