integrate.piecewise <- function(p){
   # INTEGRATE.PIECEWISE - Returns the integral of a 'piecewise' object.

   # Extract piecewise coefficients:
   beta <- p$coefficients

   # Update polynomial coefficients:
   for (i in 1:dim(beta)[1]){
      beta[i, ] <- beta[i, ] / (1:dim(beta)[2])
   }
   beta <- cbind(rep(0, dim(beta)[1]), beta)

   # Calculate polynomial constants:
   beta[1,1] <- -1*evaluate(piecewise(coef = beta[1, ]), p$knots[1])
   if (dim(beta)[1] > 1){
      for (i in 2:dim(beta)[1]){
         b <- piecewise(coef = beta[i, ])
         beta[i,1] <- -1*evaluate(b, p$knots[i]) + evaluate(piecewise(coef = beta[i-1, ]), p$knots[i])
      }
   }
   
   # Update coefficient values:
   p$coefficients <- beta
   
   # Add constant line to minus infinity:
   if (is.finite(p$knots[1])){
      v <- evaluate(piecewise(coef = beta[1, ]), p$knots[1])
      p$knots <- c(-Inf, p$knots)
      p$coefficients <- rbind(c(v, rep(0, dim(p$coefficients)[2]-1)), p$coefficients)
   }
   
   # Add constant line to infinity:
   if (is.finite(p$knots[length(p$knots)])){
      v <- evaluate(piecewise(coef = p$coefficients[dim(p$coefficients)[1], ]), p$knots[length(p$knots)]) 
      p$knots <- c(p$knots, Inf)
      p$coefficients <- rbind(p$coefficients, c(v, rep(0, dim(p$coefficients)[2]-1)))
   }
   
   return(p)
}
