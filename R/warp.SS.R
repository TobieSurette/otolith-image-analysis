warp.SS <- function(theta, t, z, zt, X, fint, plot = FALSE){
   # SS - Sum-of-squares for Warp-Trend-Periodic model.

   # Define basis functions:
   if (missing(X)){
      bw <- ispline(c(0, 0, seq(0, 1, len = 6), 1, 1))

      # Number of required parameters for each basis type:
      nw <- length(bw)

      # Design matrix:
      X <- cbind(x, evaluate(bw, x))
   }else{
      # Number of required parameters for each basis type:
      nw <- dim(X)[2]-1
   }
   
   if (missing(fint)) fint <- approxfun(x, zt)
   
   # Parse parameter vector:
   theta.w <- c(0, theta[1:nw])
   theta.w <- exp(theta.w) / sum(exp(theta.w))
   
   xw <- X %*% theta.w  # Map to new 'x' coordinates.
   SS <- sum((z-fint(xw))^2)

   if (plot){
      plot(range(t), range(z), type = "n")
      points(t, z)
      points(t, f(xw), pch = 21, bg = "red")
   }

   return(SS)
}
