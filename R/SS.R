SS <- function(theta, x, y, plot = TRUE){
   # SS - Sum-of-squares for Warp-Trend-Periodic model.

   # Define basis functions:
   bw <- ispline(c(0, seq(0, 1, len = 6), 1, 1))        # Warp basis.
   bt <- bspline(c(0, 0, 0, seq(0, 1, len = 6), 1, 1, 1))  # Trend basis:
   bp <- bspline(c(0, 0,  seq(0, 1, len = 6), 1, 1))       # Periodic basis:

   # Number of required parameters for each basis type:
   nw <- length(bw)-1
   nt <- length(bt)
   np <- length(bp)-3
   
   # Parse parameter vector:
   theta.w <- c(0, theta[1:nw])
   theta <- theta[(nw+1):length(theta)]
   theta.t <- theta[1:nt]
   theta <- theta[(nt+1):length(theta)]
   n <- exp(theta[1]);
   theta <- theta[2:length(theta)]
   theta.p <- theta[1:np]
   theta <- theta[(np+1):length(theta)]
   sigma <- abs(theta[1])
   
   # Warp coordinates:
   theta.w <- exp(theta.w) / sum(exp(theta.w))
   xw <- evaluate(bw, x) %*% theta.w  # Map to new 'x' coordinates.

   # Predicted trend line component:
   yt <- evaluate(bt, xw) %*% theta.t
   #if (plot) lines(x, yt, col = "blue")
   
   # Define period spline component:
   theta.p <- c(theta.p, -theta.p[3:1])
   yp <- evaluate(bp, (xw*n) %% 1) %*% theta.p

   # Calculate sum-of-squares:
   mu <- yt + yp

   sigma <- exp(sigma * x)
   ll <- -0.5*length(z)*log(2*pi) - sum(log(sigma)) - sum((z - mu)^2 / (2*sigma^2))

   if (plot){
      plot(range(x), range(z), type = "n")
      points(xw, z)
      lines(xw, mu, col = "purple", lwd = 3)
      lines(xw, xw, col = "red")
      lines(xw, yp, col = "green")
      lines(xw, 10*sigma, col = "gold")
   }
   
   return(-ll)
}
