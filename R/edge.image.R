edge.image <- function(I, x, y, angles, n, alpha = 0.5, threshold, plot = FALSE){
   # EDGE.IMAGE - Find edge of an object determined by a threshold.
   
   # Flag which specifies whether to redistribute the points more evenly:
   if (missing(angles)) flag <- TRUE else flag <- FALSE
   
   # Determine center of mass of image:
   if (missing(x) & missing(y)){
      wx <- apply(I$z, 1, sum, na.rm = TRUE)
      wx <- wx / sum(wx)
      wy <- apply(I$z, 2, sum, na.rm = TRUE)
      wy <- wy / sum(wy)
      x <- sum(wx * I$x)
      y <- sum(wy * I$y)
   }

   # Convert 'x' and 'y' values to pixel coordinates:
   x <- (((x - I$x[1]) / (I$x[length(I$x)] - I$x[1])) * (length(I$x)-1)) + 1
   y <- (((y - I$y[1]) / (I$y[length(I$y)] - I$y[1])) * (length(I$y)-1)) + 1
   
   # Check input arguments:
   if (!xor(missing(n), missing(angles))) stop("Either 'n' or 'angles' must be specified.")

   # Define 'angles' if none is specified:
   if (missing(angles)) angles <- seq(0, 2*pi, length = n + 1)[1:n]

   # Plot original image:
   if (plot) plot(I)

   # Initialize edge list variable:
   edges <- data.frame(x = rep(NA, length(angles)),
                       y = rep(NA, length(angles)))

   # Loop over angles:
   for (i in 1:length(angles)){
      # Determine image edge coordinates ar angle i:
      if (abs(cos(angles[i])) >= cos(pi/4)){
         if (cos(angles[i]) > 0) xt <- round(x):dim(I$z)[1]
         if (cos(angles[i]) < 0) xt <- round(x):1
         yt <- round(y) + (sin(angles[i])/sin(pi/4))*(0:(length(xt)-1))
      }
      if (abs(cos(angles[i])) < cos(pi/4)){
         if (sin(angles[i]) > 0) yt <- round(y):dim(I$z)[2]
         if (sin(angles[i]) < 0) yt <- round(y):1
         xt <- round(x) + (cos(angles[i])/cos(pi/4))*(0:(length(yt)-1))
      }

      # Determine coordinates along ray and lookup values:
      xt <- round(xt)
      yt <- round(yt)
      index <- (xt >= 1) & (yt >= 1) & (xt <= dim(I$z)[1]) & (yt <= dim(I$z)[2])
      xt <- xt[index]
      yt <- yt[index]
      zt <- I$z[yt * dim(I$z)[1] + xt]  # Lookup 'z' values.

      # Prepare vector for logistic analysis:
      if (!missing(threshold)) zt[zt < threshold] <- NA
      if (any(is.na(zt))) zi <- is.na(zt) + 1 - 1 else zi <- zt
      edge <- switch.point(zi, alpha = alpha)
      w <- edge - floor(edge)

      # Calculate edge values:
      edges$x[i] <- (1-w)*xt[floor(edge)] + w*xt[ceil(edge)]
      edges$y[i] <- (1-w)*yt[floor(edge)] + w*yt[ceil(edge)]

      # Plot ray:
      if (plot) points(I$x[xt[1:edge]], I$y[yt[1:edge]], cex = zt[1:edge], col = "blue")
      
      # Convert 'x' and 'y' and 'edges' to native coordinates:
      edges$x[i] <- (edges$x[i] / dim(I$z)[1]) * (I$x[length(I$x)] - I$x[1]) + I$x[1]
      edges$y[i] <- (edges$y[i] / dim(I$z)[2]) * (I$y[length(I$y)] - I$y[1]) + I$y[1]
   }

   # Plot edge values:
   if (plot) points(edges$x, edges$y, pch = 21, bg = "red")
      
   # Redistribute points more evenly along the contour:
   if (flag){
      d <- sqrt((edges$x[2:n] - edges$x[1:(n-1)])^2 + (edges$y[2:n] - edges$y[1:(n-1)])^2)
      d <- c(0, cumsum(d))
      di <- seq(0, d[n], length = n)
      edges$x <- approx(d, edges$x, di)$y
      edges$y <- approx(d, edges$y, di)$y
   }

   return(edges)
}
