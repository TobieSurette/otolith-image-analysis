polar.image <- function(I, x, y, n, angles, edges){
   # POLAR.IMAGE - Polar transform of an image.

   # Check input arguments:
   if (missing(x) | missing(y)) stop("'x' and 'y' coordinates must be supplied.")
   if (!xor(missing(n), missing(angles))) stop("Either 'n' or 'angles' must be specified.")
   if (missing(n)) n <- length(angles)
   if (missing(angles)){
      if (length(n) > 1) na <- n[2] else na <- n[1]
      angles <- seq(0, 2*pi, length = n + 1)[1:n]
   }
   n  <- n[1]

   # Look colour values:
   if (missing(edges)){
      Z <- matrix(NA, ncol = length(angles), nrow = n)

      # Determine maximum radius of the image:
      b <- bbox(I)
      r <- sqrt(diff(b$x)^2 + diff(b$y)^2)
      
      # Lookup rdial values:
      for (i in 1:length(angles)){
         Z[,i] <- lookup(I, c(x, x + r * cos(angles[i])), c(y, y + r * sin(angles[i])), n = n)$z
      }
      Z <- as.image(Z[dim(Z)[1]:1, ], y = angles, x = 1:n)
   }else{
      Z <- matrix(NA, ncol = length(edges$x), nrow = n)
      names(edges) <- tolower(names(edges))
      for (i in 1:length(edges$x)){
         Z[,i] <- lookup(I, c(x, edges$x[i]), c(y, edges$y[i]), n = n)$z
      }
      Z <- as.image(Z[dim(Z)[1]:1, ], y = 1:length(edges$x), x = 1:n)
   }

   return(Z)
}
