nucleus.otolith <- function(I){
   # NUCLEUS - Returns the centroid of an image using the colour intensities as weights.
   
   # If 'nucleus' is defined, return it:
   if (!is.null(I$nucleus)) return(I$nucleus)
   
   # Initial nucleus estimate:
   wx <- apply(I$image$z, 1, sum, na.rm = TRUE)
   wy <- apply(I$image$z, 2, sum, na.rm = TRUE)
   wx <- wx / sum(wx)
   wy <- wy / sum(wy)
   p <- NULL
   p[1] <- sum(wx * I$x)
   p[2] <- sum(wy * I$y)

   # Calculate gradient:
   G <- as.data.frame(gradient(I$image, polar = FALSE))
   
   # Nucleus sum-of-squares function:
   SS <- function(p, x, y, dx, dy){
      # Define vector from proposed nucleus to each pixel:
      xn <- x - p[1]
      yn <- y - p[2]
      
      # Calculate distance vector:
      dd <- sqrt(xn * xn + yn * yn)
      
      # Convert to unit vectors:
      xu <- xn / dd
      yu <- yn / dd
      
      # Perform projection:
      v <- (dx * xu) + (dy * yu)
      
      # Calculate sum-of-squares value:
      SS <- sum(v^2, na.rm = TRUE)
      
      # points(p[1], p[2], col = "red")
      return(SS)
   }

   p <- optim(p, SS, x = G$x, y = G$y, dx = G$dx, dy = G$dy, control = list(trace = 3, fnscale = -1))$par

   return(p)
}
