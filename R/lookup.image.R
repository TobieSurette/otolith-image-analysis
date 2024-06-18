lookup.image <- function(I, x, y, n, interpolate = TRUE){
   # LOOKUP.IMAGE - Look-up image values along a specified linear trajectory.
   
   # Choose first if there are multiple images:
   image.index <- as.numeric(which(unlist(lapply(I, class)) == "image"))
   if (length(image.index) > 0){
      I.copy <- as.list(I)[image.index]
      I <- I[[image.index[1]]]
   }else{
      I.copy <- list(I)
   }
   
   # Define horizontal and vertical pixel dimensions:
   dx <- as.numeric(names(sort(table(diff(I$x)), decreasing = TRUE)[1]))
   dy <- as.numeric(names(sort(table(diff(I$y)), decreasing = TRUE)[1]))
   
   # Build lookup vectors from start and end points:
   if ((length(x) == 2) & (length(y) == 2)){
       box <- bbox(I)
       temp <- line.clip(x, y, box$x, box$y)
       x <- temp$x
       y <- temp$y
       if (missing(n)) n <- ceil(max(abs(diff(x)/dx), abs(diff(y)/dy))) + 1
       x <- seq(x[1], x[2], len = n)
       y <- seq(y[1], y[2], len = n)
   }

   # Expand length of 'x' or 'y' to the length of the other if one is a scalar:
   if (xor(length(x) == 1, length(y) == 1)){
      if (length(x)==1) x <- rep(x, length(y))
      if (length(y)==1) y <- rep(y, length(x))
   }
   
   # Image dimensions:
   dz <- dim(I$z)
   
   # Convert 'x' and 'y' to pixel coordinates:
   xp <- (((x - I$x[1]) / (I$x[length(I$x)] - I$x[1])) * (length(I$x)-1)) + 1
   yp <- (((y - I$y[1]) / (I$y[length(I$y)] - I$y[1])) * (length(I$y)-1)) + 1
   
   # Initialize result variable:
   z <- matrix(NA, nrow = length(xp), ncol = max(1, length(image.index)))
      
   # Linear average if pixel coordinates are not integers:
   if (interpolate){
      fx <- floor(xp)
      fy <- floor(yp)
      
      # Index of points which lie within the image bounds:
      index <- (fx >= 1) & (fx < dz[1]) & (fy >= 1) & (fy < dz[2])
      
      # Remove exterior points:
      fx <- fx[index]
      fy <- fy[index] 
           
      # Calculate pixel weights:
      wx <- 1 - (xp[index] - fx)
      wy <- 1 - (yp[index] - fy)
      
      for (i in 1:length(I.copy)){
         # Calculate weighted 'z' value:
         z[index, i] <- wx * wy * I.copy[[i]]$z[fy * dz[1] + fx] +
                     (1-wx) * wy * I.copy[[i]]$z[fy * dz[1] + fx + 1] + 
                     wx * (1-wy) * I.copy[[i]]$z[(fy+1) * dz[1] + fx] + 
                     (1-wx) * (1-wy) * I.copy[[i]]$z[fy * dz[1] + fx + 1]
      }
   }else{
      xp <- round(xp)
      yp <- round(yp)
      index <- (xp >= 1) & (xp <= dz[1]) & (yp >= 1) & (yp <= dz[2])
      for (i in 1:length(I.copy)){
         z[index, i] <- I.copy[[i]]$z[y[index] * dz[1] + xp[index]]
      }
   }

   # Prepare result:
   v <- data.frame(x = x, y = y, z = z)
   if (!is.null(names(I.copy))) names(v) <- c("x", "y", names(I.copy))
   
   return(v)
}
