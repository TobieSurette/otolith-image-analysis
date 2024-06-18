resize.image <- function(I, dx, dy, xlim, ylim){
   # RESIZE - Resize an image.
   
   # Check input arguments:
   if (!missing(dx) & missing(dy)) dy <- dx

   # Truncate image:
   if (missing(dx) & !missing(xlim)){
      index <- (I$x >= xlim[1]) & (I$x <= xlim[2])
      I$x <- I$x[index]
      I$z <- I$z[index, ]
   }
   if (missing(dy) & !missing(ylim)){
      index <- (I$y >= ylim[1]) & (I$y <= ylim[2])
      I$y <- I$y[index]
      I$z <- I$z[, index]
   }
   
   if (missing(dx)) return(I)

   # Create coordinate matrices:
   xx <- floor((I$x - I$x[1]) / dx) * dx
   yy <- floor((I$y - I$y[1]) / dy) * dy
   xx <- outer(xx, rep(1, length(yy)))
   yy <- outer(rep(1, dim(xx)[1]), yy)

   # Calculate mean pixel values:
   temp <- aggregate(list(z = as.vector(I$z)),
                     by = list(x = as.vector(xx), y = as.vector(yy)),
                     mean, na.rm = TRUE)
   
   # Define image vectors:
   xi <- sort(unique(temp$x))
   yi <- sort(unique(temp$y))

   # Define dimensions of image matrix:
   zi <- temp$z
   dim(zi) <- c(length(xi), length(yi))
   zi <- t(zi)
   zi <- zi[dim(zi)[1]:1, ]

   # Create image:
   v <- as.image(zi, x = yi, y = xi)
   
   return(v)
}
