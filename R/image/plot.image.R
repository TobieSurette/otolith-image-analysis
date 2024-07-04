plot.image <- function(I, col, breaks, useRaster = TRUE, adjust = FALSE, xlab = "", ylab = "", ...){
   # PLOT.IMAGE - Image default plotting method.

   if ("image" %in% class(I)){
      if (missing(col)) col = gray(seq(0, 1, len = 1000))
      if (missing(breaks)) breaks = seq(min(I$z, na.rm = TRUE), max(I$z, na.rm = TRUE), len = length(col)+1)
      z <- I$z
      y <- I$y
      x <- I$x
      
      # Adjust aspect ratio:
      if (adjust){
         dx <- dim(I$z)[1]
         dy <- dim(I$z)[2]
         r <- dy / dx
         if (r < 1) par(pin = c(par("pin")[1], r*par("pin")[1]))
         if (r > 1) par(pin = c((1/r)*par("pin")[1], par("pin")[1]))
      }
      
      # Display image:
      step <- c(mean(diff(I$x)), mean(diff(I$y)))
      xlim <- c(min(I$x)-step[1], max(I$x)+step[1])
      ylim <- c(min(I$y)-step[2], max(I$y)+step[2])
      plot(xlim, ylim, type = "n", xlab = xlab, ylab = ylab, ...)
      
      graphics::image.default(x, y, z, useRaster = useRaster, breaks = breaks, col = col, add = TRUE, ...)
      box()
   }else{
      I <- as.image(I, ...)
      image(I, ...)
   }
}

