plot.bspline <- function(b, legend = TRUE, xlim, ...){
   # PLOT.BSPLINE - Graphically display a 'bspline' object.

   # Extract full vector of knots:
   knots <- b$knots
   knots <- sort(unique(knots[is.finite(knots)]))
   
   # Define plot interval:
   if (missing(xlim)){
      interval <- NULL
      if (length(knots) == 0) interval <- c(0, 1)
      if (length(knots) == 1) interval <- c(knots-0.5, knots+0.5)
      if (length(knots) > 1)  interval <- c(min(knots), max(knots))
   }else{
      interval = xlim
   }
   
   # Define plotting values:
   x <- seq(interval[1], interval[2], len = 1000)

   # Plot basis functions:
   y <- evaluate(b, x)
  
   plot(interval, c(min(y, na.rm = TRUE), max(y, na.rm = TRUE)), 
        xlab = "x", ylab = "y", type = "n", ...)
   lines(x, y, lwd = 2)
   
   # Plot knot positions:
   for (i in 1:length(knots)){
      lines(c(knots[i], knots[i]), par("usr")[3:4], col = "red", lwd = 2, lty = "dashed")
   }
}

