plot.basis <- function(b, legend = TRUE){
   # PLOT.BASIS - Graphically display a 'basis' object.

   # Extract full vector of knots:
   knots <- NULL
   for (i in 1:length(b)){
      knots <- c(knots, b[[i]]$knots)
   }
   knots <- sort(unique(knots[is.finite(knots)]))

   # Define plot interval:
   interval <- NULL
   if (length(knots) == 0) interval <- c(0, 1)
   if (length(knots) == 1) interval <- c(knots-0.5, knots+0.5)
   if (length(knots) > 1)  interval <- c(min(knots), max(knots))

   # Define plotting values:
   x <- seq(interval[1], interval[2], len = 1000)

   # Plot basis functions:
   y <- evaluate(b, x)
   plot(interval, c(min(y, na.rm = TRUE), max(y, na.rm = TRUE)),
        xlab = "x", ylab = "y", type = "n")
   for (i in  1:length(b)){
      lines(x, y[, i], lwd = 2, col = rainbow(length(b))[i])
   }

   # Plot knot positions:
   for (i in 1:length(knots)){
      lines(c(knots[i], knots[i]), par("usr")[3:4], col = "red", lwd = 2, lty = "dashed")
   }

   # Display legend:
   if ((length(b) <= 10) & (legend == TRUE)){
      legend("topleft",
             paste("Basis", 1:length(b)),
             lwd = 2, col = rainbow(length(b)),
             bg = "white")
   }
}

