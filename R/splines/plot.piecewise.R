plot.piecewise <- function(p){
   # PLOT.PIECEWISE - Graphically display a 'piecewise' object.
   
   # Define plot interval:
   if (all(is.finite(p$knots))){
      x <- seq(min(p$knots), max(p$knots), len = 1000) 
   }else{
      k <- p$knots[is.finite(p$knots)]
      if (length(k) > 1){
         x <- seq(min(k) - diff(range(k))/10, max(k) + diff(range(k))/10, len = 1000)
      }else{
         x <- seq(0, 1, len = 1000) 
      }
   }
   
   # Plot piecewise polynomial evaluation:
   plot(x, evaluate(p, x), 
        xlab = "x", ylab = "y", 
        type = "l", lwd = 2)
   
   # Plot finite knot points:
   k <- unique(p$knots)
   for (i in 1:length(k)){
      lines(c(k[i], k[i]), par("usr")[3:4], col = "red", lwd = 2, lty = "dashed")
   }
}
