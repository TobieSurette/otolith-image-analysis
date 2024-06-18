plot.otolith <- function(O){
   # PLOT.OTOLITH - Graphically display an 'otolith' object.
   
   # Display image:
   plot(O$image)
   
   # Display nucleus:
   if (!is.null(O$nucleus)) points(O$nucleus$x, O$nucleus$y, pch = 21, bg = "red", cex = 2)
   
   # Display contour:
   if (!is.null(O$edges)){
      if (!is.null(O$edges$bspline)){
         t <- seq(0, 1, len = 1000)
         t <- t[1:(length(t)-1)]
         x <- evaluate(O$edges$bspline$x, t)
         y <- evaluate(O$edges$bspline$y, t)
         lines(x, y, lwd = 3, col = "red")
      }else{
         points(O$edges$x, O$edges$y, pch = 21, bg = "red", cex = 1)
      }
   }
}
