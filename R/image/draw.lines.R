draw.lines <- function(n = 1, ...){
   # DRAW.LINES - Draw lines interactively on a plot.

   # Create default plot:
   if (is.null(dev.list())){
      plot(c(0, 1), c(0, 1), type = "n", xlab = "", ylab = "")
   }

   # Initialize result point vectors:
   xx <- rep(NA, n + 1)
   yy <- rep(NA, n + 1)

   # Initialize point counter variable:
   i <- 0
   
   # Define point plotting method:
   plot.point <- function(buttons, x, y){
      # Convert to user coordinates:
      x <- grconvertX(x, "ndc", "user")
      y <- grconvertY(y, "ndc", "user")
      xx[i+1] <<- x
      yy[i+1] <<- y

      # Draw line if there are at least two points specified:
      if (i > 0) lines(c(xx[i],  xx[i+1]), c(yy[i],  yy[i+1]), ...)
         
      # Draw points:
      points(x, y, ...)

      # Increment counter:
      i <<- i + 1

      if (i == (n+1)) return(4) else return(NULL)
   }

   setGraphicsEventHandlers(prompt="Click in plot to draw points.",
                            onMouseDown = plot.point)

   getGraphicsEvent()

   return(data.frame(x = xx, y = yy))
}
