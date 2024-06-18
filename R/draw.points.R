draw.points <- function(n = 1, ...){
   # DRAW.POINTS - Draw points interactively on a plot.

   # Create default plot:
   if (is.null(dev.list())){
      plot(c(0, 1), c(0, 1), type = "n", xlab = "", ylab = "")
   }

   # Point vectors to be returned:
   xx <- NULL
   yy <- NULL

   # Initialize point counter variable:
   i <- 0
   
   plot.point <- function(buttons, x, y){
      # Convert to user coordinates:
      x <- grconvertX(x, "ndc", "user")
      y <- grconvertY(y, "ndc", "user")
      xx[i+1] <<- x
      yy[i+1] <<- y

      # Draw points:
      points(x, y, ...)

      # Increment counter:
      i <<- i + 1
      
      if (i == n) return(4) else return(NULL)
   }

   setGraphicsEventHandlers(prompt="Click in plot to draw points.",
                            onMouseDown = plot.point)

   getGraphicsEvent()

   return(data.frame(x = xx, y = yy))
}
