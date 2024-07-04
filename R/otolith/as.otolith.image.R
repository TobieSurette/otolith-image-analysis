as.otolith.image <- function(I, species, threshold, n = 100, rotate = FALSE){
   # AS.OTOLITH.IMAGE - Convert an image to an otolith object.
   
   # Extract otolith edge polygon:
   edges <- edge(I, threshold = threshold, n = n)

   # Screen out background data:
   temp <- as.data.frame(I)
   index <- in.polygon(as.polygon(edges$x, edges$y), temp$x, temp$y)
   temp <- temp[index, ]
   temp <- temp[!is.na(temp$z), ]
   temp <- as.image(temp)  
   I <- temp                               
   
   # Determine location of nucleus:
   nucleus <- nucleus.otolith(list(image = I))
      
   # Rotate image to orient along the major axis of the edge contour:
   if (rotate){
      temp <- as.data.frame(edges)
      D <- dist(temp)
      i <- which(as.matrix(D) == max(D), arr.ind = TRUE)   
      theta <- atan(diff(temp$y[i[1,]]) / diff(temp$x[i[1,]]))
      I <- rotate(I, -theta)
   
      # Rotate edges:
      temp <- edges
      edges$x <- temp$x * cos(-theta) - temp$y * sin(-theta)
      edges$y <- temp$x * sin(-theta) + temp$y * cos(-theta)
   
      # Rotate nucleus:
      temp <- nucleus
      nucleus[1] <- temp[1] * cos(-theta) - temp[2] * sin(-theta) 
      nucleus[2] <- temp[1] * sin(-theta) + temp[2] * cos(-theta)
   }
   
   # Fit B-spline to contour data:
   bx <- by <- bspline(75, type = "periodic")
   t <- seq(0, 1, length = length(edges$x))
   T <- evaluate(bx$basis, t)
   mx <- lm(edges$x ~ T + 0)
   my <- lm(edges$y ~ T + 0)
   coef(bx) <- coef(mx)
   coef(by) <- coef(my)
   edges <- as.list(edges)
   edges$bspline <- list(x = bx, y = by)
   
   # Create 'otolith' object:
   v <- as.otolith.default(image = I, edges = edges, nucleus = nucleus)

   return(v)
}
