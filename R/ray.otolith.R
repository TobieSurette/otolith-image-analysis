ray.otolith <- function(O, angle){
   # RAY.OTOLITH - Create a B-spline ray from the nucleus to the edge of an 'otlith'.

   # Map 'angle' to [0,1]:
   angles <- (angles %% (2*pi)) / (2*pi)

   plot(O)
   for (i in 1:(length(angles)-1)){
      # Find corresponding edge coordinate and perpendicular gradient:
      xe <- evaluate(O$edges$bspline$x, angles[i])[1]
      ye <- evaluate(O$edges$bspline$y, angles[i])[1]
      dxe <- evaluate(derive(O$edges$bspline$x), angles[i])[1]
      dye <- evaluate(derive(O$edges$bspline$y), angles[i])[1]

      # Build B-spline with constrained end-points:
      xn <- O$nucleus$x
      yn <- O$nucleus$y

      dxn <- xe - xn
      dyn <- ye - yn
      #rx <- bspline(6, start = c(xn, dxn), end = c(xe, dye))
      #ry <- bspline(6, start = c(yn, dyn), end = c(ye, -dxe))
   
      rx <- bspline(8, start = c(xn), end = c(xe, dye))
      ry <- bspline(8, start = c(yn), end = c(ye, -dxe))
      
      t <- seq(0, 1, len = 1000)
      lines(evaluate(rx, t), evaluate(ry, t), lwd = 1, col = "blue")
   }
   
   # Lookup gradient values:
   G <- gradient(I, polar = FALSE)
   v <- lookup(G, evaluate(rx, t), evaluate(ry, t))
   dy <- v$dy
   dx <- v$dx
   cx <- evaluate(derive(rx), t)[, 1]
   cy <- evaluate(derive(ry), t)[, 1]
   
   # Nucleus sum-of-squares function:
   SS <- function(p){
      coef(rx) <- p[1:(length(p)/2)]
      coef(ry) <- p[((length(p)/2)+1):length(p)]

      v <- lookup(G, evaluate(rx, t), evaluate(ry, t))
      dy <- v$dy; dx <- v$dx
      cx <- evaluate(derive(rx), t)[, 1]
      cy <- evaluate(derive(ry), t)[, 1]
      
      # Calculate distance vector:
      dd <- sqrt(cx * cx + cy * cy)

      # Convert to unit vectors:
      xu <- cx / dd
      yu <- cy / dd

      # Perform projection:
      v <- (dx * xu) + (dy * yu)

      # Calculate sum-of-squares value:
      v[is.na(v)] <- sqrt(mean(v[!is.na(v)]^2, na.rm = TRUE))
      SS <- sum(v^2, na.rm = TRUE)

      return(SS)
   }
   
   p <- c(coef(rx, fixed = FALSE), coef(ry, fixed = FALSE))
   SS(p)
   p <- optim(p, SS, control = list(trace = 3))$par
   coef(rx) <- p[1:(length(p)/2)]
   coef(ry) <- p[((length(p)/2)+1):length(p)]
   cx <- evaluate(rx, t)[, 1]
   cy <- evaluate(ry, t)[, 1]
   
   lines(cx, cy, col = "red")
   
   return(list(x = rx, y = ry))
}
