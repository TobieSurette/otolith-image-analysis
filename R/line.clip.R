line.clip <- function(x, y, xbox, ybox){
   # LINE.CLIP - Clip line to a box.

   inbox <- function(x, y, xbox, ybox){
      # INBOX - Determine whether points lies inside a box.
      index <- (x >= xbox[1]) & (x <= xbox[2]) & (y >= ybox[1]) & (y <= ybox[2])
      return(index)
   }

   # If points are inside, then return them:
   inside <- inbox(x, y, xbox, ybox)
   if (all(inside)) return(list(x = x, y = y))

   # Calculate slope and intercept:
   slope <- diff(y) / diff(x)
   intercept <- y[1] - slope * x[1]

   # Slope is infinite:
   if (!is.finite(slope)){
      if ((x[1] >= xbox[1]) & (x[1] <= xbox[2])){
         return(list(x = x,
                     y = c(max(c(min(y), min(ybox))), min(c(max(y), max(ybox))))))
      }else{
         return(NULL)
      }
   }

   # Slope is zero:
   if (slope == 0){
      print(1)
      if ((y[1] >= ybox[1]) & (y[1] <= ybox[2])){
         return(list(x = c(max(c(min(x), min(xbox))), min(c(max(x), max(xbox)))),
                     y = y))
      }else{
         return(NULL)
      }
   }

   # Intersection points with the four box boundaries:
   p <- c(xbox[1], slope * xbox[1] + intercept)           # First 'x'.
   p <- rbind(p, c(xbox[2], slope * xbox[2] + intercept)) # Last  'x'.
   p <- rbind(p, c((ybox[1]-intercept)/slope, ybox[1]))   # First 'y'.
   p <- rbind(p, c((ybox[2]-intercept)/slope, ybox[2]))   # Last  'y'.

   xp <- p[, 1]; yp <- p[, 2]
   index <- which(inbox(xp, yp, xbox, ybox))
   xp <- xp[index]
   yp <- yp[index]

   # Return the two-points at the intersections:
   if (!any(inside)){
      # Define start and end coordinates:
      xp <- xp[1:2]
      yp <- yp[1:2]
      return(list(x = xp, y = yp))
   }

   # Calculate t value:
   xt <- c(x[inside], x[!inside])
   yt <- c(y[inside], y[!inside])
   #t <- (xt[2] - xp) / (xt[2] - xt[1])
   t <- (yt[2] - yp) / (yt[2] - yt[1])
   index <- (t >= 0) & (t <= 1)
   x[!inside] <- xp[index]
   y[!inside] <- yp[index]

   return(list(x = x, y = y))
}
