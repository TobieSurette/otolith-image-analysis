"+.piecewise" <- function(p, x){
   # Addition operator for 'piecewise' objects.

   # Exchange arguments:
   if (("piecewise" %in% class(x)) & !("piecewise" %in% class(p))){
      temp <- x
      x <- p
      p <- temp
   }
   
   # Increment intercept coefficient:
   if (is.numeric(x) & (length(x) == 1)){
      p$coefficients[,1] <- p$coefficients[,1] + x
   }

   # Add two piecewise objects together:
   if (("piecewise" %in% class(p)) & ("piecewise" %in% class(x))){
      # Set degree to maximum degree:
      deg <- max(degree(p), degree(x))
      degree(p) <- deg
      degree(x) <- deg
      
      # Knots are identical - simply add coefficients:
      if (all(p$knots == x$knots)){
         p$coefficients <- p$coefficients + x$coefficients
         return(p)
      }
      
      # Add with possibly irregular set of intervals:
      knots <- unique(sort(c(p$knots, x$knots)))
      intervals <- cbind(knots[1:(length(knots)-1)], knots[2:length(knots)])
      Ip <- cbind(p$knots[1:(length(p$knots)-1)], p$knots[2:length(p$knots)])
      Ix <- cbind(x$knots[1:(length(x$knots)-1)], x$knots[2:length(x$knots)])
      coefficients <- matrix(0, nrow = dim(intervals)[1], ncol = deg + 1)
      for (i in 1:dim(intervals)[1]){
         # If interval lies entirely within 'p' or 'x' interval:
         index <- which((intervals[i,1] >= Ip[, 1]) & (intervals[i,2] <= Ip[, 2]))
         if (length(index) == 1){
            coefficients[i, ] <- coefficients[i, ] + p$coefficients[index, ]
         }
         index <- which((intervals[i,1] >= Ix[, 1]) & (intervals[i,2] <= Ix[, 2]))
         if (length(index) == 1){
            coefficients[i, ] <- coefficients[i, ] + x$coefficients[index, ]         
         }
      }

      # Create result object:
      p <- piecewise(knots = knots, coefficients = coefficients)  
   }

   return(p)
}
