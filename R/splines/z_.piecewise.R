"*.piecewise" <- function(x, y){
   # Addition operator for 'piecewise' objects.

   # Exchange arguments if 'x' is numeric:
   if (("piecewise" %in% class(y)) & !("piecewise" %in% class(x))){
      temp <- y
      y <- x
      x <- temp
   }

   # Multiply coefficients:
   if (is.numeric(y) & (length(y) == 1)){
      x$coefficients <- x$coefficients * y
   }

   # Multiply two piecewise objects together:
   if (("piecewise" %in% class(x)) & ("piecewise" %in% class(y))){
      # Exchange so that 'y' is of length 2:
      if ((length(x$knots) == 2) | (all(x$knots == c(-Inf, Inf)))){
         temp <- x
         x <- y
         y <- temp
      }
      
      # Multiplication by a single polynomial:
      if (all(y$knots == c(-Inf, Inf))){
         y$coefficients <- repvec(y$coefficients, nrow = dim(x$coefficients)[1])
         y$knots <- x$knots
      } 

      # Piecewise polynomial multiplication:
      if (all(x$knots == y$knots)){
         temp <- list()
         for (i in 1:dim(x$coefficients)[1]){
            temp[[i]] <- polymult(x$coefficients[i, ], y$coefficients[i, ])
         }
         n <- max(unlist(lapply(temp, length)))
         x$coefficients <- matrix(0, nrow = length(temp), ncol = n)
         for (i in 1:length(temp)){
            x$coefficients[i, 1:length(temp[[i]])] <- temp[[i]]
         }
      }
   }

   return(x)
}