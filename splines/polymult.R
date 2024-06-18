polymult <- function(x, y){
   # POLYMULT - Polynomial multiplication.
   
   # Check if one of the arguments is zero:
   if (all(x == 0) | all(y == 0)) return(0)
   
   # Build index matrix:
   ix <- repvec(0:(length(x)-1), nrow = length(y))
   iy <- repvec(0:(length(y)-1), ncol = length(x))
   index <- as.vector(ix + iy + 1)
   
   # Build pairwise coefficient matrix:
   xx <- repvec(x, nrow = length(y))
   yy <- repvec(y, ncol = length(x))
   zz <- xx * yy
   z <- as.vector(zz)

   # Remove zero coefficients:
   index <- index[z != 0]
   z <- z[z != 0]

   # Sum up coefficient values with same exponent:
   temp <- aggregate(z, by = list(index), sum)
   
   # Build result vector:
   r <- rep(0, max(temp[, 1]))
   r[temp[, 1]] <- temp[, 2]
   
   return(r)
}
