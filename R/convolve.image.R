convolve.image <- function(I, x, pad = TRUE){
   # CONVOLVE - Discrete image convolution opterator.
   
   # Check input arguments:
   if (!is.matrix(x)) stop("'x' must be a matrix.")
   if (any(dim(x) %% 2 == 0)) stop("'x' have odd-numbered dimensions.")
   
   # Initialize result matrix:
   R <- matrix(0, nrow = dim(I$z)[1]-(dim(x)[1]-1), ncol = dim(I$z)[2]-(dim(x)[2]-1))

   # Perform convolution:
   for (i in 1:dim(x)[1]){
      ix <- i:(dim(I$z)[1]-dim(x)[1]+i)
      for (j in 1:dim(x)[2]){
         if (x[i,j] != 0){
            iy <- j:(dim(I$z)[2]-dim(x)[2]+j)     
            R <- R + x[i,j] * I$z[ix, iy]
         }
      }
   }
   
   # Pad matrix edges with NA values:
   if (pad){
      temp <- matrix(NA, nrow = (dim(x)[1]-1)/2, ncol = dim(R)[2])
      R <- rbind(temp, R, temp)
      temp <- matrix(NA, ncol = (dim(x)[2]-1)/2, nrow = dim(R)[1])
      R <- cbind(temp, R, temp)
   }

   # Update image:
   I$z <- R
   
   return(I)
}
