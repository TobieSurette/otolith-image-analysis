as.data.frame.image <- function(I, na.rm = TRUE){
   # AS.DATA.FRAME - Converts an image to a data frame format.
   
   #X <- repvec(I$x, ncol = dim(I$z)[2])
   #Y <- repvec(I$y, nrow = dim(I$z)[1])
   
   X <- outer(I$x, rep(1, dim(I$z)[2]))
   Y <- outer(rep(1, dim(I$z)[1]), I$y)
   
   v <- data.frame(x = as.vector(X),
                   y = as.vector(Y),
                   z = as.vector(I$z))
                   
   if (na.rm){
      index <- !is.na(v$z)
      v <- v[index, ]
   }

   return(v)
}
