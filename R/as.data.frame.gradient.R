as.data.frame.gradient <- function(I, na.rm = TRUE){
   # AS.DATA.FRAME - Converts a gradient to a data frame format.

   # Create coordinate matrices:
   d <- dim(I[[1]]$z) # Image dimensions.
   x <- outer(I[[1]]$x, rep(1, d[2]))
   y <- outer(rep(1, d[1]), I[[1]]$y)

   # Contruct data frame:
   v <- data.frame(x = as.vector(x),
                   y = as.vector(y),
                   z = as.vector(I[[1]]$z),
                   w = as.vector(I[[2]]$z))

   # Rename variables:
   names(v) <- c("x", "y", names(I)[1:2])
   
   if (na.rm){
      index <- !is.na(v[, 3]) & !is.na(v[, 4])
      v <- v[index, ]
   }

   return(v)
}
