quantile.image <- function(I){
   # QUANTILE.IMAGE - Balances image colour intensities using quantiles.

   a <- as.vector(I$z)
   index <- !is.na(a)
   temp <- a[index]
   temp <- match(temp, sort(temp))
   a[index] <- temp / length(temp)
   dim(a) <- dim(I$z)
   I$z <- a
   
   return(I)
}
