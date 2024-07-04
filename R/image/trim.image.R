trim.image <- function(I){
   # TRIM.IMAGE - Remove image edges with only NA values.
   
   # Find minimum and maximum rows and solumns with data:
   index <- is.na(I$z)
   ix <- which(apply(I$z, 1, function(x) !all(is.na(x))))
   iy <- which(apply(I$z, 2, function(x) !all(is.na(x))))
   
   # Define index vectors:
   ix <- min(ix):max(ix)
   iy <- min(iy):max(iy)
   
   # Extract subsets:
   I$x <- I$x[ix]
   I$y <- I$y[iy]
   I$z <- I$z[ix, iy]
   
   return(I)
}
