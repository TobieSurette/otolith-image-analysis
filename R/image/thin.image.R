thin.image <- function(I, n){
   # THIN - Thin an 'image' object.

   # Check input aguments:
   if ((length(n) != 1) | !is.numeric(n) | (n < 1) | ((n %% 1) != 0))
      stop("'n' must be a positive integer.")

   # Choose first if there are multiple images:
   index <- which(unlist(lapply(I, class)) == "image")
   if (any(index)) d <- dim(I[[which(index)[1]]]$z) else d <- dim(I$z)
   
   # Create index vectors:
   ix <- seq(1, d[1], by = n)
   iy <- seq(1, d[2], by = n)

   return(I[ix, iy])
}
