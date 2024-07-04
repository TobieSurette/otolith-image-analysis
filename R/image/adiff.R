adiff <- function(x, y, oriented = TRUE){
   # ADIFF - Absolute angular differences.

   # If 'y' is missing, define as zero:
   if (missing(x)) y <- 0

   # Calculate angle between specified vector and 0:
   if (oriented) theta <- atan2(y, x) else theta <- atan2(y, x)

   return(theta)
}
