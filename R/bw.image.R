bw.image <- function(x){
   # BW - Convert an image to black and white.

   if (dim(x)[3] == 3){
      w <- c(0.2989, 0.5870, 0.1141)
      x <- w[1]*x[,,1] + w[2]*x[,,2] + w[3]*x[,,3]
   }
   
   return(x)
}
