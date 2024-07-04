convolve.default <- function(I, ...){
   # CONVOLVE.DEFAULT - Default convolution operator.

   if (!is(I, "image")){
      v <- stats::convolve(I, ...)
   }else{
      v <- convolve.image(I, ...)
   }
   
   return(v)
}
