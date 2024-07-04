smooth.default <- function(I, ...){
   # SMOOTH.DEFAULT - Default smoothing method.

   if (!is(I, "image")){
      v <- stats::smooth(I, ...)
   }else{
      v <- smooth.image(I, ...)
   }

   return(v)
}
