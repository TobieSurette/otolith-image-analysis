plot.gradient <- function(I){
   # PLOT.GRADIENT - Plot a 'gradient' object.

   layout(matrix(1:2, ncol = 1))
   par(mar = c(3, 4, 2, 1))

   if (I$polar){
      plot(I$r, main = "Magnitude")
      plot(I$theta, main = "Angle")
   }else{
      plot(I$dx, main = "Horizontal gradient")
      plot(I$dy, main = "Vertical gradient")
   }
}
