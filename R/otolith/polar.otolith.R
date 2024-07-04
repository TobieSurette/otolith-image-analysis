polar.otolith <- function(O, ...){
   # POLAR.OTOLITH - Polar transformation for an 'otolith' object.

   P <- polar.image(O$image, x = O$nucleus$x, y = O$nucleus$y, edges = O$edges[c("x", "y")], ...)

   return(P)
}
