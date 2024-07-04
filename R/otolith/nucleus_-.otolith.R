"nucleus<-.otolith" <- function(x, value, ...){
   # NUCLEUS<- - Nucleus assignment method for an 'otolith object.
   
   # Check assignment value:
   if ((length(value) != 2) | (!is.numeric(value)))
      stop("Nucleus must be a two-element numeric vector.")
      
   x$nucleus <- list(x = value[1], y = value[2])
      
   return(x)
}
