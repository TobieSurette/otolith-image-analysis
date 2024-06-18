as.otolith.default <- function(image, edges, nucleus, species = NULL, size = 1){
   # AS.OTOLITH.DEFAULT - Convert object to an 'otolith' object.
   
   # Check input arguments:
   if ((length(size) != 1) | (!is.numeric(size)) | any(size <= 0))
      stop("'size' argument must be a positive scalar.")

   # Create an 'otolith' object:
   v <- list(image = image,
             edges = list(x = edges$x,
                          y = edges$y,
                          bspline = edges$bspline),
             nucleus = NULL,
             species = NULL,
             size = NA)
        
   class(v) <- c("otolith")
   
   if (!missing(nucleus)) nucleus(v) <- nucleus
   
   return(v)
}
