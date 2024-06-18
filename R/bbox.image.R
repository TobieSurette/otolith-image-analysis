bbox.image <- function(I){
   # BBOX - Image bounding box.
   
   # Choose first if there are multiple images:
   index <- unlist(lapply(I, class)) == "image"
   if (any(index)) I <- I[[which(index)[1]]]
   
   return(list(x = c(I$x[1], I$x[length(I$x)]),
               y = c(I$y[1], I$y[length(I$y)])))
}
