"[.image" <- function(x, i, j, ...){
   # [.IMAGE - Image indexing function.

   # Choose first if there are multiple images:
   image.index <- as.numeric(which(unlist(lapply(x, class)) == "image"))
   if (length(image.index) == 0){ 
      image.index <- 1
      x <- list(x)
   }

   for (k in 1:length(image.index)){
      x[[image.index[k]]]$x <- "["(x[[image.index[k]]]$x, i, ...)
      x[[image.index[k]]]$y <- "["(x[[image.index[k]]]$y, j, ...)
      x[[image.index[k]]]$z <- "["(x[[image.index[k]]]$z, i, j, drop = FALSE, ...)
   }
   if (length(image.index) == 1) x <- x[[1]]
   
   return(x)
}
