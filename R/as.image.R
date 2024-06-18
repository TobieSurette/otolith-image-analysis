as.image <- function(z, x, y, xlim, ylim, fun = mean, flip = TRUE, na.rm = TRUE, ...){
   # AS.IMAGE - Create an 'image' object.
   
   # Check input arguments:
   if (!missing(xlim)){
      if (length(xlim) != 2) stop("'xlim' must be a two-element numeric vector.")
      xlim <- sort(xlim)
      if (diff(xlim) == 0) stop("'xlim' must have different values.")
   }
   if (!missing(ylim)){
      if (length(ylim) != 2) stop("'ylim' must be a two-element numeric vector.")
      ylim <- sort(ylim)
      if (diff(ylim) == 0) stop("'ylim' must have different values.")
   }
   if (!missing(x)){
      if (length(dim(x)) != 0) stop("'x' must be a vector.")
      if ((length(unique(diff(x))) != 1) & any(diff(x) <= 0))
         stop("'x' must be a constantly increasing vector.")
   }
   if (!missing(y)){
      if (length(dim(y)) != 0) stop("'y' must be a vector.")
      if ((length(unique(diff(y))) != 1) & any(diff(y) <= 0))
         stop("'y' must be a constantly increasing vector.")
   }
    
   # Convert a set of 'x', 'y' and 'z' values to an image:
   if (is.data.frame(z)){
      if (all(c("x", "y", "z") %in% tolower(names(z)))){
         names(z) <- tolower(names(z))
         
         # Remove irrelevant coordinate data:
         z <- z[!is.na(z$x) & !is.na(z$y) & is.finite(z$x) & is.finite(z$y), ]
         if (!missing(xlim)) z <- z[(z$x >= xlim[1]) & (z$x <= xlim[2]), ]
         if (!missing(ylim)) z <- z[(z$y >= ylim[1]) & (z$y <= ylim[2]), ] 
         
         # Define image 'x' vector:
         ux <- sort(unique(z$x))
         if (missing(xlim)) xlim <- c(min(ux), max(ux))
         dx <- min(diff(ux))
         if ((diff(xlim) / dx) > 5000){
            temp <- table(diff(ux))
            dx <- as.numeric(names(temp[temp == max(temp)])[1])
         } 
         if ((diff(xlim) / dx) > 5000) dx <- (diff(xlim) / 999) 
         ix <- seq(xlim[1], xlim[2], by = dx)
         
         # Define image 'y' vector:
         uy <- sort(unique(z$y))
         if (missing(ylim)) ylim <- c(min(uy), max(uy))
         dy <- min(diff(uy))
         if ((diff(ylim) / dy) > 5000){
            temp <- table(diff(uy))
            dy <- as.numeric(names(temp[temp == max(temp)])[1])
         } 
         if ((diff(ylim) / dy) > 5000) dy <- (diff(ylim) / 999)  
         iy <- seq(ylim[1], ylim[2], by = dy)
         
         # Create two-dimensional coordinate and value matrices:
         xx <- as.vector(repvec(ix, ncol = length(iy)))
         yy <- as.vector(repvec(iy, nrow = length(ix)))
         zz <- rep(NA, length(xx))
         
         # Calculate coordinate index vectors and grouped 'z' values:
         index.x <- floor((z$x - xlim[1]) / dx) + 1
         index.y <- floor((z$y - ylim[1]) / dy) + 1
         if (length(formals(fun)) > 1){
            res <- aggregate(list(z = z$z), by = list(x = index.x, y = index.y), fun, na.rm = na.rm)
         }else{
            res <- aggregate(list(z = z$z), by = list(x = index.x, y = index.y), fun)
         }
         
         # Assign grouped 'z' values:
         zz[(res$y-1) * length(ix) + res$x] <- res$z
         dim(zz) <- c(length(ix), length(iy))
         
         # Prepare 'image' object:
         r <- list(x = ix, y = iy, z = zz)
         class(r) <- "image"
         
         return(r)
      }
   }
   
   # Check input arguments:
   if (!((is.numeric(z) & (length(dim(z)) >= 2)))) stop("'x' must be a numeric matrix or array.")
   
   # Convert RGB image to black and white:
   if (length(dim(z)) == 3) z <- bw.image(z)

   # Check that 'x' is now a numeric matrix:
   if (length(dim(z)) != 2) stop("'z' has inconsistent dimensions.")
   
   # Define 'x' and 'y' if they are not specified:
   if (missing(x)){
      if (!missing(xlim)) x <- seq(xlim[1], xlim[2], len = dim(z)[1]) else x <- 1:dim(z)[1]
   }
   if (missing(y)){
      if (!missing(ylim)) y <- seq(ylim[1], ylim[2], len = dim(z)[1]) else y <- 1:dim(z)[2]
   }

   # Flip image so that 'x' and 'y' are the in the first and second dimensions:
   if (flip){
      z <- t(z)
      z <- z[, dim(z)[2]:1]
      temp <- x
      x <- y
      y <- temp
   }
   
   # Truncate image if both 'x' or 'y' and limits are specified:
   if (!missing(xlim) & !missing(x)){
      index <- (x >= xlim[1]) & (x <= xlim[2])
      x <- x[index]
      z <- z[index, ]
   }
   if (!missing(ylim) & !missing(y)){
      index <- (y >= ylim[1]) & (y <= ylim[2])
      y <- y[index]
      z <- z[, index]
   }   
   
   # Prepare 'image' object:
   r <- list(x = x, y = y, z = z)
   class(r) <- "image"
   
   return(r)
}
