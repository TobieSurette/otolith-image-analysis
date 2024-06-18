gradient.image <- function(I, method = "Sobel", polar = TRUE){
   # GRADIENT.IMAGE - Calculate gradient of an image.

   # Parse 'method' argument:
   method <- match.arg(tolower(method), c("simple", "sobel"))
   
   # Check input arguments:
   if ((length(polar) != 1) | !is.logical(polar)) 
      stop("'polar' must be a single logical value.")
   
   Gh <- NULL
   x <- I$z
   
   if (method == "simple"){
      print("simple")
      # Horizontal gradient:
      Gh <- x[3:dim(x)[1], 1:dim(x)[2]] - x[1:(dim(x)[1]-2), 1:dim(x)[2]]
      
      # Vertical gradient:
      Gv <- x[1:dim(x)[1], 3:dim(x)[2]] - x[1:dim(x)[1], 1:(dim(x)[2]-2)]

      # Buffer edges of gradient matrices with NA values:
      Gh <- rbind(rep(NA, dim(Gh)[2]), Gh, rep(NA, dim(Gh)[2]))       
      Gv <- cbind(rep(NA, dim(Gv)[1]), Gv, rep(NA, dim(Gv)[1]))    
   }

   if (method == "sobel"){
      # Horizontal gradient:
      Gh <- x[3:dim(x)[1], 1:(dim(x)[2]-2)] + 2*x[3:dim(x)[1], 2:(dim(x)[2]-1)] + x[3:dim(x)[1], 3:dim(x)[2]] -
            x[1:(dim(x)[1]-2), 1:(dim(x)[2]-2)] - 2*x[1:(dim(x)[1]-2), 2:(dim(x)[2]-1)] - x[1:(dim(x)[1]-2), 3:dim(x)[2]] 

      # Vertical gradient:
      Gv <- x[1:(dim(x)[1]-2), 3:dim(x)[2]] + 2*x[2:(dim(x)[1]-1), 3:dim(x)[2]] + x[3:dim(x)[1], 3:dim(x)[2]] - 
            x[1:(dim(x)[1]-2), 1:(dim(x)[2]-2)] - 2*x[2:(dim(x)[1]-1), 1:(dim(x)[2]-2)] - x[3:dim(x)[1], 1:(dim(x)[2]-2)] 
                       
      # Buffer edges of gradient matrices with NA values:
      Gh <- cbind(rep(NA, dim(Gh)[1]), Gh, rep(NA, dim(Gh)[1]))
      Gh <- rbind(rep(NA, dim(Gh)[2]), Gh, rep(NA, dim(Gh)[2]))
      Gv <- cbind(rep(NA, dim(Gv)[1]), Gv, rep(NA, dim(Gv)[1]))
      Gv <- rbind(rep(NA, dim(Gv)[2]), Gv, rep(NA, dim(Gv)[2]))    
   }

   if (is.null(Gh))
      return(NULL)
   else{
      if (polar){
         # Calculate gradient magnitude:
         G <- I
         G$z <- sqrt(Gh^2 + Gv^2)

         # Calculate gradient angle:
         theta <- I
         theta$z <- atan2(Gv, Gh)
      
         v <- list(r = G, theta = theta)
      }else{
         dx <- I
         dx$z <- Gh
         dy <- I
         dy$z <- Gv         
         v <- list(dx = dx, dy = dy)
      }
   }
   
   # Save 'polar' variable:
   v$polar <- polar
   
   # Add 'gradient' class identifier:
   class(v) <- c("gradient", class(I))
   
   return(v)
}
