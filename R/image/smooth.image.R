smooth.image <- function(I, n = 1, method = "Gaussian"){
   # SMOOTH.IMAGE - Smooth an image.
   
   # Parse 'method' argument:
   method <- match.arg(tolower(method), c("gaussian"))
   
   # Determine Gaussian convolution opterator:
   if (method == "gaussian"){
      # Calculate distance matrix:
      d <- repvec(-n:n, ncol = 2*n+1)^2 + repvec(-n:n, nrow = 2*n+1)^2

      # Calculate weight matrix:
      m <- (1/(2*pi)) * exp(-d/2)
      m <- round(m / min(m))
      m <- m / sum(m)
   }
   
   # Apply smoothing operator:
   R <- convolve(I, m)
   
   return(R)
}
