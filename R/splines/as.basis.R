as.basis <- function(b){
   # AS.BASIS - Convert object to a 'basis' object.

   # Check input arguments:
   flag <- TRUE
   if (is.list(b)){
       for (i in 1:length(b)){
          if (!is(b[[i]], "piecewise")) flag <- FALSE
       }
   }
   if (!flag) stop("'basis' object must contain a list of 'piecewise' objects.")

   # Attach class identifier:
   class(b) <- c("basis", class(b))

   return(b)
}
