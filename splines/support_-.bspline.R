"support<-.bspline" <- function(x, value){
   # SUPPORT<-.BSPLINE - Set support value for a 'bspline' object.

   # Truncate basis functions:
   basis <- list()
   knots <- NULL
   coefficients <- NULL
   k <- 1
   for (i in 1:x$dim){
      temp <- x$basis[[i]]
      support(temp) <- value
      if (!is.null(temp)){
         basis[[k]] <- temp
         knots <- c(knots, temp$knots)
         coefficients <- c(coefficients, x$coefficients[i])
         k <- k + 1
      }
   }

   # Update 'bspline' object:
   if (length(basis) > 0){
      x$basis <- as.basis(basis)
      x$dim <- length(basis)
      x$coefficients <- coefficients
      x$knots <- sort(unique(knots))
      return(x)
   }else{
      return(NULL)
   }
}
