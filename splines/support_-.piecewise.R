"support<-.piecewise" <- function(x, value){
   # SUPPORT<-.PIECEWISE - Set support value for a 'piecewise' object.
   
   # Check value to be assigned:
   if (length(value) != 2) stop("'piecewise' support must be a two-element vector.")
   if (!is.numeric(value)) stop("'piecewise' support must be numeric.")
   if (diff(value) <= 0) stop("'piecewise' support interval must be greater than zero.")

   # Number of intervals:
   d <- dim(x$coefficients)[1]

   # Define intervals:
   intervals <- cbind(x$knots[1:d], x$knots[2:(d+1)])

   # Eliminate intervals which fall entirely out of the support:
   index <- ((intervals[, 1] <= value[1]) & (intervals[, 2] <= value[1])) |
            ((intervals[, 1] >= value[2]) & (intervals[, 2] >= value[2]))
   intervals <- intervals[!index, , drop = FALSE]
   x$coefficients <- x$coefficients[!index, , drop = FALSE]
   
   if (length(intervals) == 0) return(NULL)
   
   # Number of intervals:
   d <- dim(x$coefficients)[1]
   
   # Truncate intervals which include a support value:
   for (i in 1:d){
      intervals[i, ] <- c(max(value[1], intervals[i,1]), min(value[2], intervals[i,2]))
   }

   # Update knot vector:
   x$knots <- sort(unique(as.vector(intervals)))
   
   return(x)
}
