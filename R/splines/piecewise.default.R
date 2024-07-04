piecewise.default <- function(knots, coefficients){
   # PIECEWISE - Creates a 'piecewise' polynomial object.
   
   # Check 'knot' input argument:
   if (!missing(knots)){
      if (any(is.na(knots))) stop("'knots' must not contain NA values.")
      if (!is.null(dim(knots))) stop("'knots' must be a vector.")
      if (!is.numeric(knots)) stop("'knots' must be numeric.")
      #knots <- sort(knots)
   }else{
      knots <- c(-Inf, Inf)
   }
   
   # Check 'coefficient' input argument:
   if (!is.null(dim(coefficients))){
      if (!any(dim(coefficients) == 1) & is.null(knots))
         stop("'coefficients' must be a vector for unspecified 'knots' argument.")
   }else{
      coefficients <- matrix(coefficients, nrow = 1)
   }
   if (any(is.na(coefficients))) stop("'coefficients' must not contain NA values.")
   
   # Check consistency of input arguments:
   if (dim(coefficients)[1] != (length(knots)-1))
      stop("Number of row vectors in 'coefficient' must be one less the number of elements in 'knots'.")

   p <- list(knots = knots, coefficients = coefficients)
   class(p) <- c("piecewise", class(p))

   return(p)
}

