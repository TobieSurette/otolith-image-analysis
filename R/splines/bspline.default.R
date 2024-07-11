bspline.default <- function(df, knots, degree = 3, type = "regular", coefficients, intercept = TRUE, start, end, ...){
   # BSPLINE.DEFAULT - Create a 'bspline' object.

   # Parse input arguments:
   type <- match.arg(tolower(type), c("regular", "periodic"))
   
   # Define set of uniform knots:
   wrapped <- FALSE
   if (missing(knots) & !missing(df)){
      d <- (df - degree)
      step <- (1 / d)
      if (missing(start)) lower <- (-degree*step) else lower <- 0
      if (missing(end))   upper <- (1+degree*step) else upper <- 1 
      knots <- seq(lower, upper, by = step) 
      wrapped <- TRUE
   }

   # Sort knots vector:
   knots <- sort(knots)
     
   # Set start constraints:
   if (!missing(start)){
      k <- sum(knots == knots[1])
      knots <- c(rep(knots[1], degree-k + 1), knots)
   }

   # Set start constraints:
   if (!missing(end)){
      k <- sum(knots == knots[length(knots)])
      knots <- c(knots, rep(knots[length(knots)], degree-k + 1))
   }
   
   # Wrap knots if spline type is periodic:
   if ((type == "periodic") & !wrapped) knots <- c(knots, knots[length(knots)] + (knots[1:degree] - knots[1])) 
   
   # Build zeroth-degree basis functions:
   basis <- list()
   for (i in 1:(length(knots)-1)){
      basis[[i]] <- piecewise(knots[i:(i+1)], 1)
   }

   # Recusively define basis functions:
   for (i in 1:degree){
      temp <- list()
      for (j in 1:(length(basis)-1)){
         t1 <- basis[[j]]$knots
         w1 <- 1 / (t1[length(t1)] - t1[1])
         a1 <- piecewise(coef = c(-t1[1], 1))
         b1 <- w1 * a1 * basis[[j]]
         if (t1[1] == t1[length(t1)]) b1$coefficients[!is.finite(b1$coefficients)] <- 0 
         
         t2 <- basis[[j+1]]$knots
         w2 <- 1 / (t2[length(t2)] - t2[1])
         a2 <- piecewise(coef = c(t2[length(t2)], -1))
         b2 <- w2 * a2 * basis[[j+1]]
         if (t2[1] == t2[length(t2)]) b2$coefficients[!is.finite(b2$coefficients)] <- 0 
         
         b1$knots <- c(b1$knots, b2$knots[length(b2$knots)])
         b2$knots <- c(b1$knots[1], b2$knots)
         
         b1$coefficients <- rbind(b1$coefficients, rep(0, dim(b1$coefficients)[2]))
         b2$coefficients <- rbind(rep(0, dim(b2$coefficients)[2]), b2$coefficients)
         temp[[j]] <- b1 + b2
      }
      basis <- temp
   }
   
   # Create 'bspline' object:
   v <- list()

   v$basis <- as.basis(basis)
   v$knots <- knots
   v$dim <- length(basis)
   v$type <- type 
   v$degree <- degree 
   
   # Attach class identifier:
   class(v) <- c("bspline", class(v))
   
   # Update support:
   support(v) <- c(v$knots[degree+1], v$knots[length(v$knots)-degree])
   
   # Define periodic basis: 
   if (type == "periodic"){
      for (i in 1:degree){
         v$basis[[i]] <- v$basis[[i]] + v$basis[[v$dim-degree+i]]
      }
      v$basis <- as.basis(v$basis[1:(v$dim-degree)])
      v$dim <- length(v$basis)
   }
   
   # Create 'is.fixed' logical vector:
   v$is.fixed <- rep(FALSE, v$dim)
   
   # Remove first basis function 
   if (!intercept){
      v$basis <- as.basis(v$basis[2:v$dim])
      v$dim <- length(v$basis)
   }
   
   # Assign coefficient values:
   if (!missing(coefficients)) coef(v) <- coefficients else coef(v) <- rep(1, v$dim)
   
   # Fix coefficients to respect constraints:
   if (!missing(start)){
      if (!is.na(start[1])){
         v$coefficients[1] <- start[1]
         temp <- is.fixed(v)
         temp[1] <- TRUE
         is.fixed(v) <- temp
      }
      if (!is.na(start[2])){
         v$coefficients[2] <- (start[2] - v$coefficients[1] * evaluate(derive(v$basis[[1]]), knots[1])) / evaluate(derive(v$basis[[2]]), knots[1])
         temp <- is.fixed(v)
         temp[2] <- TRUE
         is.fixed(v) <- temp
      }
      if (!is.na(start[3])){
         v$coefficients[3] <- (start[3] - v$coefficients[1] * evaluate(derive(derive(v$basis[[1]])), knots[1]) - v$coefficients[2] * evaluate(derive(derive(v$basis[[2]])), knots[1])) / evaluate(derive(derive(v$basis[[3]])), knots[1])
         temp <- is.fixed(v)
         temp[3] <- TRUE
         is.fixed(v) <- temp
      } 
   }

   # Fix coefficients to respect constraints:
   if (!missing(end)){
      if (!is.na(end[1])){
         v$coefficients[v$dim] <- end[1]
         temp <- is.fixed(v)
         temp[length(temp)] <- TRUE
         is.fixed(v) <- temp
      }
      if (!is.na(end[2])){
         v$coefficients[v$dim-1] <- (end[2] - v$coefficients[v$dim] * evaluate(derive(v$basis[[v$dim]]), knots[length(knots)])) / evaluate(derive(v$basis[[v$dim-1]]), knots[length(knots)])
         temp <- is.fixed(v)
         temp[length(temp)-1] <- TRUE
         is.fixed(v) <- temp
      }
     # if (!is.na(start[3])){
     #    v$coefficients[3] <- (start[3] - v$coefficients[1] * evaluate(derive(derive(b$basis[[1]])), knots[1]) - v$coefficients[2] * evaluate(derive(derive(b$basis[[2]])), knots[1])) / evaluate(derive(derive(b$basis[[3]])), knots[1])
     #    temp <- is.fixed(v)
     #    temp[3] <- TRUE
     #    is.fixed(v) <- temp
     # } 
   }

   if (missing(coefficients)){
      if (missing(start)) start <- NULL
      if (missing(end)) end <- NULL
      if (missing(start)) index <- 1 else index <- (length(start)+1)
      if (missing(end)) index <- c(index, v$dim) else index <- c(index, v$dim - length(end))
      #xx <- c(1:max(1, index[1]-1), min(v$dim, index[2]+1):v$dim)
      #yy <- c(v$coefficients[1:max(1, index[1]-1)], v$coefficients[min(v$dim, index[2]+1):v$dim])
      #print(xx)
      #print(yy)
      #temp <- spline(xx, yy, xout = 1:v$dim)$y
      #print(temp)
      temp <- seq(v$coefficients[max(1, index[1]-1)], v$coefficients[min(v$dim, index[2]+1)], length = v$dim - length(start) - length(end) + 2)
      temp <- temp[2:(length(temp)-1)]
      v$coefficients[!v$is.fixed] <- temp
   }
      
   return(v)
}
