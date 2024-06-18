switch.point <- function(y, x, n, theta, alpha = 0.5, tol = 1.0E-4){
   # SWITCH.POINT - Univariate binomial regression to find the switch point along a coordinate.

   # Define design matrix if missing:
   if (missing(x)) x <- 1:length(y)

   # Check 'alpha':
   if ((length(alpha) != 1) | (alpha < 0) | (alpha > 1)) 
      stop("'alpha' must be a numeric scalar between zero and one.")
   
   # Define initial parameter values:
   if (missing(theta)){
      X <- cbind(matrix(rep(1, length(x)), ncol = 1), x)
      theta <- solve(t(X) %*% X) %*% (t(X) %*% (y-0.5))
   }

   # # Define missing 'n' vector:
   if (missing(n)){
      n <- rep(1, length(y))
   }
   
   # Apply Newton's method:
   xp <- 1
   xpold <- xp + 2*tol
   while (abs(xp-xpold) > tol){
      xpold <- xp
      a <- exp(theta[1] + (theta[2]*x))
      pr <- a / (1+a)
      d <- y - n*pr
      G <- matrix(NA, ncol = 1, nrow = 2)
      G[1] <- sum(d)
      G[2] <- sum(x*d)
      b <- pr * (1-pr)
      H <- matrix(NA, ncol = 2, nrow = 2)
      H[1,1] <- sum(n*b)
      H[1,2] <- H[2,1] <- sum(n*x*b)
      H[2,2] <- sum(n*(x^2)*b)
      theta <- theta + solve(H, G)
      xp <- (log(alpha/(1-alpha))- theta[1]) / theta[2]
      if (is.na(xp)) xp <- xpold
   }
   
   return(xp)
}
