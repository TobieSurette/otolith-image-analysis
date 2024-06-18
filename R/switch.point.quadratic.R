switch.point.quadratic <- function(y, x, theta, p = 0.5, tol = 1.0E-4){
   # SWITCH.POINT - Univariate binomial regression to find the switch point along a coordinate.

   # Define design matrix if missing:
   if (missing(x)){
      x <- cbind(rep(1, length(y)), 1:length(y), (1:length(y))^2)
   }
   
   # Define initial parameter values:
   if (missing(theta)){
      theta <- solve(t(x) %*% x) %*% (t(x) %*% (y-0.5))
   }

   # Apply Newton's method:
   xp <- c(1, 1)
   xpold <- xp + 2*tol
   while (sum(abs(xp-xpold)) > tol){
      xpold <- xp
      a <- exp(x %*% theta)
      pr <- a / (1+a)
      d <- y - pr
      G <- matrix(NA, ncol = 1, nrow = 3)
      G[1] <- sum(x[,1]*d)
      G[2] <- sum(x[,2]*d)
      G[3] <- sum(x[,3]*d)
      b <- pr * (1-pr)
      H <- matrix(NA, ncol = 3, nrow = 3)
      H[1,1] <- sum(b)
      H[2,2] <- sum((x[,2]^2)*b)
      H[3,3] <- sum((x[,3]^2)*b)
      H[1,2] <- H[2,1] <- sum(x[,2]*b)
      H[1,3] <- H[3,1] <- sum(x[,3]*b)
      H[2,3] <- H[3,2] <- sum(x[,2]*x[,3]*b)
      
      theta <- theta + solve(H, G)
   }

   # Find quadratic roots:
   xp <- Re(polyroot(theta))

   return(xp)
}