rotate.image <- function(I, angle, x = 0, y = 0){
   # ROTATE.IMAGE - Rotate an image by specified degrees.

   # Create rotation matrix:
   if (!is.numeric(angle) | length(angle) != 1)
      stop("'angle' must be a numeric scalar.")
   R <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), nrow = 2)

   # Center image about rotation point:
   I$x <- I$x - x
   I$y <- I$y - y

   # Extract image corner coordinates:
   corners <- cbind(c(I$x[1], I$x[length(I$x)], I$x[length(I$x)], I$x[1]),
                    c(I$y[1], I$y[1], I$y[length(I$y)], I$y[length(I$y)]))
   rc <- t(R %*% t(corners))
   
   # Prepare target image matrix:
   dx <- (mean(unique(diff(I$x))) * abs(cos(angle)) + mean(unique(diff(I$y))) * abs(sin(angle))) / 2
   dy <- (mean(unique(diff(I$y))) * abs(cos(angle)) + mean(unique(diff(I$x))) * abs(sin(angle))) / 2
   dmin <- apply(rc, 2, min)
   dmax <- apply(rc, 2, max)
   xi <- seq(dmin[1], dmax[1], by = dx)
   yi <- seq(dmin[2], dmax[2], by = dy)
   xx <- repvec(xi, ncol = length(yi))
   yy <- repvec(yi, nrow = length(xi))

   # Vectorized target image:
   vv <- cbind(as.vector(xx), as.vector(yy))

   # Perform inverse rotation operator:
   vi <- t(solve(R) %*% t(vv))

   # Remove coordinates which fall outside of the original image dimensions:
   index <- (vi[,1] >= I$x[1]) & (vi[,1] <= I$x[length(I$x)]) & (vi[,2] >= I$y[1]) & (vi[,2] <= I$y[length(I$y)])
   vii <- vi[index, ]

   # Convert to matrix coordinates:
   dx <- abs(mean(unique(diff(I$x))))
   dy <- abs(mean(unique(diff(I$y))))
   wx <- (I$x[length(I$x)] - I$x[1])
   wy <- (I$y[length(I$y)] - I$y[1])
   vx <- floor((length(I$x)-1)*((vii[, 1] - I$x[1]) / wx)) + 1
   vy <- floor((length(I$y)-1)*((vii[, 2] - I$y[1]) / wy)) + 1
   
   # Lookup values in original image:
   zz <- rep(NA, dim(vv)[1])
   zz[index] <- t(I$z)[(vx-1) * dim(I$z)[2] + vy]
   dim(zz) <- dim(xx)

   # print(dim(zz))
   # print(length(xi))
   # print(length(yi))

   v <- as.image(zz, x = xi, y = yi, flip = FALSE)

   return(v)
}
