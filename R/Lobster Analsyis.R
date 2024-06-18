setwd("C:/Users/SuretteTJ/Desktop/Ageing")

library(jpeg)
library(gulf)

source("source.R")
graphics.off()

I <- readJPEG("Images/Lobster.jpg", native = FALSE)
I <- as.image(I)
I$z[I$z > 0.9] <- NA
I <- quantile(I)
plot(I, adjust = FALSE)
G <- gradient(I)
G$G <- quantile(G$G)

adiff <- function(x, y){
   x <- (x %% (2*pi)) - pi
   y <- (y %% (2*pi)) - pi

   index <- !is.na(x) & (x > (pi/2))
   x[index] <- x[index] - pi
   index <- !is.na(x) & (x < -(pi/2))
   x[index] <- x[index] + pi

   index <- !is.na(y) & (y > (pi/2))
   y[index] <- y[index] - pi
   index <- !is.na(y) & (y < -(pi/2))
   y[index] <- y[index] + pi

   r <- apply(abs(x - cbind(y - pi, theta, theta + pi)), 1, min)

   return(r)
}


k <- 5000
res <- matrix(NA, nrow = k, ncol = 5)
for (i in 1:k){
   if ((i %% 100) == 0) print(i)
   flag <- TRUE
   while (flag){
      # Pick two random points:
      xr <- floor(runif(2)*dim(I$z)[1] + 1)
      yr <- floor(runif(2)*dim(I$z)[2] + 1)

      # Calculate slope and intercept:
      m <- diff(yr) / diff(xr)
      b <- yr[1] - m * xr[1]

      # Normalize vector:
      xr <- cos(atan(m))
      yr <- sin(atan(m))
      
      if (is.finite(m)) flag <- FALSE
   }
   
   r <- lookup(G$G, slope = m, intercept = b)[, 3]
   theta <- lookup(G$theta, slope = m, intercept = b)[, 3]
   x <- r * cos(theta)
   y <- r * sin(theta)

   d <- (xr * x) + (yr * y)

   res[i, ] <- c(m , b, mean(d, na.rm = TRUE), sd(d, na.rm = TRUE), sum(!is.na(d)))
}

res <- res[order(res[, 4]), ]
windows()
plot(G$, adjust = FALSE)
for (i in 500:1000){
   abline(res[i, 2], res[i, 1], col = "red")
}

r <- lookup(G$theta, x = c(0, 1600), y = c(1000, 1000))
plot(G$theta)
points(r[, 1], r[, 2], pch = 21, bg = "red")
plot(r[, 3])

plot(apply(G$theta$z[1:1600, 550:650], 1, mean))
