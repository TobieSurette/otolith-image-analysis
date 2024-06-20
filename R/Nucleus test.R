library(jpeg)

# Read image:
I <- readJPEG("photos/examples/Plaice.jpg", native = FALSE)

# Convert to greyscale:
I <- 0.2989 * I[,,1] + 0.5870 * I[,,2] +  0.1141 * I[,,3]

# Display image:
image(1:nrow(I), 1:ncol(I), I, col = grey.colors(100))

# Intensity-weighted coordinate average:
cluster <- kmeans(as.numeric(I), 2)$cluster
w <- ((cluster == 1) * as.numeric(I)) / sum(as.numeric(I)[cluster == 1], na.rm = TRUE)
x0 <- sum(w * as.numeric(repvec(1:nrow(I), ncol = ncol(I))), na.rm = TRUE)
y0 <- sum(w * as.numeric(repvec(1:ncol(I), nrow = nrow(I))), na.rm = TRUE)

points(x0, y0, pch = 21, bg = "red2")



# Filter
loglike <- function(theta, y, fixed){
   if (!missing(fixed)) theta <- c(theta, fixed)
      
   # Define index values:
   x <- 1:length(y)
   
   # Mixture proportions:
   #p <- 1 / (1 + exp(-theta["logit.p"]))
   logit.p <- theta["a"] * (x - theta["xp"]) ^2 + theta["logit.p"]
   p <- 1 / (1 + exp(-logit.p))
   
   # Mixture density:
   v <- log((1-p) * dnorm(y, theta["mu0"], exp(theta["log.sigma0"]) ) + 
                p * dnorm(y, theta["mu1"], exp(theta["log.sigma1"]) ))    
   
   return(-sum(v))
}
  
logit <- function(x) return(log(x/(1-x)))


# Sweep from point to find contour:
angle <- seq(0, pi, len = 100)[-1]
rho <- -max(c(ncol(I), nrow(I))):max(c(ncol(I), nrow(I)))
res <- NULL
for (i in 1:length(angle)){
   print(i)
   # Determine points along ray:
   xx <- unique(cbind(round(x0 + rho * cos(angle[i])), round(y0 + rho * sin(angle[i]) )))
   yy <- xx[,2]
   xx <- xx[,1]
   ix <- (xx %in% 1:nrow(I)) & (yy %in% 1:ncol(I))
   xx <- xx[ix]
   yy <- yy[ix]
   
   #lines(xx, yy)
   
   # Extract intensities in image:
   int <- rep(NA, length(xx))
   for (i in 1:length(int)) int[i] <- I[xx[i], yy[i]]
   
   theta <- c(a = 0,
              xp = which.min(abs(xx - x0)),
              logit.p = logit(sum(cluster == 2) / length(cluster)),
              mu0 = mean(I[cluster == 1]),
              mu1 = mean(I[cluster == 2]),
              log.sigma0 = log(sd(I[cluster == 1])),
              log.sigma1 = log(sd(I[cluster == 2])))
   
   # Estimate probability vector:
   fixed <- theta[c("mu0", "mu1", "log.sigma0", "log.sigma1")]
   theta <- theta[setdiff(names(theta), names(fixed))]
   loglike(theta, int, fixed = fixed)
   theta <- optim(theta, loglike, y = int, fixed = fixed)$par
   
   theta <- c(theta, fixed)
   theta <- optim(theta, loglike, y = int, control = list(maxit = 2000))$par
   theta <- optim(theta, loglike, y = int, control = list(maxit = 2000))$par
   logit.p <- theta["a"] * (x - theta["xp"]) ^2 + theta["logit.p"]
   p <- 1 / (1 + exp(-logit.p))
   ix <- c(min(which(p > 0.5)), max(which(p > 0.5)))
   
   res <- rbind(res, cbind(xx[ix], yy[ix]))
   #points(xx[ix], yy[ix], pch = 21, bg = "red2")
}

clg()
image(1:nrow(I), 1:ncol(I), I, col = grey.colors(100))
points(res[,1], res[,2], pch = 21, bg = "red2", cex = 0.5)

res <- as.data.frame(res)
names(res) <- c("x", "y")
res$angle <- as.numeric(rbind(angle, angle + pi))
res <- res[order(res$angle), ]

# Show otolith shape:
plot(range(res$x), range(res$y), type = "n")
polygon(res$x, res$y)

library(FluMoDL)
library(splines)

plot(range(res$x), range(res$y), type = "n", xlab = "", ylab = "")
ix <- seq(1, nrow(I), by = 2)
iy <- seq(1, ncol(I), by = 2)
image((1:nrow(I))[ix], (1:ncol(I))[iy], I[ix,iy], col = grey.colors(100), add = TRUE)
polygon(res$x, res$y, lwd = 2)
for (df in c(50)){
   X <- pbs(res$angle, df = df, Boundary.knots = c(0, 2*pi))
   mx <- lm(res$x ~ X)
   my <- lm(res$y ~ X)
   
   lines(predict(mx), predict(my), lwd = 1, col = "red2")
}
mtext("X", 1, 2.5, cex = 1.25, font = 2)
mtext("Y", 2, 2.5, cex = 1.25, font = 2)
box(col = "grey50")
 
library(sf)
pts = st_sfc(st_point(c(400, 400)), st_point(c(800, 600)), st_point(c(700, 600)))
pol = st_polygon(list(cbind(c(res$x, res$x[1]), c(res$y, res$y[1]))))
inside <- st_within(pts, pol)
unlist(lapply(inside, function(x) ifelse(length(x) == 0, 0, 1)))

point.in.polygon(800, 600, res$x, res$y)

plot(res$angle, predict(mx))
lines(res$angle, res$x, lwd = 2, col = "red2")

plot(res$angle, predict(my))
lines(res$angle, res$y, lwd = 2, col = "red2")


clg()
ix <- seq(1, nrow(I), by = 2)
iy <- seq(1, ncol(I), by = 2)
image((1:nrow(I))[ix], (1:ncol(I))[iy], I[ix,iy], col = grey.colors(100))
lines(xx, yy)
points(x0, y0, pch = 21, bg = "red2", cex = 1.25)
points(xx[ix[1]], yy[ix[1]], pch = 21, bg = "red2")
points(xx[ix[2]], yy[ix[2]], pch = 21, bg = "red2")

