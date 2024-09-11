library(jpeg)

# Read image:
I <- readJPEG("photos/examples/Plaice.jpg", native = FALSE)

# Convert to greyscale:
I <- 0.2989 * I[,,1] + 0.5870 * I[,,2] +  0.1141 * I[,,3]
I <- I > 0.25

# Display image:
#image(1:nrow(I), 1:ncol(I), I, col = grey.colors(100))
ix <- seq(1, nrow(I), by = 3)
iy <- seq(1, ncol(I), by = 3)
image((1:nrow(I))[ix], (1:ncol(I))[iy], (I[ix, iy] - mean(I[ix, iy])), col = grey.colors(100))

L <- cluster.image(I)

# Remove small blobs:
t <- table(L)
L[which(L != as.numeric(names(t[t == max(t)])))] <- NA

# Find pixels making up the contours:
contour <- NULL
for (i in 2:(nrow(L)-1)){
   print(i)
   for (j in 2:(ncol(L)-1)){
       if (!is.na(L[i,j]) & any(is.na(L[(i-1):(i+1),(j-1):(j+1)]))){
          contour <- rbind(contour, c(i, j))
       }
   }
}

#
M <- matrix(0, nrow = nrow(L), ncol = ncol(L))
M[contour] <- 1
M <- cluster.image(M)
M[which(M != 1)] <- NA
contour <- which(M==1, arr.ind = TRUE)

# Read image:
I <- readJPEG("photos/examples/Plaice.jpg", native = FALSE)

# Convert to greyscale:
I <- 0.2989 * I[,,1] + 0.5870 * I[,,2] +  0.1141 * I[,,3]

image(1:nrow(I), 1:ncol(I), I, col = grey.colors(100), xlim = c(200, 1200), ylim = c(350, 1300))
points(contour[,1], contour[,2], cex = 0.25, col = "green3")

# Calculate center of gravity:
xx <- repvec(1:nrow(I), ncol = ncol(I))
yy <- repvec(1:ncol(I), nrow = nrow(I))
mx <- sum(L * I * xx, na.rm = TRUE) / sum(L * I, na.rm = TRUE)
my <- sum(L * I * yy, na.rm = TRUE) / sum(L * I, na.rm = TRUE)
points(mx, my, pch = 21, bg = "red2")

# Choose a set of 4 contour points:
ix <- round(seq(1,length(theta), len = 32))
points(contour[ix,1], contour[ix,2], pch = 21, bg = "red2")

t <- seq(0, 1, len = 100)
dx <- dy <- dz <- NULL
for (i in 1:length(ix)){
   dx <- rbind(dx, mx + t * (contour[ix[i], 1] - mx))
   dy <- rbind(dy, my + t * (contour[ix[i], 2] - my))
   dz <- rbind(dz, I[round(cbind(dx[i,], dy[i,]))])
}

plot(c(0, length(t)), c(0,1), type = "n")
for (i in 1:length(ix)){
   lines(1:length(t), dz[i,])
}


# Order contour points by angle wrt center of mass:
theta <- atan2(my - contour[,2], mx - contour[,1])
contour <- contour[order(theta),]
theta   <- theta[order(theta)]

# Calculate normal vectors for each contour coordinate:
k <- 10
for (i in 1:nrow(contour)){
   ix <- (i-10):(i+10)
   ix[ix < 1] <- length(theta) + ix[ix < 1]
   ix[ix > length(theta)] <- ix[ix > length(theta)] - length(theta)
   xx <- 1:length(ix)
   yy <- contour[ix,1]
   modelx <- lm(yy ~ xx)
   xx <- 1:length(ix)
   zz <- contour[ix,2]
   modely <- lm(zz ~ xx)   
   
   m <- diff(predict(modely)[c(length(ix), 1)]) / diff(predict(modelx)[c(length(ix), 1)])
   
   pm <- -1/m
   
   b <- contour[i,2] - pm * contour[i,1]
   abline(b, pm, lwd = 0.1, col = fade("black", 0.1))
   
   print(c(pm, b))
}

# Image gradient:
Gx <- I[3:nrow(I), 2:(ncol(I)-1)] - I[1:(nrow(I)-2), 2:(ncol(I)-1)]
Gy <- I[2:(nrow(I)-1), 3:ncol(I)] - I[2:(nrow(I)-1), 1:(ncol(I)-2)]
Gx <- cbind(NA, rbind(NA, Gx, NA), NA)
Gy <- cbind(NA, rbind(NA, Gy, NA), NA)
Mag <- sqrt(Gx*Gx + Gy*Gy) 
eta <- atan2(Gy,Gx)

image(1:nrow(I), 1:ncol(I), I, col = grey.colors(100), xlim = c(200, 1200), ylim = c(350, 1300))
points(contour[,1], contour[,2], cex = 0.25, col = "green3")

# Calculate center of gravity:
xx <- repvec(1:nrow(I), ncol = ncol(I))
yy <- repvec(1:ncol(I), nrow = nrow(I))
mx <- sum(L * I * xx, na.rm = TRUE) / sum(L * I, na.rm = TRUE)
my <- sum(L * I * yy, na.rm = TRUE) / sum(L * I, na.rm = TRUE)
points(mx, my, pch = 21, bg = "red2")

library(splines)


B %*% c(1,1,1,1,1,1,1,1)

# Interpolate a curve between two coordinates:
xy <- c(0, 0)
Mag * cos(theta)

dx


t <- seq(0, 1, len = 101)
B <- bs(t, knots = c(0, 0, seq(0, 1, len = 3), 1, 1))

loglike <- function(beta){
   xx <- x0 + t * (x1-x0) + B[, 4:6] %*% beta[1:3]
   yy <- y0 + t * (y1-y0) + B[, 4:6] %*% beta[4:6]
   
   ll <- Mag[round(cbind(xx, yy))] * abs(cos(atan2(yy, xx) - eta[round(cbind(xx, yy))]))

   return(-sum(ll))
}
   
beta <- 20*rnorm(6)



image(1:nrow(I), 1:ncol(I), I, col = grey.colors(100), xlim = c(200, 1200), ylim = c(350, 1300))
points(contour[,1], contour[,2], cex = 0.25, col = "green3")

# Calculate center of gravity:
xx <- repvec(1:nrow(I), ncol = ncol(I))
yy <- repvec(1:ncol(I), nrow = nrow(I))
mx <- sum(L * I * xx, na.rm = TRUE) / sum(L * I, na.rm = TRUE)
my <- sum(L * I * yy, na.rm = TRUE) / sum(L * I, na.rm = TRUE)
points(mx, my, pch = 21, bg = "red2")

for (i in 1:nrow(dx)){
   print(i)
   x0 <- dx[i,1]
   x1 <- dx[i,ncol(dx)]
   y0 <- dy[i,1]
   y1 <- dy[i,ncol(dy)]
   
   beta <- 25 * rnorm(6)
   beta <- optim(beta, loglike, control = list(trace = 3))$par
   beta <- optim(beta, loglike, control = list(trace = 3))$par
   beta <- optim(beta, loglike, control = list(trace = 3))$par
   
   xx <- x0 + t * (x1-x0) + B[, 4:6] %*% beta[1:3]
   yy <- y0 + t * (y1-y0) + B[, 4:6] %*% beta[4:6]
   lines(x0 + t * (x1-x0), y0 + t * (y1-y0), lwd = 1, col = fade("red3"), lty = "dashed")
   lines(xx, yy, lwd = 2, col = "red3")
}

loglike <- function(phi, k = 8){
   mx <- phi[1]
   my <- phi[2]

   ix <- round(seq(1, nrow(contour), len = k))
   
   t <- seq(0, 1, len = 100)
   dx <- dy <- NULL
   for (i in 1:length(ix)){
      dx <- rbind(dx, mx + t * (contour[ix[i], 1] - mx))
      dy <- rbind(dy, my + t * (contour[ix[i], 2] - my))
   }

   
   ll <- NULL
   for (i in 1:nrow(dx)){
      xx <- dx[i,]
      yy <- dy[i,]
      ll <- rbind(ll, rev(t) * Mag[round(cbind(xx, yy))] * abs(cos(atan2(yy, xx) - eta[round(cbind(xx, yy))])))
   }
  
   points(mx, my)
   
   return(-sum(ll))   
}

clg()
image(1:nrow(I), 1:ncol(I), I, col = grey.colors(100), xlim = c(200, 1200), ylim = c(350, 1300))
points(contour[,1], contour[,2], cex = 0.25, col = "green3")

# Calculate center of gravity:
xx <- repvec(1:nrow(I), ncol = ncol(I))
yy <- repvec(1:ncol(I), nrow = nrow(I))
mx <- sum(L * I * xx, na.rm = TRUE) / sum(L * I, na.rm = TRUE)
my <- sum(L * I * yy, na.rm = TRUE) / sum(L * I, na.rm = TRUE)
points(mx, my, pch = 21, bg = "red2")

phi <- optim(c(mx, my), loglike, k = 24, control = list(trace = 3))$par
phi <- optim(phi, loglike, k = 24, control = list(trace = 3))$par
points(phi[1], phi[2], pch = 21, bg = "green3")

P <- (logit(L * I) - mean(logit(L * I), na.rm = TRUE)) / sd(logit(L * I), na.rm = TRUE)
P <- 1 / (1 + exp(-P))
I <- P

# for a given point in the otolith, map each point of the contour to a proportion of that point:
clg()
image(1:nrow(I), 1:ncol(I), I, col = grey.colors(100), xlim = c(200, 1200), ylim = c(350, 1300))
points(contour[,1], contour[,2], cex = 0.25, col = "green3")
mx <- 810
my <- 810
for (p in seq(0, 1, len = 12)){
   xx <- mx + p * (contour[,1] - mx)
   yy <- my + p * (contour[,2] - my)
   points(xx, yy, cex = 0.2)
}

# Each pixel lies along a ray from the contour to the center point
#Calculate standardized distance to center point:
# Classify pixels into bands:

   
# Cone likelihood:
# divide otolith into bands, 
# each band has a mean and standard deviation





loglike <- function(phi, g = 9){
   # Parameters are the nucleus coordinates:
   mx <- phi[1]
   my <- phi[2]
   
   # Define vectors from nucleus to contour, and group them by relative distance:
   t <- seq(0, 1, len = 600)
   groups <- seq(min(t), max(t), len = g + 1)
   b <- NA * t
   for (i in 2:length(groups)){
      b[(t >= groups[i-1]) & (t <= groups[i])] <- i-1
   }
   ddx <- (contour[,1] - mx)
   ddy <- (contour[,2] - my)
   B <- L
   for (i in 1:length(t)){
      xx <- mx + t[i] * ddx
      yy <- my + t[i] * ddy
      
      B[unique(cbind(xx, yy))] <- b[i]
   }
   #image(1:nrow(I), 1:ncol(I), B, xlim = c(200, 1200), ylim = c(350, 1300))
   print(phi)
   points(phi[1], phi[2])
   groups <- 1:g
   ll <- 0
   for (i in 1:length(groups)){
      ix <- which(B == i)
      ll <- ll + mean(dnorm(I[ix], mean(I[ix], na.rm = TRUE) , sd(I[ix], na.rm = TRUE), log = TRUE), na.rm = TRUE)
   }
   
   return(-ll)
}

loglike(c(800, 800))
optim(c(800, 800), loglike, control = list(trace = 3))

clg()
image(1:nrow(I), 1:ncol(I), I, col = grey.colors(100), xlim = c(200, 1200), ylim = c(350, 1300))
points(contour[,1], contour[,2], cex = 0.25, col = "green3")
points(892.6437, 921.5713)


loglike <- function(theta, x, v, plot = FALSE){
   xp <- theta[["xp"]]
   yp <- theta[["yp"]]
   scale <- exp(theta[grep("log.scale", names(theta))])
   sigma <- exp(theta[["log.sigma"]])
   beta <- theta[["beta"]]
   
   ix <- x >= xp
   
   mu <- rep(NA, length(v))
   mu[ix] <- beta * ((x[ix]-xp) / scale[1])^2 + yp
   mu[!ix] <- beta * ((x[!ix]-xp) / scale[2])^2 + yp
   
   if (plot){
      plot(x, v)
      lines(x, mu, col = "red")
   }
   
   ll <- dnorm(v, mu, sigma, log = TRUE)
   
   return(-sum(ll))
}




theta["gamma"] <- 2

loglike(theta, x, v, plot = TRUE)

clg()
image(1:nrow(I), 1:ncol(I), I, col = grey.colors(100), xlim = c(200, 1200), ylim = c(350, 1300))
points(contour[,1], contour[,2], cex = 0.25, col = "green3")


theta <- c(xp = 800, 
           yp = 0.9, 
           log.sigma = -2, log.scale = c(4, 4), beta = -0.01, gamma = 2)
for (xx in seq(400, 1000, by = 25)){
   x <- 1:ncol(I)
   v <- I[xx, ]
   x <- x[!is.na(v)]
   v <- v[!is.na(v)]
   
   loglike(theta, x, v)
   theta <- optim(theta, loglike, x = x, v = v, control = list(trace = 3, maxit = 2000))$par
   theta <- optim(theta, loglike, x = x, v = v, control = list(trace = 3, maxit = 2000))$par
   theta <- optim(theta, loglike, x = x, v = v, control = list(trace = 3, maxit = 2000))$par
   points(xx, theta["xp"])
}

theta <- c(xp = 800, 
           yp = 0.9, 
           log.sigma = -2, log.scale = c(4, 4), beta = -0.01, gamma = 2)

for (xx in seq(400, 1000, by = 25)){
   x <- 1:nrow(I)
   v <- I[, xx]
   x <- x[!is.na(v)]
   v <- v[!is.na(v)]
   
   loglike(theta, x, v)
   theta <- optim(theta, loglike, x = x, v = v, control = list(trace = 3, maxit = 2000))$par
   theta <- optim(theta, loglike, x = x, v = v, control = list(trace = 3, maxit = 2000))$par
   theta <- optim(theta, loglike, x = x, v = v, control = list(trace = 3, maxit = 2000))$par
   #points(xx, theta["xp"])
   points(theta["xp"], xx, col = "blue")
}


loglike <- function(theta, x, v, plot = FALSE){
   xp <- theta[["xp"]]
   yp <- theta[["yp"]]
   scale <- exp(theta[grep("log.scale", names(theta))])
   sigma <- exp(theta[["log.sigma"]])
   beta <- theta[["beta"]]
   
   ix <- x >= xp
   
   mu <- rep(NA, length(v))
   mu[ix] <- beta * ((x[ix]-xp) / scale[1])^2 + yp
   mu[!ix] <- beta * ((x[!ix]-xp) / scale[2])^2 + yp
   
   if (plot){
      plot(x, v)
      lines(x, mu, col = "red")
   }
   
   ll <- dnorm(v, mu, sigma, log = TRUE)
   
   return(-sum(ll))
}

v <- I[900, ]
x <- 1:length(v)
x <- x[!is.na(v)]
v <- v[!is.na(v)]

ss <- NULL
#plot(x - xp, v, type = "l")
for (xp in 500:1100){
   vv <- v
   xx <- -(x - xp)
   vv <- vv[order(xx)]
   xx <- xx[order(xx)]
   
   xx[xx < 0] <- xx[xx < 0] / (-min(xx))
   xx[xx >= 0] <- xx[xx >= 0] / max(xx)
   xx[xx < 0] <- (xp - min(x)) * xx[xx < 0] 
   xx[xx >= 0] <- (max(x) - xp) * xx[xx >= 0]
   
   
   #plot(x - xp, v, type = "l")
   #lines(xx, vv, col = "red")
   #vline(0, col = "red", lty = "dashed")
   ss <- c(ss, sum((approx(xx, vv, x-xp)$y - v)^2))
}


loglike <- function(xp, x, v){
   vv <- v
   xx <- -(x - xp)
   vv <- vv[order(xx)]
   xx <- xx[order(xx)]
   
   xx[xx < 0] <- xx[xx < 0] / (-min(xx))
   xx[xx >= 0] <- xx[xx >= 0] / max(xx)
   xx[xx < 0] <- (xp - min(x)) * xx[xx < 0] 
   xx[xx >= 0] <- (max(x) - xp) * xx[xx >= 0]

   ss <- sum((approx(xx, vv, x-xp)$y - v)^2)
   
   return(ss)
}

fit <- function(x, v){
   x <- x[!is.na(v)]
   v <- v[!is.na(v)]
   if ((length(x) == 0) | all(is.na(v))) return(NA) 
   
   xp <- sum((v / sum(v)) * x)
   
   #xp <- optim(xp, loglike, x = x, v = v, control = list(maxit = 2000, trace = 3))$par
   
   return(xp)
}

clg()
image(1:nrow(I), 1:ncol(I), I, col = grey.colors(100), xlim = c(200, 1200), ylim = c(350, 1300))
points(contour[,1], contour[,2], cex = 0.25, col = "green3")
for(i in 450:1100){
   xp <- fit(1:ncol(I), I[i, ])
   points(i, xp, col = "red", cex = 0.5)
   
   xp <- fit(1:nrow(I), I[,i])
   points(xp, i, col = "green3", cex = 0.5)
}
   
II <- as.numeric(I)
XX <- as.numeric(repvec(1:ncol(I), nrow = nrow(I)))
YY <- as.numeric(repvec(1:nrow(I), ncol = ncol(I)))
ix <- which(!is.na(II))
xp <- sum((II[ix] / sum(II[ix])) * XX[ix])
yp <- sum((II[ix] / sum(II[ix])) * YY[ix])
points(xp, yp)
points(yp, xp)


points(800, 810)


v <- I[800, ]
x <- 1:length(v)
x <- x[!is.na(v)]
v <- v[!is.na(v)]

xp <- 810
vv <- v
xx <- -(x - xp)
vv <- vv[order(xx)]
xx <- xx[order(xx)]

xx[xx < 0] <- xx[xx < 0] / (-min(xx))
xx[xx >= 0] <- xx[xx >= 0] / max(xx)
xx[xx < 0] <- (xp - min(x)) * xx[xx < 0] 
xx[xx >= 0] <- (max(x) - xp) * xx[xx >= 0]

approx(xx, vv, x-xp)$y

pbeta()

ss <- sum((approx(xx, vv, x-xp)$y - v)^2)

xp <- 810
plot(x-xp, v)
points(xx, vv, col = "red")
vline(xp, col = "red", lty = "dashed")
pleft <- pbeta()

xleft <- - xx[xx < 0] / min(x-xp) + 1
plot(xleft, pbeta(xleft, 2, 2))

xleft <- - xx[xx < 0] / min(x-xp) + 1
plot(xleft, pbeta(xleft, 2, 2))


ix <- range(which(apply(I, 1, function(x) !all(is.na(x)))))


plot(x, v)
xp <- x[which.min(ss)]
vv <- v
xx <- -(x - xp)
vv <- vv[order(xx)]
xx <- xx[order(xx)]

xx[xx < 0] <- xx[xx < 0] / (-min(xx))
xx[xx >= 0] <- xx[xx >= 0] / max(xx)
xx[xx < 0] <- (xp - min(x)) * xx[xx < 0] 
xx[xx >= 0] <- (max(x) - xp) * xx[xx >= 0]
lines(xx + xp, vv, col = "red")

vline(xp, col = "red", lty = "dashed")

xp <- 820
v <- I[800, ]
x <- 1:ncol(I)
x <- x[!is.na(v)]
v <- v[!is.na(v)]

# Flip data around xp:
xx <- -(x - xp) + xp
vv <- v[order(xx)]
xx <- xx[order(xx)]

# Rescale coordinates on either side:
#xx[xx < xp]  <- (xp - min(x)) * (xx[xx < xp] - xp) / (xp - min(xx)) + xp
#xx[xx >= xp] <- (max(x) - xp) * (xx[xx >= xp]  - xp) / (max(xx) - xp) + xp

#xx[xx < xp]  <- (xx[xx < xp] - xp) / (xp - min(xx)) 
#xx[xx >= xp] <- (xx[xx >= xp]  - xp) / (max(xx) - xp) 

# Rescale coordinates using beta transform:
alpha <- 0.5
beta <- 0.5
xx[xx < xp] <- (xp - min(x)) * -pbeta(-(xx[xx < xp] - xp) / (xp - min(xx)), alpha, beta) + xp
alpha <- 1.5
beta <- 0.95
xx[xx >= xp] <- (max(x) - xp) * pbeta((xx[xx >= xp] - xp) / (max(xx) - xp), alpha, beta) + xp

plot(x, v, type = "l")
lines(xx, vv, col = "red")
vline(xp, col = "red")

points(x, approx(xx, vv, x)$y, col = "green2", cex = 0.5)

loglike <- function(theta, x, v, xp, plot = FALSE){
   # Parse model parameters:
   alpha <- exp(theta[grep("alpha", names(theta))])
   beta  <- exp(theta[grep("beta", names(theta))])
   if (missing(xp)){
      if (!("xp" %in% names(theta))) stop("'xp' must be specified.")
      xp <- theta[["xp"]]
   }
   
   # Flip data around xp:
   xx <- -(x - xp) + xp
   vv <- v[order(xx)]
   xx <- xx[order(xx)]
   
   # Rescale coordinates using beta transform:
   xx[xx < xp] <- (xp - min(x)) * -pbeta(-(xx[xx < xp] - xp) / (xp - min(xx)), alpha[1], beta[1]) + xp
   xx[xx >= xp] <- (max(x) - xp) * pbeta((xx[xx >= xp] - xp) / (max(xx) - xp), alpha[2], beta[2]) + xp
   
   ss <- sum((v - approx(xx, vv, x)$y)^2)
   
   if (plot){
      plot(x, v, type = "l")
      lines(xx, vv, col = "red")
      vline(xp, col = "red")
      points(x, approx(xx, vv, x)$y, col = "green2", cex = 0.5)
   }
   
   return(ss)
}

loglike(theta, x, v, 810)
xp <- 750
parameters <- NULL
for (xp in 700:900){
   print(xp)
   theta <- c(alpha = c(0, 0), beta = c(0, 0))
   theta <- optim(theta, loglike, x = x, v = v, xp = xp, control = list(trace = 0))$par
   ss <- loglike(theta, x, v, xp, plot = FALSE)
   parameters <- rbind(parameters, c(ss, xp, theta))
}


loglike(parameters[which.min(parameters[,1]), 3:6], x, v, xp = parameters[which.min(parameters[,1]), 2], plot = TRUE)

clg()
image(1:nrow(I), 1:ncol(I), I, col = grey.colors(100), xlim = c(200, 1200), ylim = c(350, 1300))
points(contour[,1], contour[,2], cex = 0.25, col = "green3")

for (i in 600:1000){
   v <- I[i, ]
   x <- 1:ncol(I)
   x <- x[!is.na(v)]
   v <- v[!is.na(v)]
   
   xp = sum((v / sum(v)) * x)
   theta <- c(alpha = c(0, 0), beta = c(0, 0))
   theta <- optim(theta, loglike, x = x, v = v, xp = xp, control = list(trace = 3, maxit = 2000))$par
   loglike(theta, x, v, xp = xp)
   
   theta <- c(xp = xp, theta)
   #loglike(theta, x, v)
   theta <- optim(theta, loglike, x = x, v = v, control = list(trace = 3, maxit = 2000))$par
   #loglike(theta, x, v, plot = TRUE)
   #points(theta[["xp"]], i, col = "green", cex = 0.25)
   points(i, theta[["xp"]], col = "red", cex = 0.25)
}


image(I[600:1000, 600:1000])





