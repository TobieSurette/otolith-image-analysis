library(jpeg)
library(gulf.graphics)

# Read image:
I <- readJPEG("photos/examples/Plaice.jpg", native = FALSE)

# Convert to greyscale:
I <- 0.2989 * I[,,1] + 0.5870 * I[,,2] +  0.1141 * I[,,3]

# Display image:
#image(1:nrow(I), 1:ncol(I), I, col = grey.colors(100))
ix <- seq(1, nrow(I), by = 3)
iy <- seq(1, ncol(I), by = 3)
image((1:nrow(I))[ix], (1:ncol(I))[iy], (I[ix, iy] - mean(I[ix, iy], na.rm = TRUE)), col = grey.colors(100))

Gx <- I[3:nrow(I), 2:(ncol(I)-1)] - I[1:(nrow(I)-2), 2:(ncol(I)-1)]
Gy <- I[2:(nrow(I)-1), 3:ncol(I)] - I[2:(nrow(I)-1), 1:(ncol(I)-2)]
Mag <- sqrt(Gx*Gx + Gy*Gy) 
theta <- atan2(Gy,Gx)

image(theta, col = grey.colors(100))

# Identify contour points:
dx <- 2:(nrow(I)-1)
dy <- 2:(ncol(I)-1)
ix <- !is.na(I[dx,dy]) & (is.na(I[dx-1,dy]) | is.na(I[dx+1,dy]) | is.na(I[dx,dy-1]) | is.na(I[dx,dy+1]))
image(ix)
p <- which(ix, arr.ind = TRUE)
plot(p[, 1], p[, 2], cex = 0.5)
d <- dist(as.matrix(p))
H <- hclust(d)
ix <- cutree(H, 10)

plot(p[,1], p[,2], cex = 0.5)
points(p[ix==1,1], p[ix==1,2], col = "red")
points(p[ix==2,1], p[ix==2,2], col = "green2")
points(p[ix==3,1], p[ix==3,2], col = "blue")
points(p[ix==4,1], p[ix==4,2], col = "yellow4")
points(p[ix==5,1], p[ix==5,2], col = "purple")

d <- as.matrix(d)
hist(I, n = 100)

I[I < 0.3] <- NA
vline(1000)
plot(I[1000,])
plot(I[,1000], type = "l")


x0 <- 805
x1 <- 240
y0 <- 805
y1 <- 1035
t <- seq(0, 1, len = 1000)
xx <- x0 + (x1-x0) * t
yy <- y0 + (y1-y0) * t
xx <- unique(cbind(round(xx), round(yy)))
yy <- xx[,2]
xx <- xx[,1]
dd <- sqrt((xx - x0)^2 + (yy - y0)^2)

lines(xx, yy)
points(805, 805, pch = 21, bg = "red2")
points(240, 1035, pch = 21, bg = "red2")

library(gulf.graphics)
rx <- c(262, 283, 302, 342, 392, 445, 510, 660)
ry <- rep(NA, length(rx))
ri <- rep(NA, length(rx))
for (i in 1:length(rx)){
   ri[i] <- which.min(abs(xx-rx[i]))
   ry[i] <- yy[ri[i]]
}
points(rx, ry, pch = 21, bg = "red2")

zz <- rep(NA, length(xx))
for (i in 1:length(xx)) zz[i] <- I[xx[i], yy[i]]
plot(dd, zz, type = "l")
points(dd[ri], zz[ri], pch = 21, bg = "red")

plot(1:length(ri), dd[ri])
tt <- 1:length(ri)
tt2 <- tt * tt
model <- lm(dd[ri] ~ tt2 + tt + 1)
points(tt, predict(model), pch = 21, bg = "red")

dm <- predict(model, newdata = list(tt = seq(0, 9, len = 1000), tt2 = seq(0, 9, len = 1000)^2))
age <- approx(dm, seq(0, 9, len = 1000), dd)$y
plot(dd, age)

plot(age, zz, type = "l")
vline(1:8, lty = "dashed")

plot(dd, zz)
model <- lm(zz ~ poly(dd, 2))
lines(dd, predict(model), col = "red")

plot(log(max(dd) - dd), zz - predict(model), type = "l", lwd = 2, col = "red3", xlim = c(3, 7), ylim = c(-0.1, 0.05))
vline(log(max(dd) - dd[ri]), lty = "dashed")


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
x0 <- 805
y0 <- 805
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


theta = c(x0 = 805, y0 = 805,
          radius = 400, 
          mu = c(mean(I[cluster == 1]), mean(I[cluster == 2])), 
          log.sigma = log(c(sd(I[cluster == 1]), sd(I[cluster == 2]))))

x0 <- 805
y0 <- 805
xx <- repvec(1:nrow(I), ncol = ncol(I)) 
yy <- repvec(1:ncol(I), nrow = nrow(I)) 

loglike <- function(theta, xx, yy, I){
   dd <- sqrt((xx - x0)^2 + (yy - y0)^2)
   ix <- which(dd <= theta["radius"])
   mu <- theta[grep("mu", names(theta))]
   sigma <- exp(theta[grep("log.sigma", names(theta))])
   v <- sum(dnorm(I[ix], mu[1], sigma[1], log = TRUE)) + sum(dnorm(I[-ix], mu[2], sigma[2], log = TRUE))
   return(-v)
}
   
loglike(theta, xx, yy, I)
theta <- optim(theta, loglike, xx = xx, yy = yy, I = I, control = list(trace = 3))$par   
   
clg()
ix <- seq(1, nrow(I), by = 2)
iy <- seq(1, ncol(I), by = 2)
image((1:nrow(I))[ix], (1:ncol(I))[iy], I[ix,iy], col = grey.colors(100))

angle <- seq(0, 2*pi, len = 100)
points(theta["x0"] + theta["radius"] * cos(angle), theta["y0"] + theta["radius"] * sin(angle), pch = 21, bg = "red")



