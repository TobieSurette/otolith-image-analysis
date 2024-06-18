setwd("C:/Users/SuretteTJ/Desktop/Ageing")

library(jpeg)
library(gulf)

source("source.R")
graphics.off()

I <- readJPEG("Images/Plaice.jpg", native = FALSE)
I <- as.image(I)
I <- smooth(I, 3)
#I$z <- log(I$z/(1-I$z))
#I <- thin(I, n = 2)
 
O <- as.otolith(I, threshold = 0.25)
nucleus(O) <- c(800, 400)
O$image$z <- log(O$image$z/(1-O$image$z))
#O$image <- quantile(O$image)
G <- gradient(O$image, polar = FALSE)

I <- trim(O$image)
I <- trim(I)
I <- quantile(I)
I$x <- 1:dim(I$z)[1]
I$y <- 1:dim(I$z)[2]
G <- gradient(I)$G
theta <- gradient(I)$theta
xg <- G$z * cos(theta$z)
yg <- G$z * sin(theta$z)

dG <- as.data.frame(G)
dt <- as.data.frame(theta)

# Nucleus sum-of-squares:
SS <- function(p, x, y, r, theta){
   xn <- p[1]
   yn <- p[2]
   dd <- sqrt((x - xn)^2 + (y - yn)^2)
   xu <- (x-xn)/dd
   yu <- (y-yn)/dd
   v <- (r * cos(theta) * xu) + (r * sin(theta) * yu)
   SS <- sum(v, na.rm = TRUE)
   points(p[1], p[2], col = "red")
   return(SS)
}

p <- optim(c(1000, 800), SS, x = dG$x, y = dG$y, r = dG$z, theta = dt$z)$par
points(p[1], p[2], col = "red")
   
while (TRUE){
   plot(I)
   p <- draw.points()
   SS(c(p$x, p$y), dG$x, dG$y, dG$z, dt$z)
   
   xn <- p$x
   yn <- p$y
   xx <- repvec(I$x, ncol = length(I$y)) - xn
   yy <- repvec(I$y, nrow = length(I$x)) - yn
   dd <- sqrt((xx)^2 + (yy)^2)
   xu <- xx / dd
   yu <- yy / dd

   v <- (xg * xu) + (yg * yu)
   SS <- sum(v, na.rm = TRUE)
}

plot(I)
dG <- as.data.frame(G)
dt <- as.data.frame(theta)
w <- dG$z / sum(dG$z)
for (i in 1:1000){
   index <- sample(1:dim(dG)[1], 1, prob = w)
   m <- atan2(dG[index,]$z * sin(dt[index,]$z), dG[index,]$z * cos(dt[index,]$z))
   b <- dG[index,]$y - m * dG[index,]$x
   abline(b, m, lwd = 0.1)
}



