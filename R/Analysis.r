rm(list = ls())
setwd("C:/Users/SuretteTJ/Desktop/Ageing")
source("source.R")

n <- 10
p <- 3

knots <- c(seq(0, 1, len = 10), 1, 1, 1)
b <- bspline(knots = knots)

v <- -10*evaluate(derive(b$basis)[[9]], 1)/evaluate(derive(b$basis)[[8]], 1) 
coef(b) <- c(rnorm(7), 10, 10)
plot(b)


step <- 1/n
knots <- seq(-p*step, 1+p*step, by = step)
b <- bspline(knots, degree = p)
plot(derive(b))


knots <- seq(0, 1, len = n+1)
#knots <- c(knots, knots[1:(p+2)])

b <- bspline(8, type = "per")
coef(b) <- runif(b$dim) 
t <- seq(0, 1, len = 1000)
x <- evaluate(b, t)
coef(b) <- runif(b$dim) 
y <- evaluate(b, t)
plot(x, y)

       
w <- rnorm(b$dim-1)
w <- c(w, w[1])
coef(b) <- w
plot(b)

windows()
layout(matrix(1:3, ncol = 1))
x <- seq(0, 2, len = 1000)
plot(x, evaluate(b, x %% 1), type = "l", lwd = 2)
plot(x, evaluate(derive(b), x %% 1), type = "l", lwd = 2)
plot(x, evaluate(derive(derive(b)), x %% 6), type = "l", lwd = 2)

plot(b$basis)

x <- seq(0, 1, len = 100)
y <- bs(x, intercept = TRUE)
plot(c(0, 1), c(0, 1), type = "n")
for (i in 1:dim(y)[2]){
   lines(x, y[,i])
}

y <- bs(x, intercept = FALSE)
plot(c(0, 1), c(0, 1), type = "n")
for (i in 1:dim(y)[2]){
   lines(x, y[,i], col = "red", lwd = 2)
}

b <- bspline(c(0, 0, 0, 0, 1, 1, 1, 1))




