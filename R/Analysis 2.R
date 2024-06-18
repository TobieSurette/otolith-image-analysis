I <- readJPEG("T678_1237.jpg", native = FALSE)
I <- as.image(I)

I$z[I$z < 0.25] <- NA
I <- quantile.image(I)

Z <- sweep(I, 810, 400, angles = seq(0, 2*pi, len = 500), n = 500)

temp <- repvec(apply(Z$z, 2, mean, na.rm = TRUE), nrow = dim(Z$z)[1])
index <- is.na(Z$z)
Z$z[index] <- temp[index]

t <- seq(0, 1, len = dim(Z$z)[2])
#for (i in 1:dim(Z$z)[1]){
#   Z$z[i,] <- approx(t^3, Z$z[i,], t)$y
#}

# Design matrix:
bw <- ispline(c(0, 0, seq(0, 1, len = 6), 1, 1))
X <- cbind(t, evaluate(bw, t))
      
k <- 300

# Detrend signal:
R <- matrix(NA, nrow = k, ncol = length(zt))
for (i in 1:k){
   cat(paste(i, "\n"))
   R[i, ] <- mgam(t, 1-Z[i, ]$z[1,])[, 1]
}

M  <- matrix(NA, nrow = k, ncol = 7)
XW <- matrix(NA, nrow = k, ncol = dim(R)[2])
ZM <- matrix(NA, nrow = k, ncol = dim(R)[2])
f <- approxfun(t, R[1, ], rule = 2)

theta <- rep(-3, 7)
warp.SS(theta,  t = t, z = R[1, ], zt = R[1, ], fint = f, X = X, plot = TRUE)
theta <- optim(theta, warp.SS, t = t, z = R[i, ], zt = R[1, ], fint = f, X = X, control = list(trace = 0))$par
for (i in  2:k){
   res <- optim(theta, warp.SS, t = t, z = R[i, ], zt = R[1, ], fint = f, X = X, control = list(trace = 0))
   M[i, ] <- res$par
   theta <- res$par
   warp.SS(theta, t = t, z = R[i, ], zt = R[1, ], fint = f, X = X, plot = TRUE)
   title(main = i)
   
   p <- c(0, theta)
   p <- exp(p) / sum(exp(p))
   XW[i, ] <- (X %*% p)[, 1]
   ZM[i, ] <- approx(XW[i,], R[i, ], t)$y
  # index <- round(XW[i, ]*(length(zt)-1))+1
  # ZM[i, ] <- R[i,index]
}

#
# Detrend Z
windows()
image(1:dim(R)[1], 1:dim(ZM)[2], R, useRaster = TRUE,
      breaks = seq(min(R), max(R), len = 1001),
      col = gray(seq(0, 1, len = 1000)))
windows()
image(1:dim(ZM)[1], 1:dim(ZM)[2], ZM, useRaster = TRUE,
      breaks = seq(min(ZM, na.rm = TRUE), max(ZM, na.rm = TRUE), len = 1001),
      col = gray(seq(0, 0.9, len = 1000)))

