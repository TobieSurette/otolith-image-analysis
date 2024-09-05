setwd("C:/Users/THIBODEAUJU/Documents/Otolithe")

library(MASS)
library(readxl)
library(imager)

set.seed(2024) # Résultats reproductibles

image_originale = load.image("C:/Users/THIBODEAUJU/Documents/Otolithe/Plaice2.jpg")

image_grise = grayscale(image_originale)

hist(image_grise)

TrouverContour = function(x,image){
  IM = as.cimg(clean(image>x[1],x[2]))
  return((mean(image[IM==1])-1)^2+100*(mean(image[IM==0])-0)^2)
}



tc = optim(c(0.2,30),TrouverContour,image = image_grise)$par
image = as.cimg(clean(image_grise>tc[1],tc[2]))

image_grise = as.matrix(image_grise[,,1,1])

# ContourX = c()
# ContourY = c()
#
# for(i in seq(1,length(image[,1]),10)){
#   for(j in seq(2,length(image[1,]),10)){
#     a=image[i,j-1]-image[i,j]
#     if(a == 1 || a == -1){
#       ContourX = c(ContourX, i)
#       ContourY = c(ContourY, j)
#     }
#   }
# }

centre = c(x=nrow(image_grise)/2, y=ncol(image_grise)/2)

x=c()
y=c()
angles = seq(0,2*pi, length.out = 200)
anglesD = seq(0,2*pi, length.out = 1000)
for(theta in angles){
  r=0
  k=100
  while(TRUE){
    if(image[floor(centre["x"]+(r+k)*cos(theta)),floor(centre["y"]+(r+k)*sin(theta))] == 1){
      r=r+k
    }else{
      k=k*0.1
    }
    if(k<1e-10) break
  }
  x=c(x,floor(centre["x"]+r*cos(theta)))
  y=c(y,floor(centre["y"]+r*sin(theta)))
}
plot(x,y,type = "l")

Contour=matrix(data=c(x,y), nrow = length(x), ncol = 2)

noyau = Trouver_noyau(tc[1], tc[2], as.cimg(image_grise))
axe_principal = Trouver_axe_principal(noyau,Contour)

Tranche = Tranche_Contour(image, noyau, Contour, n_points = 1000)

SX = spline(angles, x)

SY = spline(angles, y)

plot(angles, x, type = "o", xlab = "angle", ylab = "", col = "red", lwd = 2,
     ylim = c(200,1300))
lines(SX$x,SX$y,lwd=2)
lines(angles, y, type = "o", col = "blue", lwd = 2)
lines(SY$x,SY$y,lwd=2)
grid()
legend(x="bottomleft", legend = c("x","y"), fill = c("red","blue"),cex=0.7)

plot(x,y,type="o",col="green",lwd=2,xlab="x",ylab="y")
lines(SX$y,SY$y,lwd=2)

plot(image_originale)
lines(SX$y,SY$y,col="red",lwd=3)
lines(axe_principal[c(1,3)],axe_principal[c(2,4)],col="purple",lwd=2)

#noyau = Trouver_noyau(Grain_filter(image_grise,3),Contour,100)

points(noyau[1],noyau[2],col="red")

Contraste = function(x,c=1){
  A = log(x/(1-x))
  A = c*(A-mean(A))/sd(A)
  return(1/(1+exp(-A)))
}

Grain_filter = function(image, a){
  new_image = image
  for(i in 1:length(image[,1])){
    for(j in 1:length(image[1,])){
      a1 = i-a
      a2 = i+a
      a3 = j-a
      a4 = j+a

      if(i-a < 0) a1 = 1
      if(i+a > length(image[,1])) a2 = length(image[,1])
      if(j-a < 0) a3 = 1
      if(j+a > length(image[1,])) a4 = length(image[1,])

      new_image[i,j] = mean(image[a1:a2, a3:a4])
    }
  }
  return(new_image)
}

Sweep = t(Sweep_anneaux(noyau,image_grise,Tranche,1000))

#f = ecdf(Sweep)

#Sweep2 = matrix(f(Sweep),nrow(Sweep),ncol(Sweep))

Sweep2 = Contraste(Sweep,1)


plot(as.cimg(Sweep2),xlim = c(0,length(Sweep2[,1])), ylim=c(0,length(Sweep2[1,])),axes = F)
axis(1, at = seq(1,length(Sweep2[,1]),len=5), labels = c(0, "pi/2", "pi", "3pi/2", "2pi"))
axis(2, at = seq(0,length(Sweep2),5))
title(main = "Sweep de l'otolithe", xlab = "Angle", ylab = "Points pris des rayons")


# On trouve l'image qui a été warpé
Transforme = Transformation(noyau,Tranche,Sweep2)

param=Transforme$param

Sweep_Warp = Transforme$Transforme

plot(as.cimg(Sweep_Warp),xlim = c(0,length(Sweep_Warp[,1])), ylim=c(0,length(Sweep_Warp[1,])), axes=F)
axis(1, at = seq(1,length(Sweep_Warp[,1]),len=5), labels = c(0, "pi/2", "pi", "3pi/2", "2pi"))
axis(2, at = seq(0,length(Sweep2),5))
title(main = "Warp du sweep de l'otolithe (init = param du dernier)", xlab = "Angle", ylab = "Points pris des rayons")


plot(colMeans(Sweep_Warp),type="l",lwd=3, xaxt = "n", ylim=c(0,1),
     xlab = "Noyau à contour", ylab = "Intensité de gris (0 = noir, 1 = blanc)")
axis(1, at = seq(1,length(colMeans(Sweep_Warp)),len=11), labels = seq(0,1,len=11))
grid()


dx <- imgradient(as.cimg(Sweep_Warp),"x")
dy <- imgradient(as.cimg(Sweep_Warp),"y")
grad.mag <- sqrt(dx^2+dy^2)
plot(grad.mag,main="Gradient magnitude")

plot(c(0,100), c(0,100),type="n")


for(i in 1:nrow(param)){
  T = Fonction_Warp(param[i,1]^2, param[i,2]^2, 100)
  lines(1:100,T)
}



