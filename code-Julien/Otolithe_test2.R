library(MASS)
library(readxl)
library(imager)

#rm(list=ls()) # effacer la mémoire de l'environnement de travail
#ls() # vérifier la liste des objets
set.seed(2024) # Résultats reproductibles

image_originale = load.image("C:/Users/THIBODEAUJU/Documents/Otolithe/Plaice2.jpg")

image_grise = grayscale(image_originale)

hist(image_grise)



tc = optim(c(0.2,30),TrouverContour,image = image_grise)$par
image = as.cimg(clean(image_grise>0.2,tc[2]))

image_grise = as.matrix(image_grise[,,1,1])

noyau = Trouver_noyau2(0.2, tc[2], as.cimg(image_grise))

centre = noyau

x1=x2=c()
y1=y2=c()
angles = seq(0,2*pi, length.out = 200)
anglesD = seq(0,2*pi, length.out = 1000)
I=as.matrix(image)
for(theta in angles){
  r1=r2=0
  k1=k2=100
  while(TRUE){
    J1 = I[cbind(floor(centre[1,1]+(r1+k1)*cos(theta)), floor(centre[1,2]+(r1+k1)*sin(theta)))]
    J2 = I[cbind(floor(centre[2,1]+(r2+k2)*cos(theta)), floor(centre[2,2]+(r2+k2)*sin(theta)))]
    if(J1 == 1){
      r1=r1+k1
    }else{
      k1=k1*0.1
    }
    if(J2 == 1){
      r2=r2+k2
    }else{
      k2=k2*0.1
    }
    if(k1<1e-8 && k2<1e-8) break
  }
  x1=c(x1,floor(centre[1,1]+r1*cos(theta)))
  y1=c(y1,floor(centre[1,2]+r1*sin(theta)))
  x2=c(x2,floor(centre[2,1]+r2*cos(theta)))
  y2=c(y2,floor(centre[2,2]+r2*sin(theta)))
}
plot(x1,y1,type = "l")
lines(x2,y2)

Contour1=matrix(data=c(x1,y1), nrow = length(x1), ncol = 2)
Contour2=matrix(data=c(x2,y2), nrow = length(x2), ncol = 2)

#noyau = Trouver_noyau(tc[1], tc[2], as.cimg(image_grise))
axe_principal1 = Trouver_axe_principal(noyau[1,],Contour1)
axe_principal2 = Trouver_axe_principal(noyau[2,],Contour2)

Tranche1 = Tranche_Contour(image, noyau[1,], Contour1, n_points = 1000)
Tranche2 = Tranche_Contour(image, noyau[2,], Contour2, n_points = 1000)


#derive_extremiteX = ((x[2]-x[1])/(angles[2]-angles[1])+(x[length(x)]-x[length(x)-1])/(angles[length(angles)]-angles[length(angles)-1]))
SX1 = spline(angles, x1, method = "natural")

#derive_extremiteY = ((y[2]-y[1])/(angles[2]-angles[1])+(y[length(y)]-y[length(y)-1])/(angles[length(angles)]-angles[length(angles)-1]))
SY1 = spline(angles, y1, method = "natural")

#derive_extremiteX = ((x[2]-x[1])/(angles[2]-angles[1])+(x[length(x)]-x[length(x)-1])/(angles[length(angles)]-angles[length(angles)-1]))
SX2 = spline(angles, x2, method = "natural")

#derive_extremiteY = ((y[2]-y[1])/(angles[2]-angles[1])+(y[length(y)]-y[length(y)-1])/(angles[length(angles)]-angles[length(angles)-1]))
SY2 = spline(angles, y2, method = "natural")

plot(angles, x1, type = "o", xlab = "angle", ylab = "", col = "red", lwd = 2,
     ylim = c(0,3000))
lines(SX1$x,SX1$y,lwd=2)
lines(angles, y1, type = "o", col = "blue", lwd = 2)
lines(SX1$x,SY1$y,lwd=2)
grid()
legend(x="bottomleft", legend = c("x1","y1"), fill = c("red","blue"),cex=0.7)

plot(angles, x2, type = "o", xlab = "angle", ylab = "", col = "red", lwd = 2,
     ylim = c(0,3000))
lines(SX2$x,SX2$y,lwd=2)
lines(angles, y2, type = "o", col = "blue", lwd = 2)
lines(SY2$x,SY2$y,lwd=2)
grid()
legend(x="bottomleft", legend = c("x2","y2"), fill = c("red","blue"),cex=0.7)

plot(x1,y1,type="o",col="green",lwd=2,xlab="x1",ylab="y1")
lines(SX1$y,SY1$y,lwd=2)
plot(x2,y2,type="o",col="green",lwd=2,xlab="x2",ylab="y2")
lines(SX2$y,SY2$y,lwd=2)

plot(image_originale)
lines(SX1$y,SY1$y,col="red",lwd=3)
lines(SX2$y,SY2$y,col="red",lwd=3)
lines(axe_principal1[c(1,3)],axe_principal1[c(2,4)],col="purple",lwd=2)
lines(axe_principal2[c(1,3)],axe_principal2[c(2,4)],col="purple",lwd=2)

#noyau = Trouver_noyau(Grain_filter(image_grise,3),Contour,100)

points(noyau[1,1],noyau[1,2],col="red")
points(noyau[2,1],noyau[2,2],col="red")

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

SweepD = t(Sweep_anneaux(noyau[1,],image_grise,Tranche1,1000))
SweepG = t(Sweep_anneaux(noyau[2,],image_grise,Tranche2,1000))

fG = ecdf(SweepG)
fD = ecdf(SweepD)

Sweep2G = matrix(fG(SweepG),nrow(SweepG),ncol(SweepG))
Sweep2D = matrix(fG(SweepD),nrow(SweepD),ncol(SweepD))


plot(as.cimg(Sweep2G),xlim = c(0,length(Sweep2G[,1])), ylim=c(0,length(Sweep2G[1,])),axes = F)
axis(1, at = seq(1,length(Sweep2G[,1]),len=5), labels = c(0, "pi/2", "pi", "3pi/2", "2pi"))
axis(2, at = seq(0,length(Sweep2G),5))
title(main = "Sweep de l'otolithe", xlab = "Angle", ylab = "Points pris des rayons")

plot(as.cimg(Sweep2D),xlim = c(0,length(Sweep2D[,1])), ylim=c(0,length(Sweep2D[1,])),axes = F)
axis(1, at = seq(1,length(Sweep2G[,1]),len=5), labels = c(0, "pi/2", "pi", "3pi/2", "2pi"))
axis(2, at = seq(0,length(Sweep2G),5))
title(main = "Sweep de l'otolithe", xlab = "Angle", ylab = "Points pris des rayons")


TransformeG = Transformation(noyau[2,],Tranche2,Sweep2G)
TransformeD = Transformation(noyau[1,],Tranche1,Sweep2D)

paramG=TransformeG$param
paramD=TransformeD$param

Sweep_WarpG = TransformeG$Transforme
Sweep_WarpD = TransformeD$Transforme

plot(as.cimg(Sweep_WarpG),xlim = c(0,length(Sweep_WarpG[,1])), ylim=c(0,length(Sweep_WarpG[1,])), axes=F)
axis(1, at = seq(1,length(Sweep_WarpG[,1]),len=5), labels = c(0, "pi/2", "pi", "3pi/2", "2pi"))
axis(2, at = seq(0,length(Sweep2G),5))
title(main = "Warp du sweep de l'otolithe (init = param du dernier)", xlab = "Angle", ylab = "Points pris des rayons")

plot(as.cimg(Sweep_WarpD),xlim = c(0,length(Sweep_WarpD[,1])), ylim=c(0,length(Sweep_WarpD[1,])), axes=F)
axis(1, at = seq(1,length(Sweep_WarpD[,1]),len=5), labels = c(0, "pi/2", "pi", "3pi/2", "2pi"))
axis(2, at = seq(0,length(Sweep2G),5))
title(main = "Warp du sweep de l'otolithe (init = param du dernier)", xlab = "Angle", ylab = "Points pris des rayons")


# plot(Sweep2[axe_principal[5],],type = "l",lwd=3)
# couleur = c("blue","red","green","orange")
# c=1
# for(i in seq(1,length(Sweep_Warp[,1]),len=4)){
#   lines(Sweep_Warp[i,],col=couleur[c],lwd=1)
#   c=c+1
# }
# legend(x="bottomleft",legend = seq(1,length(Sweep_Warp[,1]),len=4), fill = couleur, cex=0.7)

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



