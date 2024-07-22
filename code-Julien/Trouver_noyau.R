
library(zoo)

euclidien = function(a, b) sqrt(sum((a - b)^2))

# Fait un sweep de tout le tour de l'otolithe à partir du noyau et 
# l'étend 
Sweep_anneaux = function(point, image, Contour, n_points){
  L = matrix(,nrow = n_points, ncol = length(Contour[,1]))
  
  for(j in 1:length(Contour[,1])){
    indiceX = floor(seq(point[1],Contour[j,1],length.out=n_points))
    indiceY = floor(seq(point[2],Contour[j,2],length.out=n_points))
    L[,j] = image[cbind(indiceX, indiceY)]
  }
  
  return(L)
}

Tranche_Contour = function(image,noyau,Contour,angle = pi/4, n_points = 100){
  axe_p = Trouver_axe_principal(noyau, Contour)[5]
  angle_p = seq(0,2*pi,len=length(Contour[,1]))[axe_p]
  
  angles = seq(0,2*pi, length.out = length(Contour[,1]))
  
  angles_tranche = seq(angle_p - angle/2, angle_p + angle/2, length.out = n_points)
  for(i in 1:n_points){
    if(angles_tranche[i] < 0){
      angles_tranche[i] = angles_tranche[i]+2*pi
    }
  }
  
  TrancheX = floor(spline(angles, Contour[,1], xout = angles_tranche)$y)
  TrancheY = floor(spline(angles, Contour[,2], xout = angles_tranche)$y)
  
  return(matrix(c(TrancheX,TrancheY), length(TrancheX), 2))
}

Enrouler = function(Sweep, Contour, point){
  n_points=length(Sweep[1,])
  
  L = matrix(,nrow = max(Contour[,1]), ncol = max(Contour[,2]))
  
  for(j in 1:length(Contour[,1])){
    indiceX = floor(seq(point[1],Contour[j,1],length.out=n_points))
    indiceY = floor(seq(point[2],Contour[j,2],length.out=n_points))
    L[cbind(indiceX,indiceY)] = Sweep[j,]
  }
  return(L)
}

Similitude_des_rayons = function(point, image, Contour, n_points){
  points(point[1],point[2])
  L = Sweep_anneaux(point,image,Contour,n_points)
  
  #L = TransformationVieux(point,Contour,Sweep)
  
  moy=c()
  # Prend la moyenne de gris pour chaques anneaux.
  for(i in 1:length(L[,1])){
    moy = c(moy, mean(L[i,]))
  }
  
  
  S=0
  
  for(i in 1:length(L[1,])){
    S=(S+sum((L[,i]-moy)^2))
  }
  return(S)
}

Trouver_noyau = function(i, g, image){
  while(sum(clean(image>i+0.01, g))>1){
    i=i+0.01
  }
  noyau = which(matrix(clean(image>i, g),nrow(image),ncol(image)) == 1, arr.ind = TRUE)
  
  centre = c(mean(noyau[,1]), mean(noyau[,2]))
  
  
  return(round(centre))
}




TrouverFormes = function(x,image){
  IM = as.cimg(clean(image>x[1],x[2]))
  return((mean(image[IM==1])-1)^2+100*(mean(image[IM==0])-0)^2)
}

# TrouverContour = function(image_grise){
#   #tc = optim(c(0.2,30),TrouverFormes,image = as.cimg(image_grise))$par
#   image = as.vector(clean(as.cimg(image_grise)>0.2,30))
#   image = matrix(image, nrow(image_grise), ncol(image_grise))
# 
#   centre = c(x=nrow(image_grise)/2, y=ncol(image_grise)/2)
# 
#   x=c()
#   y=c()
#   angles = seq(0,2*pi, length.out = 200)
#   for(theta in angles){
#     r=0
#     k=100
#     while(TRUE){
#       indice = c(floor(centre["x"]+(r+k)*cos(theta)),floor(centre["y"]+(r+k)*sin(theta)))
#       if(image[indice[1], indice[2]] == 1){
#         r=r+k
#       }else{
#         k=k*0.1
#       }
#       if(k<1e-10) break
#     }
#     x=c(x,floor(centre["x"]+r*cos(theta)))
#     y=c(y,floor(centre["y"]+r*sin(theta)))
#   }
#   Contour=matrix(data=c(x,y), nrow = length(x), ncol = 2)
# 
#   return(list(coord = Contour, param = c(0.2,30)))
# }

TrouverContour = function(image_grise,noyau){
  #tc = optim(c(0.2,30),TrouverFormes,image = as.cimg(image_grise))$par
  image = as.vector(clean(as.cimg(image_grise)>0.2,30))
  image = matrix(image, nrow(image_grise), ncol(image_grise))
  
  centre = noyau
  
  nb_centre = 3
  NP = which(image == 1, arr.ind = TRUE)
  rand_i = runif(nb_centre,1,nrow(NP))
  
  x=c()
  y=c()
  angles = seq(0,2*pi, length.out = 400)
  
  for(j in 1:nb_centre){
    for(theta in angles){
      r=0
      k=100
      while(TRUE){
        indice = c(floor(centre[1]+(r+k)*cos(theta)),floor(centre[2]+(r+k)*sin(theta)))
        if(indice[1] > nrow(image) || indice[1] < 1 || indice[2] > ncol(image) || indice[2] < 1){
          k=k*0.1
        }
        else {
          if(image[indice[1], indice[2]] == 1){
            r=r+k
          }else{
            k=k*0.1
          }
        }
        if(k<1e-10) break
      }
      x=c(x,floor(centre[1]+r*cos(theta)))
      y=c(y,floor(centre[2]+r*sin(theta)))
    }
    
    centre = c(NP[rand_i[j], 1], NP[rand_i[j], 2])
  }
  
  Contour=matrix(data=c(x,y), nrow = length(x), ncol = 2)
  
  return(list(coord = Contour, param = c(0.2,30)))
}

# TrouverContour = function(image_grise){
#   #tc = optim(c(0.2,30),TrouverFormes,image = as.cimg(image_grise))$par
#   image = as.vector(clean(as.cimg(image_grise)>0.2,30))
#   image = matrix(image, nrow(image_grise), ncol(image_grise))
# 
#   centre = c(x=nrow(image_grise)/2, y=ncol(image_grise)/2)
# 
#   x=c()
#   y=c()
#   jump = c()
#   angles = seq(0,2*pi, length.out = 200)
#   for(i in 1:length(angles)){
#     theta = angles[i]
#     r=0
#     k=100
#     while(TRUE){
#       while(TRUE){
#         indice = c(floor(centre["x"]+(r+k)*cos(theta)),floor(centre["y"]+(r+k)*sin(theta)))
#         if(image[indice[1], indice[2]] == 1){
#           r=r+k
#         }else{
#           k=k*0.1
#         }
#         if(k<1e-10) break
#       }
#       if (i>1){
#         if(abs(jump[length(jump)] < 30)){
#           jump = c(jump, euclidien(c(x[length(x)], y[length(y)]), c(floor(centre["x"]+r*cos(theta)), floor(centre["y"]+r*sin(theta)))))
#           break
#         } else{
#           theta = theta-0.001*(angles[i] - angles[i-1])
#           if(theta < angles[i-1]) break
#         }
#       }else{
#         jump = c(jump, euclidien(c(x[length(x)], y[length(y)]), c(floor(centre["x"]+r*cos(theta)), floor(centre["y"]+r*sin(theta)))))
#         break
#       }
#     }
# 
# 
#     x=c(x,floor(centre["x"]+r*cos(theta)))
#     y=c(y,floor(centre["y"]+r*sin(theta)))
#   }
#   Contour=matrix(data=c(x,y), nrow = length(x), ncol = 2)
# 
#   return(list(coord = Contour, param = c(0.2,30)))
# }


