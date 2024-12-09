library(MASS)
library(readxl)
library(imager)
library(splines)
library(SpatialPack)
library(mmand)
library(e1071)

set.seed(2024)

######################################################################
#
# Quelques fonctions utilisées
#
######################################################################
Affichage_rayon = function(Ref, R, avant){
  plot(Ref, ylim = c(0,1), lwd = 3, type= "l", 
       main = "Comparaison d'un rayon avec celui de référence",
       xlab = "# du point pris du rayon", ylab = "Intensité de gris")
  lines(R, lwd = 2, col = "cyan3")
  lines(avant, col = "orange3")
  legend(x="topright", legend = c("Référence", "Warp", "Avant warp"), 
         fill = c("black", "cyan3", "orange3"), cex = 0.7)
}

Affichage_otolithe = function(image, Contour, SX, SY, noyau, Sweep, Sweep_Warp, angles, param, angle_tranche = 1/4, nb_melange=1){
  # Graphique des courbes x et y du contour séparé avec leurs splines
  plot(angles, Contour[,1], type = "o", xlab = "angle", ylab = "", col = "red", lwd = 2,
       ylim = c(min(c(min(Contour[,1]), min(Contour[,2]))), c(max(c(max(Contour[,1]), max(Contour[,2]))))))
  lines(SX$x,SX$y,lwd=2)
  lines(angles, Contour[,2], type = "o", col = "blue", lwd = 2)
  lines(SY$x,SY$y,lwd=2)
  grid()
  legend(x="bottomleft", legend = c("x","y"), fill = c("red","blue"),cex=0.7)
  
  # Graphique du contour de l'otolithe avec spline
  plot(Contour[,1],Contour[,2],type="o",col="green",lwd=2,xlab="x",ylab="y")
  lines(SX$y,SY$y,lwd=2)
  
  # Graphique de la spline du contour superposé à l'image de l'otolithe
  plot(image)
  lines(SX$y,SY$y,col="red",lwd=1)
  axe_principal = Trouver_axe_principal(noyau,Contour)
  lines(axe_principal[c(1,3)],axe_principal[c(2,4)],col="purple",lwd=2)
  points(noyau[1],noyau[2],col="red")
  
  par(mfrow = c(1,2))
  # Graphique du sweep (avec contraste) avant la transformation
  plot(as.cimg(Sweep),xlim = c(0,length(Sweep[,1])), ylim=c(0,length(Sweep[1,])),axes = F)
  axis(1, at = seq(1,length(Sweep[,1]),len=5), labels = c(paste0("-pi/",1/angle_tranche), paste0("-pi/",2/angle_tranche), 0, paste0("pi/",2/angle_tranche), paste0("pi/",1/angle_tranche)))
  axis(2, at = seq(0,length(Sweep),5))
  title(main = "Sweep de l'otolithe", xlab = "Angle", ylab = "Points pris des rayons")
  
  # Graphique du sweep transformé pour que les anneaux soient alignés
  plot(as.cimg(Sweep_Warp),xlim = c(0,length(Sweep_Warp[,1])), ylim=c(0,length(Sweep_Warp[1,])), axes=F)
  axis(1, at = seq(1,length(Sweep_Warp[,1]),len=5), labels = c(paste0("-pi/",1/angle_tranche), paste0("-pi/",2/angle_tranche), 0, paste0("pi/",2/angle_tranche), paste0("pi/",1/angle_tranche)))
  axis(2, at = seq(0,length(Sweep),5))
  title(main = "Warp du sweep de l'otolithe (init = param du dernier)", xlab = "Angle", ylab = "Points pris des rayons")
  par(mfrow = c(1,1))
  
  # Graphique de la moyenne des intensités de gris du sweep transformé
  # du noyau au contour
  plot(colMeans(Sweep_Warp),type="l",lwd=3, xaxt = "n", ylim=c(0,1),
       xlab = "Noyau à contour", ylab = "Intensité de gris (0 = noir, 1 = blanc)")
  axis(1, at = seq(1,length(colMeans(Sweep_Warp)),len=11), labels = seq(0,1,len=11)) 
  grid()
  
  # Graphique des transformations subies par le sweep
  plot(c(0,ncol(Sweep_Warp)), c(0,ncol(Sweep_Warp)),type="n",xlab = "Coordonnée d'entrée", ylab = "Coordonnée de sortie")
  for(i in 1:nrow(param)){
    if (nb_melange == 1){
      T = Fonction_Warp(param[i,paste0("a",1:nb_melange)], param[i,paste0("a",1:nb_melange)], ncol(Sweep_Warp))
    }else{
      T = Fonction_Warp(param[i,paste0("a",1:nb_melange)], param[i,paste0("a",1:nb_melange)], ncol(Sweep_Warp), c(param[i,paste0("w",1:(nb_melange-1))], 1-sum(param[i,paste0("w",1:(nb_melange-1))])))
    }
    lines(1:ncol(Sweep_Warp),T)
    title(main = "Fonctions Warp")
  }
}

normalisation = function(Sweep_Warp){
  Sweep_Warp_norm = c()
  n = nrow(Sweep_Warp)
  for(i in 1:n){
    JT = lm(y~bs(1:n, df=5), data = data.frame(x=1:n, y=Sweep_Warp[i,]))
    smooth_SW = smooth.spline(1:n, Sweep_Warp[i,], spar = 0.2)$y
    
    JT_R = predict(JT, data.frame(x=1:n))
    
    # Sorte de normalisation des données pour qu'elles oscillent autour
    # d'un point fixe (1) au lieu d'osciller en descendant
    Sweep_Warp_norm = rbind(Sweep_Warp_norm, Sweep_Warp[i,]-JT_R)
  }
  
  plot(as.cimg(Sweep_Warp_norm), ylim = c(0,n))
  
  return(Sweep_Warp_norm)
}

Contraste = function(x,c=1){
  A = log(x/(1-x))
  A = c*(A-mean(A))/sd(A)
  return(1/(1+exp(-A)))
}

Contraste2 = function(image){
  
  f = ecdf(image)
  new_image = matrix(f(image), nrow(image), ncol(image))
  return(new_image)
}

adaptive_equalize <- function(img, window_size = 64) {
  equalized_img <- img  # Initialize equalized image
  
  # Apply adaptive equalization
  for (i in 1:nrow(img)) {
    for (j in 1:ncol(img)) {
      xmin <- max(1, i - window_size)
      xmax <- min(nrow(img), i + window_size)
      ymin <- max(1, j - window_size)
      ymax <- min(ncol(img), j + window_size)
      
      local_ecdf <- ecdf(img[xmin:xmax, ymin:ymax])
      
      equalized_img[i, j] <- local_ecdf(img[i,j])
    }
  }
  
  return(equalized_img)
}

Grain_filter = function(image, a, b=a){
  new_image = image
  for(i in 1:length(image[,1])){
    for(j in 1:length(image[1,])){
      a1 = i-a
      a2 = i+a
      a3 = j-b
      a4 = j+b
      
      if(i-a < 0) a1 = 1
      if(i+a > length(image[,1])) a2 = length(image[,1])
      if(j-b < 0) a3 = 1
      if(j+b > length(image[1,])) a4 = length(image[1,])
      
      new_image[i,j] = mean(image[a1:a2, a3:a4])
    }
  }
  return(new_image)
}


Division_image = function(image_originale){
  image_grise = grayscale(image_originale)
  #tc = optim(c(0.2,30),TrouverFormes,image = as.cimg(image_grise))$par
  image = as.vector(clean(as.cimg(image_grise)>0.2,30))
  image = matrix(image, nrow(image_grise), ncol(image_grise))
  
  div = round(nrow(image)/2)
  
  i=0
  while(sum(image[div+i,]) > 0 && sum(image[div-i,]) > 0){
    i=i+1

    if(sum(image[div+i,]) == 0) div = div + i
    
    if(sum(image[div-i,]) == 0) div = div - i
  }
  
  i = 0
  d=div
  g=div
  
  while(sum(image[div+i,]) == 0 || sum(image[div-i,]) == 0){
    i=i+1

    if(sum(image[div+i,]) == 0) d = div + i
    
    if(sum(image[div-i,]) == 0) g = div - i
  }

  div = round((g+d)/2)
  
  imageG = image_originale[1:div,,,]
  imageD = image_originale[(div+1):nrow(image),,,]
  
  return(list(gauche = imageG, droite = imageD))
}


Compte_anneau = function(Sweep_Warp, angle_tranche = 1/4){
  SW_moy = colMeans(Sweep_Warp)
  
  plot(SW_moy,type = "l")
  
  JT = lm(y~bs(1:length(SW_moy), df=5), data = data.frame(x=1:length(SW_moy), y=SW_moy))
  smooth_SW = smooth.spline(1:length(SW_moy),SW_moy,spar = 0.2)$y
  
  # Superpose SW_moy après être lissé
  lines(smooth_SW, col = "blue")
  
  JT_R = predict(JT, data.frame(x=1:length(SW_moy)))
  
  # Ajoute la spline qui suit la tendance de SW_moy sans les oscillations
  lines(1:length(SW_moy), JT_R, col="red")
  
  # Sorte de normalisation des données pour qu'elles oscillent autour
  # d'un point fixe (1) au lieu d'osciller en descendant
  SW_moy = smooth_SW-JT_R
  
  # Trouve les minimums locaux de SW_moy
  dif = c()
  anneau = c()
  pente = 1
  for(i in 1:(length(SW_moy)-1)){
    dif = c(dif, SW_moy[i+1] - SW_moy[i])
    
    if(dif[i] < 0){
      pente = -1
    }
    else{
      if(pente == -1){
        anneau = c(anneau, i)
      }
      pente = 1
    }
  }
  
  # Graphiques qui montre les anneaux sur la courbe des moyennes d'intensité
  # et sur l'image Sweep_Warp
  plot(seq(0,1,len=length(SW_moy)), SW_moy, type="l", lwd=2)
  abline(v=seq(0,1,len=length(SW_moy))[anneau], col="red")
  grid()
  
  plot(as.cimg(Sweep_Warp),xlim = c(0,length(Sweep_Warp[,1])), ylim=c(0,length(Sweep_Warp[1,])), axes=F)
  axis(1, at = seq(1,length(Sweep_Warp[,1]),len=5), labels = c(paste0("-pi/",1/angle_tranche), paste0("-pi/",2/angle_tranche), 0, paste0("pi/",2/angle_tranche), paste0("pi/",1/angle_tranche)))
  axis(2, at = seq(0,length(Sweep_Warp),5))
  title(main = "Warp du sweep de l'otolithe (init = param du dernier)", xlab = "Angle", ylab = "Points pris des rayons")
  abline(h=anneau, col="red")
  
  # Affiche le nombre d'anneaux trouvé
  library("emoji")
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, paste("Le nombre d'anneaux est \n définitivement, absolument,\n 100% sûr : ",
                               length(anneau), "\n", emoji("grinning face"), emoji("thumbs up")),
       cex = 1.6, col = "black")
  
  return(anneau)
}



######################################################################
#
# Initialisation de variables
#
######################################################################
nb_melange = 3 # Prend plus longtemps plus il y a de mélanges

n_points = 400 # Prend plus longtemps plus il y a de points
n_rayons = 400

grain = 2
######################################################################
#
# Début du main code
#
######################################################################

#image_originale = load.image("C:/Users/THIBODEAUJU/Documents/Otolithe/1994-COMM-40-0001-photo-0001-image-scale-250.jpg")
#image_originale = load.image("C:/Users/THIBODEAUJU/Documents/Otolithe/1993-RES-40-1087-photo-0001-image-scale-150.jpg")
#image_originale = load.image("C:/Users/THIBODEAUJU/Documents/Otolithe/1993-RES-40-1333-photo-0001-image-scale-150.jpg")
#image_originale = load.image("C:/Users/THIBODEAUJU/Documents/Otolithe/1994-COMM-40-0034-photo-0001-image-scale-150.jpg")
#image_originale = load.image("C:/Users/THIBODEAUJU/Documents/Otolithe/1994-RES-40-1268-photo-0001-image-scale-250.jpg")
#image_originale = load.image("C:/Users/THIBODEAUJU/Documents/Otolithe/1993-RES-40-1338-photo-0002-image-scale-150.jpg")
#image_originale = load.image("C:/Users/THIBODEAUJU/Documents/Otolithe/1994-COMM-40-0001-photo-0001-image-scale-250.jpg")
image_originale = load.image("C:/Users/THIBODEAUJU/Documents/Otolithe/plaice.jpg")

# image_divise = Division_image(image_originale)
# rm(image_originale)
# image_originale = as.cimg(image_divise$droite)
# rm(image_divise)

image_grise = grayscale(image_originale) # Transforme l'image en gris
image_grise = as.matrix(image_grise[,,1,1])

# Donne les paramètres de seuil de gris et de nettoyage de l'image
#tc = TrouverContour(image_grise)$param

# Donne les coordonnées du noyau de l'otolithe
noyau = Trouver_noyau(image_grise, 0.2, 30)

#noyau = optim(noyau,Similitude_des_rayons, image = image_grise, Contour = Contour, n_points=n_rayons)$par

# Donne les coordonnées en désordre du contour de l'otolithe
Contour = TrouverContour(image_grise,n_rayons, nb_centre = 5)

pointsContour = Contour$centres
decoupe = Contour$decoupe
Contour = Contour$coord

t = c()
for(i in 1:nrow(Contour)){
  if(Contour[i,2]-noyau[2] < 0){
    t = c(t, 2*pi-acos((Contour[i,1]-noyau[1])/euclidien(noyau, Contour[i,])))
  }else{
    t = c(t, acos((Contour[i,1]-noyau[1])/euclidien(noyau, Contour[i,])))
  }
}
# Replace les coordonnées en ordre
Contour = unique(Contour[order(t), ])

# Spline des x et y du contour de l'otolithe
angles = seq(0,2*pi, length.out = length(Contour[,1])+1)
SX = spline(angles, c(Contour[,1], Contour[1,1]), n=length(Contour[,1]), method = "periodic")
SY = spline(angles, c(Contour[,2], Contour[1,2]), n=length(Contour[,1]), method = "periodic")
angles = seq(0,2*pi, length.out = length(Contour[,1]))

# On prend une tranche de l'otolithe
Tranche = Tranche_Contour(image, noyau, Contour, n_points = n_rayons, angle = pi/4)
  
# On fait le sweep de l'otolithe de 0 à 2pi et on l'étend de sorte
# que le noyau est à y=0 et que le contour soit à 
# y = (# points pris par rayon). Les distances sont normalisées par le
# nombre de points pris par rayon
Sweep = t(Sweep_anneaux(noyau,as.matrix(medianblur(as.cimg(image_grise),10)),Tranche,n_points))
  
# On ajoute plus de contraste
Sweep2 = Contraste(Grain_filter(Sweep, 0))
  
# On trouve l'image qui a été transformé pour faire en sorte que les
# anneaux sont alignés
#k = shapeKernel(c(21,7), type = "box")
Transforme = Transformation(noyau,Tranche,Sweep2, nb_melange = nb_melange)
param=Transforme$param
Sweep_Warp = Transforme$Transforme
  
# On affiche les graphiques
Affichage_otolithe(image_originale, Contour, SX, SY, noyau, Sweep2, Sweep_Warp, angles, param, 1/4, nb_melange)

anneau = Compte_anneau(Sweep_Warp, 1/4)

# Test de la fonction d'inversion du warp
######################################################################
IW1 = matrix(,400,400)
for(i in 1:400){
  indice = (399)*Inverse_Warp((1:400-1)/399, param[i,paste0("a",1:nb_melange)], param[i,paste0("b",1:nb_melange)], c(param[i,paste0("w",1:(nb_melange-1))], 1-sum(param[i,paste0("w",1:(nb_melange-1))])))+1
  IW1[i,] = spline(1:400,Sweep_Warp[i,],xout=indice)$y
}

par(mfrow = c(1,2))
plot(as.cimg(Sweep2),ylim=c(0,400), main = "Sweep avant warp")
plot(as.cimg(IW1), ylim = c(0,400), main = "Ramené à image initiale\n après le warp")
par(mfrow = c(1,1))

IW2 = matrix(,400,length(anneau))
for(i in 1:400){
  IW2[i,] = (399)*Inverse_Warp((anneau-1)/399, param[i,paste0("a",1:nb_melange)], param[i,paste0("b",1:nb_melange)], c(param[i,paste0("w",1:(nb_melange-1))], 1-sum(param[i,paste0("w",1:(nb_melange-1))])))+1
}

plot(as.cimg(Sweep2),ylim=c(0,400))
for(i in 1:ncol(IW2)){
  lines(1:400, IW2[,i], type="l",col="red")
}
#######################################################################


# 
# #Align = function(image, Contour, noyau){
#   nom = c(paste0("x",1:nrow(Sweep2)), paste0("y",1:nrow(Sweep2)))
#   init = c(seq(100,101,len = nrow(Sweep2)), 1:400)
#   init = setNames(init, nom)
#   
#   S = optim(init,Warp_spline,image=Sweep2, control = list("maxit" = 10000))
# #}


Sweep_smooth = Sweep2
for(i in 1:nrow(Sweep2)){
  Sweep_smooth[i,] = smooth.spline(Sweep2[i,])$y
}
# for(j in 1:ncol(Sweep2)){
#   Sweep_smooth[,j] = smooth.spline(Sweep2[,j])$y
# }
plot(as.cimg(matrix(smooth(Sweep2),nrow(Sweep2),ncol(Sweep2))), ylim = c(0,400))


ET=matrix(,ncol(Sweep2),length(anneau))
for(j in 1:length(anneau)){
  i=1
  while(anneau[j]-i > 0 && i+anneau[j] < ncol(Sweep2)){
    ET[i,j] = sd(Sweep_Warp[100,(anneau[j]-i):(anneau[j]+i)])
    i=i+1
  }
}






